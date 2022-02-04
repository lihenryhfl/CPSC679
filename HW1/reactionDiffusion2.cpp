#include "FIELD_2D.h"
#include <cmath>
#include <iostream>

#if _WIN32
#include <gl/glut.h>
#elif __APPLE__
#include <GLUT/glut.h>
#elif __linux__
#include <GL/glut.h>
#endif

using namespace std;

// using FHN? (or else GS)
bool fhn;

// resolution of the field
int xRes = 200;
int yRes = 200;

// the field being drawn and manipulated
FIELD_2D field(xRes, yRes);
FIELD_2D field_b(xRes, yRes);

// the resolution of the OpenGL window -- independent of the field resolution
int xScreenRes = 850;
int yScreenRes = 850;

// Text for the title bar of the window
string windowLabel("Reaction Diffusion");

// mouse tracking variables
int xMouse = -1;
int yMouse = -1;
int mouseButton = -1;
int mouseState = -1;
int mouseModifiers = -1;

// current grid cell the mouse is pointing at
int xField = -1;
int yField = -1;

// animate the current runEverytime()?
bool animate = true;

// draw the grid over the field?
bool drawingGrid = true;

// print out what the mouse is pointing at?
bool drawingValues = true;

// currently capturing frames for a movie?
bool captureMovie = false;

// the current viewer eye position
float eyeCenter[] = {0.5, 0.5, 1};

// current zoom level into the field
float zoom = 1.0;

// forward declare the caching function here so that we can
// put it at the bottom of the file
void runOnce();

// forward declare the timestepping function here so that we can
// put it at the bottom of the file
void runEverytime();

// forward declare the method-specific timestepping function
//void runForGS(float d_a = 2e-4, float d_b = 1e-5, float dt = 0.1, float f = 0.05, float k = 0.0675);
void runForGS(float d_a = 2e-5, float d_b = 1e-5, float dt = 0.61, float f = 0.04, float k = 0.06);
void runForFHN(float d_a = 0.75, float d_b = 0.0, float dt = 0.02, float alpha = 0.75, float beta = 0.01, float epsilon = 0.02);

///////////////////////////////////////////////////////////////////////
// Figure out which field element is being pointed at, set xField and
// yField to them
///////////////////////////////////////////////////////////////////////
void refreshMouseFieldIndex(int x, int y)
{
  // make the lower left the origin
  y = yScreenRes - y;

  float xNorm = (float)x / xScreenRes;
  float yNorm = (float)y / yScreenRes;

  float halfZoom = 0.5 * zoom;
  float xWorldMin = eyeCenter[0] - halfZoom;
  float xWorldMax = eyeCenter[0] + halfZoom;

  // get the bounds of the field in screen coordinates
  //
  // if non-square textures are ever supported, change the 0.0 and 1.0 below
  float xMin = (0.0 - xWorldMin) / (xWorldMax - xWorldMin);
  float xMax = (1.0 - xWorldMin) / (xWorldMax - xWorldMin);

  float yWorldMin = eyeCenter[1] - halfZoom;
  float yWorldMax = eyeCenter[1] + halfZoom;

  float yMin = (0.0 - yWorldMin) / (yWorldMax - yWorldMin);
  float yMax = (1.0 - yWorldMin) / (yWorldMax - yWorldMin);

  float xScale = 1.0;
  float yScale = 1.0;

  if (xRes < yRes)
    xScale = (float)yRes / xRes;
  if (xRes > yRes)
    yScale = (float)xRes / yRes;

  // index into the field after normalizing according to screen
  // coordinates
  xField = xScale * xRes * ((xNorm - xMin) / (xMax - xMin));
  yField = yScale * yRes * ((yNorm - yMin) / (yMax - yMin));

  // clamp to something inside the field
  xField = (xField < 0) ? 0 : xField;
  xField = (xField >= xRes) ? xRes - 1 : xField;
  yField = (yField < 0) ? 0 : yField;
  yField = (yField >= yRes) ? yRes - 1 : yField;
}

///////////////////////////////////////////////////////////////////////
// Print a string to the GL window
///////////////////////////////////////////////////////////////////////
void printGlString(string output)
{
  glColor4f(1.0f, 0.0f, 0.0f, 1.0f);
  for (unsigned int x = 0; x < output.size(); x++)
    glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, output[x]);
}

///////////////////////////////////////////////////////////////////////
// dump the field contents to a GL texture for drawing
///////////////////////////////////////////////////////////////////////
void updateTexture(FIELD_2D &texture)
{
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glTexImage2D(GL_TEXTURE_2D, 0, 3, texture.xRes(), texture.yRes(), 0, GL_LUMINANCE, GL_FLOAT, texture.data());

  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
  glEnable(GL_TEXTURE_2D);
}

///////////////////////////////////////////////////////////////////////
// draw a grid over everything
///////////////////////////////////////////////////////////////////////
void drawGrid()
{
  glColor4f(0.1, 0.1, 0.1, 1.0);

  float dx = 1.0 / xRes;
  float dy = 1.0 / yRes;

  if (xRes < yRes)
    dx *= (float)xRes / yRes;
  if (xRes > yRes)
    dy *= (float)yRes / xRes;

  glBegin(GL_LINES);
  for (int x = 0; x < field.xRes() + 1; x++) {
    glVertex3f(x * dx, 0, 1);
    glVertex3f(x * dx, 1, 1);
  }
  for (int y = 0; y < field.yRes() + 1; y++) {
    glVertex3f(0, y * dy, 1);
    glVertex3f(1, y * dy, 1);
  }
  glEnd();
}

///////////////////////////////////////////////////////////////////////
// GL and GLUT callbacks
///////////////////////////////////////////////////////////////////////
void glutDisplay()
{
  // Make ensuing transforms affect the projection matrix
  glMatrixMode(GL_PROJECTION);

  // set the projection matrix to an orthographic view
  glLoadIdentity();
  float halfZoom = zoom * 0.5;

  glOrtho(-halfZoom, halfZoom, -halfZoom, halfZoom, -10, 10);

  // set the matrix mode back to modelview
  glMatrixMode(GL_MODELVIEW);

  // set the lookat transform
  glLoadIdentity();
  gluLookAt(eyeCenter[0], eyeCenter[1], 1, // eye
            eyeCenter[0], eyeCenter[1], 0, // center
            0, 1, 0);                      // up

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  float xLength = 1.0;
  float yLength = 1.0;

  if (xRes < yRes)
    xLength = (float)xRes / yRes;
  if (yRes < xRes)
    yLength = (float)yRes / xRes;

  glEnable(GL_TEXTURE_2D);
  glBegin(GL_QUADS);
  glTexCoord2f(0.0, 0.0);
  glVertex3f(0.0, 0.0, 0.0);
  glTexCoord2f(0.0, 1.0);
  glVertex3f(0.0, yLength, 0.0);
  glTexCoord2f(1.0, 1.0);
  glVertex3f(xLength, yLength, 0.0);
  glTexCoord2f(1.0, 0.0);
  glVertex3f(xLength, 0.0, 0.0);
  glEnd();
  glDisable(GL_TEXTURE_2D);

  // draw the grid, but only if the user wants
  if (drawingGrid)
    drawGrid();

  // if there's a valid field index, print it
  if (xField >= 0 && yField >= 0 && xField < field.xRes() && yField < field.yRes()) {
    glLoadIdentity();

    // must set color before setting raster position, otherwise it won't take
    glColor4f(1.0f, 0.0f, 0.0f, 1.0f);

    // normalized screen coordinates (-0.5, 0.5), due to the glLoadIdentity
    float halfZoom = 0.5 * zoom;
    glRasterPos3f(-halfZoom * 0.95, -halfZoom * 0.95, 0);

    // build the field value string
    char buffer[256];
    string fieldValue("(");
    sprintf(buffer, "%i", xField);
    fieldValue = fieldValue + string(buffer);
    sprintf(buffer, "%i", yField);
    fieldValue = fieldValue + string(", ") + string(buffer) + string(") = ");
    sprintf(buffer, "%f", field(xField, yField));
    fieldValue = fieldValue + string(buffer);

    // draw the grid, but only if the user wants
    if (drawingValues)
      printGlString(fieldValue);
  }

  glutSwapBuffers();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void printCommands()
{
  cout << "=============================================================== " << endl;
  cout << " Field viewer code CPSC 679" << endl;
  cout << "=============================================================== " << endl;
  cout << " q           - quit" << endl;
  cout << " v           - type the value of the cell under the mouse" << endl;
  cout << " g           - throw a grid over everything" << endl;
  cout << " w           - write out a PPM file " << endl;
  cout << " left mouse  - pan around" << endl;
  cout << " right mouse - zoom in and out " << endl;
  cout << " shift left mouse - draw on the grid " << endl;
}

///////////////////////////////////////////////////////////////////////
// Map the arrow keys to something here
///////////////////////////////////////////////////////////////////////
void glutSpecial(int key, int x, int y)
{
  switch (key) {
  case GLUT_KEY_LEFT:
    break;
  case GLUT_KEY_RIGHT:
    break;
  case GLUT_KEY_UP:
    break;
  case GLUT_KEY_DOWN:
    break;
  default:
    break;
  }
}

///////////////////////////////////////////////////////////////////////
// Map the keyboard keys to something here
///////////////////////////////////////////////////////////////////////
void glutKeyboard(unsigned char key, int x, int y)
{
  switch (key) {
  case 'a':
    animate = !animate;
    break;
  case 'g':
    drawingGrid = !drawingGrid;
    break;
  case '?':
    printCommands();
    break;
  case 'v':
    drawingValues = !drawingValues;
    break;
  case 'w': {
    static int count = 0;
    char buffer[256];
    sprintf(buffer, "output_%i.ppm", count);
    field.writePPM(buffer);
    count++;
  } break;
  case 'q':
    exit(0);
    break;
  default:
    break;
  }
}

///////////////////////////////////////////////////////////////////////
// Do something if the mouse is clicked
///////////////////////////////////////////////////////////////////////
void glutMouseClick(int button, int state, int x, int y)
{
  int modifiers = glutGetModifiers();
  mouseButton = button;
  mouseState = state;
  mouseModifiers = modifiers;

  if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN && modifiers & GLUT_ACTIVE_SHIFT) {
    // figure out which cell we're pointing at
    refreshMouseFieldIndex(x, y);

    // set the cell
    field(xField, yField) = 1;

    // make sure nothing else is called
    return;
  }

  xMouse = x;
  yMouse = y;
}

///////////////////////////////////////////////////////////////////////
// Do something if the mouse is clicked and moving
///////////////////////////////////////////////////////////////////////
void glutMouseMotion(int x, int y)
{
  if (mouseButton == GLUT_LEFT_BUTTON && mouseState == GLUT_DOWN && mouseModifiers & GLUT_ACTIVE_SHIFT) {
    // figure out which cell we're pointing at
    refreshMouseFieldIndex(x, y);

    // set the cell
    field(xField, yField) = 1;

    // make sure nothing else is called
    return;
  }

  float xDiff = x - xMouse;
  float yDiff = y - yMouse;
  float speed = 0.001;

  if (mouseButton == GLUT_LEFT_BUTTON) {
    eyeCenter[0] -= xDiff * speed;
    eyeCenter[1] += yDiff * speed;
  }
  if (mouseButton == GLUT_RIGHT_BUTTON)
    zoom -= yDiff * speed;

  xMouse = x;
  yMouse = y;
}

///////////////////////////////////////////////////////////////////////
// Do something if the mouse is not clicked and moving
///////////////////////////////////////////////////////////////////////
void glutPassiveMouseMotion(int x, int y)
{
  refreshMouseFieldIndex(x, y);
}

///////////////////////////////////////////////////////////////////////
// animate and display new result
///////////////////////////////////////////////////////////////////////
void glutIdle()
{
  if (animate) {
    runEverytime();
  }
  updateTexture(field);
  glutPostRedisplay();
}

//////////////////////////////////////////////////////////////////////////////
// open the GLVU window
//////////////////////////////////////////////////////////////////////////////
int glvuWindow()
{
  glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
  glutInitWindowSize(xScreenRes, yScreenRes);
  glutInitWindowPosition(10, 10);
  glutCreateWindow(windowLabel.c_str());

  // set the viewport resolution (w x h)
  glViewport(0, 0, (GLsizei)xScreenRes, (GLsizei)yScreenRes);

  // set the background color to gray
  glClearColor(0.1, 0.1, 0.1, 0);

  // register all the callbacks
  glutDisplayFunc(&glutDisplay);
  glutIdleFunc(&glutIdle);
  glutKeyboardFunc(&glutKeyboard);
  glutSpecialFunc(&glutSpecial);
  glutMouseFunc(&glutMouseClick);
  glutMotionFunc(&glutMouseMotion);
  glutPassiveMotionFunc(&glutPassiveMouseMotion);

  // enter the infinite GL loop
  glutMainLoop();

  // Control flow will never reach here
  return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  // In case the field is rectangular, make sure to center the eye
  if (xRes < yRes) {
    float xLength = (float)xRes / yRes;
    eyeCenter[0] = xLength * 0.5;
  }
  if (yRes < xRes) {
    float yLength = (float)yRes / xRes;
    eyeCenter[1] = yLength * 0.5;
  }

  printCommands();
  const string GS("GS");
  const string FHN("FHN");
  if (argc >= 2)
  {
    string argv1(argv[1]);

    if (argv1.compare(GS) == 0) {
      cout << " Simulating Gray-Scott reaction-diffusion " << endl;
      fhn = false;
    } else if (argv1.compare(FHN) == 0) {
      cout << " Simulating FitzHugh-Nagumo reaction-diffusion " << endl;
      fhn = true;
    }
  }

  runOnce();

  // initialize GLUT and GL
  glutInit(&argc, argv);

  // open the GL window
  glvuWindow();
  return 1;
}


///////////////////////////////////////////////////////////////////////
// These functions contain FHN / GS specific logic
///////////////////////////////////////////////////////////////////////

bool fieldUnstable(FIELD_2D &field, string printstr="") {
  float max_val = field.max();
  bool nan_bool = field.isnan();
  if (max_val > 1e3) {
    cout << "Over stability threshold for field " << printstr << " ! max_val = " << max_val << endl;
    return true;
  } else if (nan_bool) {
    cout << "NaNs in field " << printstr <<"!" << endl;
    return true;
  } else {
    return false;
  }
}

FIELD_2D &computeLaplacian(const FIELD_2D &field_, FIELD_2D &laplacian)
{
  // compute convolution
  for (int y = 1; y < yRes - 1; y++) {
    for (int x = 1; x < xRes - 1; x++) {
      laplacian(x,y) = -4 * field_(x,y);
      laplacian(x,y) += field_(x-1,y);
      laplacian(x,y) += field_(x+1,y);
      laplacian(x,y) += field_(x,y-1);
      laplacian(x,y) += field_(x,y+1);
    }
  }

  return laplacian;
}

void zeroBoundary(FIELD_2D &field_)
{
  // zero xth column
  for (int y = 0; y < yRes; y++) {
      field_(0,y) = 0.0;
      field_(xRes - 1,y) = 0.0;
    }
  for (int x = 0; x < xRes; x++) {
      field_(x,0) = 0.0;
      field_(x,yRes - 1) = 0.0;
  }
}
void runForGS(float d_a, float d_b, float dt, float f, float k)
{
  // initialize time step field
  FIELD_2D field_step(xRes, yRes);
  FIELD_2D field_step_b(xRes, yRes);
  FIELD_2D field_tmp(xRes, yRes);
  FIELD_2D laplacian(xRes, yRes);
  float dx;

  // set d_a and d_b properly
  //dx = 4.0 / xRes;
  dx = 2.0 / xRes;
  d_a = d_a / dx / dx;
  d_b = d_b / dx / dx;

  //// add initial condition
  //for (int y = 0.45 * yRes; y < 0.55 * yRes; y++) {
    //for (int x = 0.45 * xRes; x < 0.55 * xRes; x++) {
      //field_b(x,y) += 1.0;
    //}
  //}

  // react chemicals
  field_tmp = field;
  field_tmp *= field_b;
  field_tmp *= field_b;

  fieldUnstable(field_tmp, "field_tmp");

  field_step = 1.0;
  field_step -= field;
  field_step *= f;
  field_step -= field_tmp;

  fieldUnstable(field_step, "field_step");

  field_step_b = field_b;
  field_step_b *= -(f + k);
  field_step_b += field_tmp;

  fieldUnstable(field_step_b, "field_step_b");

  // zero the boundaries
  zeroBoundary(field);
  zeroBoundary(field_b);

  // compute laplacian
  computeLaplacian(field, laplacian);
  //string s = string("./outputs/GS_laplacian_a.MAT");
  //laplacian.writeMatlab(s, string("x"));
  field_step += d_a * laplacian;
  field += dt * field_step;

  computeLaplacian(field_b, laplacian);
  //s = string("./outputs/GS_laplacian_b.MAT");
  //laplacian.writeMatlab(s, string("x"));
  field_step_b += d_b * laplacian;
  field_b += dt * field_step_b;

  fieldUnstable(field_b, "field_b");
  fieldUnstable(field, "field");

  // zero the boundaries
  zeroBoundary(field);
  zeroBoundary(field_b);
}

void runForFHN(float d_a, float d_b, float dt, float alpha, float beta, float epsilon)
{
  // initialize time step field
  FIELD_2D field_step(xRes, yRes);

  // add initial condition
  for (int y = 0.45 * yRes; y < 0.55 * yRes; y++) {
    for (int x = 0.45 * xRes; x < 0.55 * xRes; x++) {
      field(x,y) += 1.0;
      field_b(x,y) += alpha / 2;
    }
  }

  // react chemicals


  // zero the boundaries
  zeroBoundary(field);
  zeroBoundary(field_b);
}

///////////////////////////////////////////////////////////////////////
// This function is called every frame -- do something interesting
// here.
///////////////////////////////////////////////////////////////////////
void runEverytime()
{
  static int counter = 0;

  if (counter == 0) {
    cout << "Is it unstable to begin with?" << endl;
    fieldUnstable(field, "field");
    fieldUnstable(field_b, "field_b");
    cout << "Now starting simulation." << endl;
  }

  //// add initial condition
  //for (int y = 0.45 * yRes; y < 0.55 * yRes; y++) {
    //for (int x = 0.45 * xRes; x < 0.55 * xRes; x++) {
      ////field_b(x,y) += 1.0;
      //field(x,y) += 1.0;
    //}
  //}

  if (!fieldUnstable(field, "if field") && !fieldUnstable(field_b, "if field_b")) {
  //if (counter < 1 && !fieldUnstable(field, "if field") && !fieldUnstable(field_b, "if field_b")) {
    //string s = string("./outputs/GS_a_") + to_string(counter) + string(".MAT");
    //field.writeMatlab(s, string("x"));
    //s = string("./outputs/GS_b_") + to_string(counter) + string(".MAT");
    //field_b.writeMatlab(s, string("x"));
    if (fhn) {
      for (int iters = 0; iters < 10; iters++)
        runForFHN();
    } else {
      for (int iters = 0; iters < 100; iters++)
      //for (int iters = 0; iters < 1; iters++)
        runForGS();
    }
  }

    //field.normalize();
    //field_b.normalize();

  counter++;

}

///////////////////////////////////////////////////////////////////////
// This is called once at the beginning so you can precache
// something here
///////////////////////////////////////////////////////////////////////
void runOnce()
{
  float x_, y_;
  float fx, fy;
  //// let's insert gray everywhere
  //for (int y = 0; y < yRes; y++)
    //for (int x = 0; x < xRes; x++)
      //field(x,y) = 0.25;

  //// let's insert gray everywhere
  //for (int y = 0; y < yRes; y++)
    //for (int x = 0; x < xRes; x++)
      //field_b(x,y) = 0.25;

  for (int y = 0; y < yRes; y++) {
    for (int x = 0; x < xRes; x++) {
      fx = (float) x;
      fy = (float) y;
      x_ = fx / xRes * 2 - 1.0;
      y_ = fy / yRes * 2 - 1.0;
      field(x,y) = 1 - exp(-80 * (pow(x_ + 0.05, 2) + pow(y_ + 0.05, 2)));
    }
  }

  for (int y = 0; y < yRes; y++) {
    for (int x = 0; x < xRes; x++) {
      fx = (float) x;
      fy = (float) y;
      x_ = fx / xRes * 2 - 1.0;
      y_ = fy / yRes * 2 - 1.0;
      field_b(x,y) = exp(-80 * (pow(x_ - 0.05, 2) + pow(y_ - 0.05, 2)));
      //if (x > 80 && x < 100 && y > 80 && y < 100)
        //cout << "x: " << x << ", y:" << y << ", x_: " << x_ << ", y_: " << y_ << ", field(x,y): " <<  field(x,y) << endl;
    }
  }

  // zero the boundaries
  zeroBoundary(field);
  zeroBoundary(field_b);
  //field.normalize();
  //field_b.normalize();
}
