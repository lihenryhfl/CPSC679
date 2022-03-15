#include "SETTINGS.h"
#include <cmath>
#include <iostream>
#include <vector>

#if _WIN32
#include <gl/glut.h>
#elif __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "Geometry/TRIANGLE_MESH.h"
#include "Geometry/SQUARE.h"
#include "Hyperelastic/STVK.h"
#include "Hyperelastic/SNH.h"
#include "Timestepper/TIMESTEPPER.h"
//#include "Timestepper/FORWARD_EULER.h"
#include "Timestepper/BACKWARD_EULER.h"
#include "util/FFMPEG_MOVIE.h"

using namespace std;

// the resolution of the OpenGL window -- independent of the field resolution
int xScreenRes = 800;
int yScreenRes = 800;

// Text for the title bar of the window
string windowLabel("Pancake Simulator");

// mouse tracking variables
int xMouse         = -1;
int yMouse         = -1;
int mouseButton    = -1;
int mouseState     = -1;
int mouseModifiers = -1;

// animate the current runEverytime()?
bool animate = false;
bool singleStep = false;

// the current viewer eye position
VECTOR3 eyeCenter(0, -1.25, 1);

// current zoom level into the field
float zoom = 3.75;

//REAL poissonsRatio = 0.0;
REAL poissonsRatio = 0.45;
//REAL youngsModulus = 1.0;
//REAL youngsModulus = 5.0;
REAL youngsModulus = 10.0;
//REAL youngsModulus = 25.0;
//REAL youngsModulus = 50.0;
REAL gravityConstant = -1.0;
REAL eps = 0.02;

TRIANGLE_MESH* triangleMesh = NULL;
TIMESTEPPER* integrator = NULL;
MATERIAL* material = NULL;
VECTOR2 bodyForce;

FFMPEG_MOVIE movie;

// pin one ear
SQUARE ceiling(VECTOR2(0.35,0.95), 1.0);
SQUARE below(VECTOR2(-1,-1), 0.5);

vector<SQUARE> squares;
vector<VECTOR2> centers;

// pin both ears
//SQUARE ceiling(VECTOR2(0,0.95), 1.0);

// pin the tail
//SQUARE ceiling(VECTOR2(0.85,0), 1.0);

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void drawSquare(const SQUARE& square, const VECTOR4& color)
{
  vector<VECTOR2> vertices = square.vertices();

  //glColor4f(0, 0, 1.0, 0.5);
  glColor4f(color[0], color[1], color[2], color[3]);

  glBegin(GL_TRIANGLES);
    glVertex2f(vertices[0][0], vertices[0][1]);
    glVertex2f(vertices[1][0], vertices[1][1]);
    glVertex2f(vertices[2][0], vertices[2][1]);

    glVertex2f(vertices[2][0], vertices[2][1]);
    glVertex2f(vertices[3][0], vertices[3][1]);
    glVertex2f(vertices[0][0], vertices[0][1]);
  glEnd();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void drawTriangle(const int index, const VECTOR4& color)
{
  const vector<VECTOR3I>& triangles = triangleMesh->triangles();
  const vector<VECTOR2>& nodes = triangleMesh->vertices();

  const VECTOR3I& t = triangles[index];
  //glColor4f(254.0 / 255.0, 240.0 / 255, 217.0 / 255.0, 1.0);
  glColor4f(color[0], color[1], color[2], color[3]);
  glBegin(GL_TRIANGLES);
    for (int x = 0; x < 3; x++)
      glVertex2f(nodes[t[x]][0], nodes[t[x]][1]);
  glEnd();

  glColor4f(0,0,0,1);
  glBegin(GL_LINE_LOOP);
    for (int x = 0; x < 3; x++)
      glVertex2f(nodes[t[x]][0], nodes[t[x]][1]);
  glEnd();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void drawMesh(const VECTOR4& color)
{
  for (unsigned int x = 0; x < triangleMesh->triangles().size(); x++)
    drawTriangle(x, color);
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
  gluLookAt(eyeCenter[0], eyeCenter[1], 1,  // eye
            eyeCenter[0], eyeCenter[1], 0,  // center
            0, 1, 0);   // up

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // draw the simulation mesh
  drawMesh(VECTOR4(254.0 / 255.0, 240.0 / 255, 217.0 / 255.0, 1.0));
  //drawMesh(VECTOR4(0.0 / 255.0, 0.0 / 255, 0.0 / 255.0, 0.5));

  // draw the kinematic colliders
  //glColor4f(0, 0, 1.0, 0.5);
  VECTOR4 squareColor(0, 0, 1.0, 0.5);
  drawSquare(ceiling, squareColor);
  //drawSquare(below);
  //
  for (unsigned int x = 0; x < squares.size(); x++) {
    drawSquare(squares[x], squareColor);
  }

  glutSwapBuffers();
}

///////////////////////////////////////////////////////////////////////
// Map the keyboard keys to something here
///////////////////////////////////////////////////////////////////////
void glutKeyboard(unsigned char key, int x, int y)
{
  switch (key)
  {
    case 'Q':
      exit(0);
      break;
    case 'q':
      movie.streamWriteMovie("simulator.mov");
      exit(0);
      break;

    case 'a':
      animate = !animate;
      break;

    case ' ':
      animate = true;
      singleStep = true;
      break;

    case 'v':
      cout << "  eyeCenter  = VECTOR3(" << eyeCenter.transpose() << "); " << endl;
      cout << "  zoom = " << zoom << ";" << endl;
      break;

    case 's':
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

  if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN && modifiers & GLUT_ACTIVE_SHIFT)
  {
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
  if (mouseButton == GLUT_LEFT_BUTTON &&
      mouseState == GLUT_DOWN &&
      mouseModifiers & GLUT_ACTIVE_SHIFT)
  {
    // make sure nothing else is called
    return;
  }

  float xDiff = x - xMouse;
  float yDiff = y - yMouse;
  float speed = 0.001;

  if (mouseButton == GLUT_LEFT_BUTTON)
  {
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
}

///////////////////////////////////////////////////////////////////////
// animate and display new result
///////////////////////////////////////////////////////////////////////
void glutIdle()
{
  if (animate)
  {
    static int frame = 0;
    integrator->clearExternalForces();
    VECTOR2 gravity(0.0, gravityConstant);
    integrator->addGravity(gravity);
    integrator->findCollidedVertices(squares);
    integrator->solve(true);
    frame++;

    if (singleStep)
    {
      animate = false;
      singleStep = false;
    }

    movie.addFrameGL();
  }
  glutPostRedisplay();
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
int glutWindow()
{
  glutInitDisplayMode(GLUT_DOUBLE| GLUT_RGBA);
  glutInitWindowSize(xScreenRes, yScreenRes);
  glutInitWindowPosition(10, 10);
  glutCreateWindow(windowLabel.c_str());

  // set the viewport resolution (w x h)
  glViewport(0, 0, (GLsizei) xScreenRes, (GLsizei) yScreenRes);

  // set the background color to gray
  //glClearColor(0.1, 0.1, 0.1, 0);
  glClearColor(1,1,1,1);

  // register all the callbacks
  glutDisplayFunc(&glutDisplay);
  glutIdleFunc(&glutIdle);
  glutKeyboardFunc(&glutKeyboard);
  glutMouseFunc(&glutMouseClick);
  glutMotionFunc(&glutMouseMotion);
  glutPassiveMotionFunc(&glutPassiveMouseMotion);

  //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  //glEnable(GL_MULTISAMPLE);
  glLineWidth(1.0);
  glEnable(GL_LINE_SMOOTH);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_BLEND);
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

  // enter the infinite GL loop
  glutMainLoop();

  // Control flow will never reach here
  return EXIT_SUCCESS;
}

string toUpper(const string& input)
{
  string copy(input);
  for (unsigned int x = 0; x < input.length(); x++)
    copy[x] = std::toupper(copy[x]);
  return copy;
}

///////////////////////////////////////////////////////////////////////
// Read in a *.node file
///////////////////////////////////////////////////////////////////////
void readNodes(const string& filename, vector<VECTOR2>& nodes, vector<int>& indices)
{
  // read the nodes file
  FILE* file = NULL;
  file = fopen(filename.c_str(), "r");

  if (file == NULL)
  {
    cout << " File " << filename.c_str() << " does not exist! " << endl;
    exit(0);
  }

  // read in the total number of nodes and other attributes
  int totalNodes = -1;
  int dimension = -1;
  int totalAttributes = -1;
  int totalBoundaryMarkers = -1;
  fscanf(file, "%i %i %i %i", &totalNodes, &dimension, &totalAttributes, &totalBoundaryMarkers);

  cout << " Total nodes: " << totalNodes << endl;
  cout << " Dimension: " << dimension << endl;
  cout << " Attributes: " << totalAttributes << endl;
  cout << " Boundary markers: " << totalBoundaryMarkers << endl;

  assert(dimension == 2);

  for (int x = 0; x < totalNodes; x++)
  {
    // get the vertex position
    int index = -1;
    double position[2];
    fscanf(file, "%i %lf %lf", &index, &(position[0]), &(position[1]));
    //cout << " index: " << index << "\t node: " << position[0] << " " << position[1] << endl;

    // store it as a node
    VECTOR2 node(position[0], position[1]);
    nodes.push_back(node);
    indices.push_back(index);

    // parser is pretty brittle, just going to assume indices are sequential
    assert(index == x + 1);

    // strip off the attributes
    double throwAway;
    for (int y = 0; y < totalAttributes; y++)
      fscanf(file, "%lf", &throwAway);

    // strip off the boundary markers
    for (int y = 0; y < totalBoundaryMarkers; y++)
      fscanf(file, "%lf", &throwAway);
  }
  fclose(file);
}

///////////////////////////////////////////////////////////////////////
// Read in a *.ele file
///////////////////////////////////////////////////////////////////////
void readElements(const string& filename, const int offset, vector<VECTOR3I>& triangles)
{
  FILE* file = NULL;
  file = fopen(filename.c_str(), "r");

  if (file == NULL)
  {
    cout << " File " << filename.c_str() << " does not exist! " << endl;
    exit(0);
  }

  int totalTriangles = -1;
  int totalNodesPerTriangle = -1;
  int totalAttributes = -1;
  fscanf(file, "%i %i %i", &totalTriangles, &totalNodesPerTriangle, &totalAttributes);

  cout << " Total triangles: " << totalTriangles << endl;
  cout << " Total nodes in each triangle: " << totalNodesPerTriangle << endl;
  cout << " Total attributes per triangle: " << totalAttributes << endl;
  assert(totalNodesPerTriangle == 3);

  for (int x = 0; x < totalTriangles; x++)
  {
    int triangleIndex = -1;
    int nodeIndices[3];

    fscanf(file, "%i %i %i %i", &triangleIndex, &nodeIndices[0], &nodeIndices[1], &nodeIndices[2]);

    VECTOR3I triangle(nodeIndices[0], nodeIndices[1], nodeIndices[2]);
    triangle -= VECTOR3I(offset, offset, offset);
    triangles.push_back(triangle);

    //cout << " Triangle " << triangleIndex << ": " << triangle[0] << " " << triangle[1] << " " << triangle[2] << endl;
  }
}

///////////////////////////////////////////////////////////////////////
// Load up a 2D triangle mesh
///////////////////////////////////////////////////////////////////////
void loadTriangles2D(const string& prefix,
                     vector<VECTOR2>& nodes,
                     vector<VECTOR3I>& triangles)
{
  // this actually gets thrown away in the end.
  vector<int> indices;

  // read the nodes file
  string nodeFile = prefix + string(".node");
  readNodes(nodeFile, nodes, indices);

  // did the file indices start at 1 or 0? If 1, we need to subtract off one
  // from all the triangles' node indices.
  int offset = (indices[0] == 0) ? 0 : 1;

  // read in the elements file
  string elementFile = prefix + string(".ele");
  readElements(elementFile, offset, triangles);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void readSquare()
{
  vector<VECTOR2> nodes;
  vector<VECTOR3I> triangles;

  nodes.resize(4);
  nodes[0] = VECTOR2(0,0);
  nodes[1] = VECTOR2(1,0);
  nodes[2] = VECTOR2(0,1);
  nodes[3] = VECTOR2(1,1);

  triangles.resize(2);
  triangles[0] = VECTOR3I(0,1,2);
  triangles[1] = VECTOR3I(1,3,2);

  triangleMesh = new TRIANGLE_MESH(nodes, triangles, eps);

  const REAL mu     = MATERIAL::computeMu(youngsModulus, poissonsRatio);
  const REAL lambda = MATERIAL::computeLambda(youngsModulus, poissonsRatio);
  //material = new STVK(mu, lambda);
  material = new SNH(mu, lambda);
  //integrator = new FORWARD_EULER(*triangleMesh, *material);
  integrator = new BACKWARD_EULER(*triangleMesh, *material);
  //integrator->dt() = 1.0 / 30.0;

  //eyeCenter = VECTOR3(0.421, 0.446, 1);
  //zoom = 3.75;
  eyeCenter  = VECTOR3(0.544, -3.182,1);
  zoom = 10.0;

  //eyeCenter  = VECTOR3(0.206, -2.231, 1);
  //zoom = 5.516;

  REAL area = triangleMesh->triangleArea(0) + triangleMesh->triangleArea(1);
  cout << " total area: " << area << endl;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void readSquare2()
{
  vector<VECTOR2> nodes;
  vector<VECTOR3I> triangles;

  REAL xbias = -0.4;
  nodes.resize(4);
  nodes[0] = VECTOR2(0 + xbias,0);
  nodes[1] = VECTOR2(0.2 + xbias,0);
  nodes[2] = VECTOR2(0 + xbias,0.2);
  nodes[3] = VECTOR2(0.2 + xbias,0.2);

  triangles.resize(2);
  triangles[0] = VECTOR3I(0,1,2);
  triangles[1] = VECTOR3I(1,3,2);

  triangleMesh = new TRIANGLE_MESH(nodes, triangles, eps);

  const REAL mu     = MATERIAL::computeMu(youngsModulus, poissonsRatio);
  const REAL lambda = MATERIAL::computeLambda(youngsModulus, poissonsRatio);
  //material = new STVK(mu, lambda);
  material = new SNH(mu, lambda);
  //integrator = new FORWARD_EULER(*triangleMesh, *material);
  integrator = new BACKWARD_EULER(*triangleMesh, *material);
  //integrator->dt() = 1.0 / 30.0;

  squares.reserve(10);
  centers.reserve(10);

  VECTOR2 center(0.35, 0.95);
  SQUARE square(center, 0.75);
  square.rotation() = Eigen::Rotation2D<REAL>(M_PI / 4.0).toRotationMatrix();

  center = VECTOR2(-0.5, -1.0);
  centers.push_back(center);
  square.translation() = centers.back();
  squares.push_back(square);

  center = VECTOR2(0.75, -1.625);
  centers.push_back(center);
  square.translation() = centers.back();
  squares.push_back(square);

  center = VECTOR2(-0.5, -2.25);
  centers.push_back(center);
  square.translation() = centers.back();
  squares.push_back(square);

  center = VECTOR2(0.75, -2.875);
  centers.push_back(center);
  square.translation() = centers.back();
  squares.push_back(square);

  center = VECTOR2(-0.5, -3.5);
  centers.push_back(center);
  square.translation() = centers.back();
  squares.push_back(square);

  center = VECTOR2(0.75, -4.125);
  centers.push_back(center);
  square.translation() = centers.back();
  squares.push_back(square);

  center = VECTOR2(0.0, -7.125);
  centers.push_back(center);
  square.translation() = centers.back();
  square.scale() = MATRIX2::Identity() * 6.0;
  square.scaleInverse() = square.scale().inverse();
  square.rotation() = Eigen::Rotation2D<REAL>(0.0).toRotationMatrix();
  squares.push_back(square);
  //eyeCenter = VECTOR3(0.421, 0.446, 1);
  //zoom = 3.75;
  eyeCenter  = VECTOR3(0.544, -3.182,1);
  zoom = 10.0;

  //eyeCenter  = VECTOR3(0.206, -2.231, 1);
  //zoom = 5.516;

  REAL area = triangleMesh->triangleArea(0) + triangleMesh->triangleArea(1);
  cout << " total area: " << area << endl;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void readTriangle()
{
  vector<VECTOR2> nodes;
  vector<VECTOR3I> triangles;

  nodes.resize(3);
  nodes[0] = VECTOR2(0,0);
  nodes[1] = VECTOR2(1,0);
  nodes[2] = VECTOR2(0,1);

  triangles.resize(1);
  triangles[0] = VECTOR3I(0,1,2);

  triangleMesh = new TRIANGLE_MESH(nodes, triangles, eps);
  const REAL mu     = MATERIAL::computeMu(youngsModulus, poissonsRatio);
  const REAL lambda = MATERIAL::computeLambda(youngsModulus, poissonsRatio);
  //material = new STVK(mu, lambda);
  material = new SNH(mu, lambda);
  //integrator = new FORWARD_EULER(*triangleMesh, *material);
  integrator = new BACKWARD_EULER(*triangleMesh, *material);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void readBunny()
{
  vector<VECTOR2> nodes;
  vector<VECTOR3I> triangles;

  loadTriangles2D("./data/bunny/bunny.1", nodes, triangles);
  triangleMesh = new TRIANGLE_MESH(nodes, triangles, eps);
  const REAL mu     = MATERIAL::computeMu(youngsModulus, poissonsRatio);
  const REAL lambda = MATERIAL::computeLambda(youngsModulus, poissonsRatio);
  cout << " Mu: " << mu << " Lambda: " << lambda << endl;
  //material = new STVK(mu, lambda);
  material = new SNH(mu, lambda);
  //integrator = new FORWARD_EULER(*triangleMesh, *material);
  integrator = new BACKWARD_EULER(*triangleMesh, *material);
  //integrator->dt() = 1.0 / 30.0;

  // pin down the vertices inside the ceiling
  //integrator->applyPinConstraints(ceiling);
  below.rotation() = Eigen::Rotation2D<REAL>(M_PI / 4.0).toRotationMatrix();

  squares.reserve(10);
  centers.reserve(10);

  VECTOR2 center(0.35, 0.95);
  SQUARE square(center, 0.75);
  square.rotation() = Eigen::Rotation2D<REAL>(M_PI / 4.0).toRotationMatrix();

  center = VECTOR2(-0.5, -1.0);
  centers.push_back(center);
  square.translation() = centers.back();
  squares.push_back(square);

  center = VECTOR2(0.75, -1.625);
  centers.push_back(center);
  square.translation() = centers.back();
  squares.push_back(square);

  center = VECTOR2(-0.5, -2.25);
  centers.push_back(center);
  square.translation() = centers.back();
  squares.push_back(square);

  center = VECTOR2(0.75, -2.875);
  centers.push_back(center);
  square.translation() = centers.back();
  squares.push_back(square);

  center = VECTOR2(-0.5, -3.5);
  centers.push_back(center);
  square.translation() = centers.back();
  squares.push_back(square);

  center = VECTOR2(0.75, -4.125);
  centers.push_back(center);
  square.translation() = centers.back();
  squares.push_back(square);

  center = VECTOR2(0.0, -7.125);
  centers.push_back(center);
  square.translation() = centers.back();
  square.scale() = MATRIX2::Identity() * 6.0;
  square.scaleInverse() = square.scale().inverse();
  square.rotation() = Eigen::Rotation2D<REAL>(0.0).toRotationMatrix();
  squares.push_back(square);

  eyeCenter  = VECTOR3(0.206, -2.231, 1);
  zoom = 5.516;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void readBunny2()
{
  vector<VECTOR2> nodes;
  vector<VECTOR3I> triangles;

  loadTriangles2D("./data/bunny/bunny.1", nodes, triangles);
  triangleMesh = new TRIANGLE_MESH(nodes, triangles, eps);
  const REAL mu     = MATERIAL::computeMu(youngsModulus, poissonsRatio);
  const REAL lambda = MATERIAL::computeLambda(youngsModulus, poissonsRatio);
  cout << " Mu: " << mu << " Lambda: " << lambda << endl;
  //material = new STVK(mu, lambda);
  material = new SNH(mu, lambda);
  //integrator = new FORWARD_EULER(*triangleMesh, *material);
  integrator = new BACKWARD_EULER(*triangleMesh, *material);
  //integrator->dt() = 1.0 / 30.0;

  // pin down the vertices inside the ceiling
  integrator->applyPinConstraints(ceiling);
  below.rotation() = Eigen::Rotation2D<REAL>(M_PI / 4.0).toRotationMatrix();

  squares.reserve(10);
  centers.reserve(10);

  REAL leftcenter = -0.6;
  REAL rightcenter = 1.1;

  // this first square is never actually used
  VECTOR2 center(0.35, 0.95);
  SQUARE square(center, 1.0);
  square.rotation() = Eigen::Rotation2D<REAL>(M_PI / 4.0).toRotationMatrix();

  center = VECTOR2(leftcenter, -1.0);
  centers.push_back(center);
  square.translation() = centers.back();
  squares.push_back(square);

  center = VECTOR2(rightcenter, -2.0);
  centers.push_back(center);
  square.translation() = centers.back();
  squares.push_back(square);

  center = VECTOR2(leftcenter, -3.0);
  centers.push_back(center);
  square.translation() = centers.back();
  squares.push_back(square);

  center = VECTOR2(rightcenter, -4.125);
  centers.push_back(center);
  square.translation() = centers.back();
  squares.push_back(square);

  center = VECTOR2(0.0, -7.125);
  centers.push_back(center);
  square.translation() = centers.back();
  square.scale() = MATRIX2::Identity() * 6.0;
  square.scaleInverse() = square.scale().inverse();
  square.rotation() = Eigen::Rotation2D<REAL>(0.0).toRotationMatrix();
  squares.push_back(square);

  eyeCenter  = VECTOR3(0.206, -2.231, 1);
  zoom = 5.516;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void readCustom()
{
  vector<VECTOR2> nodes;
  vector<VECTOR3I> triangles;

  loadTriangles2D("./data/custom/custom.1", nodes, triangles);
  triangleMesh = new TRIANGLE_MESH(nodes, triangles, eps);
  const REAL mu     = MATERIAL::computeMu(youngsModulus, poissonsRatio);
  const REAL lambda = MATERIAL::computeLambda(youngsModulus, poissonsRatio);
  cout << " Mu: " << mu << " Lambda: " << lambda << endl;
  //material = new STVK(mu, lambda);
  material = new SNH(mu, lambda);
  //integrator = new FORWARD_EULER(*triangleMesh, *material);
  integrator = new BACKWARD_EULER(*triangleMesh, *material);
  //integrator->dt() = 1.0 / 30.0;

  // pin down the vertices inside the ceiling
  //integrator->applyPinConstraints(ceiling);
  below.rotation() = Eigen::Rotation2D<REAL>(M_PI / 4.0).toRotationMatrix();

  squares.reserve(10);
  centers.reserve(10);

  VECTOR2 center(0.35, 0.95);
  SQUARE square(center, 0.75);
  square.rotation() = Eigen::Rotation2D<REAL>(M_PI / 4.0).toRotationMatrix();

  center = VECTOR2(-0.5, -1.0);
  centers.push_back(center);
  square.translation() = centers.back();
  squares.push_back(square);

  center = VECTOR2(0.75, -1.625);
  centers.push_back(center);
  square.translation() = centers.back();
  squares.push_back(square);

  center = VECTOR2(-0.5, -2.25);
  centers.push_back(center);
  square.translation() = centers.back();
  squares.push_back(square);

  center = VECTOR2(0.75, -2.875);
  centers.push_back(center);
  square.translation() = centers.back();
  squares.push_back(square);

  center = VECTOR2(-0.5, -3.5);
  centers.push_back(center);
  square.translation() = centers.back();
  squares.push_back(square);

  center = VECTOR2(0.75, -4.125);
  centers.push_back(center);
  square.translation() = centers.back();
  squares.push_back(square);

  center = VECTOR2(0.0, -7.125);
  centers.push_back(center);
  square.translation() = centers.back();
  square.scale() = MATRIX2::Identity() * 6.0;
  square.scaleInverse() = square.scale().inverse();
  square.rotation() = Eigen::Rotation2D<REAL>(0.0).toRotationMatrix();
  squares.push_back(square);

  eyeCenter  = VECTOR3(0.206, -2.231, 1);
  zoom = 5.516;
}


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  cout << " Usage: " << argv[0] << endl;

  if (argc > 1 && !strcmp(argv[1], "CUSTOM")) {
    cout << "LOADING CUSTOM" << endl;
    readCustom();
  }
  else {
    cout << "LOADING BUNNY" << endl;
    readBunny();
    //readBunny2();
    //readTriangle();
    //readSquare();
    //readSquare2();
  }

  // initialize GLUT and GL
  glutInit(&argc, argv);

  // open the GL window
  glutWindow();
  return 1;
}
