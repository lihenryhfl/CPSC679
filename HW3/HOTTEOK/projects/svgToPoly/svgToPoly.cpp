#include "SETTINGS.h"
#include <cmath>
#include <iostream>

#define NANOSVG_IMPLEMENTATION
#include "nanosvg.h"

#include <vector>

using namespace std;

// the actual triangle data
vector<VECTOR2> nodes;
vector<int> indices;
vector<VECTOR3I> triangles;

VECTOR2 mins, maxs;

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
    cout << " index: " << index << "\t node: " << position[0] << " " << position[1] << endl;

    // store it as a node
    VECTOR2 node(position[0], position[1]);
    nodes.push_back(node);
    indices.push_back(index);

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

    cout << " Triangle " << triangleIndex << ": " << triangle[0] << " " << triangle[1] << " " << triangle[2] << endl;
  }
}

///////////////////////////////////////////////////////////////////////
// Load up a 2D triangle mesh
///////////////////////////////////////////////////////////////////////
void loadTriangles2D(const string& prefix)
{
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
// Get the bounding box of the nodes
///////////////////////////////////////////////////////////////////////
void getBoundingBox(VECTOR2& localMins, VECTOR2& maxs)
{
  mins = nodes[0];
  maxs = nodes[0];

  for (unsigned int x = 1; x < nodes.size(); x++)
  {
    if (nodes[x][0] < mins[0])
      mins[0] = nodes[x][0];
    if (nodes[x][1] < mins[1])
      mins[1] = nodes[x][1];
    if (nodes[x][0] > maxs[0])
      maxs[0] = nodes[x][0];
    if (nodes[x][1] > maxs[1])
      maxs[1] = nodes[x][1];
  }

  cout << " Bounding box mins: " << mins[0] << " " << mins[1] << endl;
  cout << " Bounding box maxs: " << maxs[0] << " " << maxs[1] << endl;
}

///////////////////////////////////////////////////////////////////////
// write the SVG path out the a POLY file for Triangle
///////////////////////////////////////////////////////////////////////
void writePoly(const char* filename)
{
  FILE* file = NULL;
  file = fopen(filename, "w");
  if (file == NULL)
  {
    cout << " Could not open file " << filename << "!!!" << endl;
    exit(0);
  }
  cout << " Writing out " << filename << endl;

  // output vertices
  fprintf(file, "%i 2 0 0\n", (int)nodes.size());
  for (unsigned int x = 0; x < nodes.size(); x++)
    fprintf(file, "%i %f %f\n", x + 1, (float)nodes[x][0], (float)nodes[x][1]);

  // output edges
  fprintf(file, "%i 0\n", (int)nodes.size());
  for (unsigned int x = 0; x < nodes.size() - 1; x++)
    fprintf(file, "%i %i %i\n", x + 1, x + 1, x + 2);
  fprintf(file, "%i %i 1\n", (int)nodes.size(), (int)nodes.size());

  // output holes
  fprintf(file, "0\n");

  fclose(file);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  if (argc < 3)
  {
    cout << " USAGE: " << argv[0] << " <SVG filename> <POLY filename>" << endl;
    return 0;
  }

	NSVGimage* image;
	image = nsvgParseFromFile(argv[1], "px", 96);
	printf("size: %f x %f\n", image->width, image->height);

	NSVGshape* shape;
	NSVGpath* path;
	for (shape = image->shapes; shape != NULL; shape = shape->next) {
		for (path = shape->paths; path != NULL; path = path->next) {
			for (int i = 0; i < path->npts-1; i += 3) {
				float* p = &path->pts[i*2];
        cout << " point " << i << ": " << p[0] << " " << p[1] << endl;

        VECTOR2 node(p[0], p[1]);
        nodes.push_back(node);
      }
    }
  }

	// Delete
	nsvgDelete(image);
  getBoundingBox(mins, maxs);

  VECTOR2 center = (mins + maxs) * 0.5;
  VECTOR2 diff = maxs - mins;
  REAL scale = (diff[0] > diff[1]) ? diff[0] : diff[1];
  for (unsigned int x = 0; x < nodes.size(); x++)
  {
    nodes[x] -= center;
    nodes[x] *= 1.0 / scale;

    // Triangle seems to flip things
    nodes[x][1] *= -1;
  }

  writePoly(argv[2]);
  return 1;
}
