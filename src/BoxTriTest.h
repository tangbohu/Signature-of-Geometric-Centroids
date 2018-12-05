///////////////////////////////////////////////////////////////
//
// BoxTriTest.h
//
//   Intersection Test between Cube and Triangle
//
// by Song Peng ( song0083@ntu.edu.sg )
// 
// 11/Oct/2013
//
//
///////////////////////////////////////////////////////////////

#ifndef _CUBE_TRI_TEST_H
#define _CUBE_TRI_TEST_H

#define  TRIANGLE_OUT_BOX       0
#define  TRIANGLE_CROSS_BOX     1
#define  TRIANGLE_IN_BOX        2

#include "Vec.h"
#include <Vector>

using namespace std;
using namespace trimesh;


///////////////////////////////////////////////////////////////
// Function Declaration
///////////////////////////////////////////////////////////////

int BoxTriangleIntersect(Vector3f boxCenter, Vector3f boxSize, Vector3f triVertices[3]);
bool IsPointInsideBox(Vector3f boxMinPt, Vector3f boxMaxPt, Vector3f point);
bool BoxTriangleSATest(Vector3f boxMinPt, Vector3f boxMaxPt, Vector3f triVertices[3]);
float getmin(const vector<Vector3f> &points, Vector3f axis);
float getmax(const vector<Vector3f> &points, Vector3f axis);
bool isect(const vector<Vector3f> &points1, const vector<Vector3f> &points2, Vector3f axis);


#endif