///////////////////////////////////////////////////////////////
//
// BoxTriTest.h
//
//   Intersection Test between Cube and Triangle
//
// by Song Peng ( song0083@ntu.edu.sg )
// 
// 13/Oct/2012
//
///////////////////////////////////////////////////////////////

#include "BoxTriTest.h"


//**************************************************************************************//
//                             Vector, Matrix and Array
//**************************************************************************************//

int BoxTriangleIntersect(Vector3f boxCenter, Vector3f boxSize, Vector3f triVertices[3])
{
	//printf("center: [%.3f %.3f %.3f] \n", boxCenter[0], boxCenter[1], boxCenter[2]);
	//printf("size:   [%.3f %.3f %.3f] \n", boxSize[0], boxSize[1], boxSize[2]);
	//printf("v0:     [%.3f %.3f %.3f] \n", triVertices[0][0], triVertices[0][1], triVertices[0][2]);
	//printf("v1:     [%.3f %.3f %.3f] \n", triVertices[1][0], triVertices[1][1], triVertices[1][2]);	
	//printf("v2:     [%.3f %.3f %.3f] \n", triVertices[2][0], triVertices[2][1], triVertices[2][2]);

	Vector3f boxMinPt = boxCenter - 0.5f * boxSize;
	Vector3f boxMaxPt = boxCenter + 0.5f * boxSize;

	if ( IsPointInsideBox(boxMinPt, boxMaxPt, triVertices[0]) && 
		 IsPointInsideBox(boxMinPt, boxMaxPt, triVertices[1]) &&
		 IsPointInsideBox(boxMinPt, boxMaxPt, triVertices[2]) )
	{
		return TRIANGLE_IN_BOX;
	}

	else
	if ( IsPointInsideBox(boxMinPt, boxMaxPt, triVertices[0]) || 
		 IsPointInsideBox(boxMinPt, boxMaxPt, triVertices[1]) ||
		 IsPointInsideBox(boxMinPt, boxMaxPt, triVertices[2]) )
	{
		return TRIANGLE_CROSS_BOX;
	}

	else
	{
		if( BoxTriangleSATest(boxMinPt, boxMaxPt, triVertices) )
			return TRIANGLE_CROSS_BOX;
		else
			return TRIANGLE_OUT_BOX;
	}
}

bool IsPointInsideBox(Vector3f boxMinPt, Vector3f boxMaxPt, Vector3f point)
{
	if ( point[0] > boxMinPt[0] && point[0] < boxMaxPt[0] &&
		 point[1] > boxMinPt[1] && point[1] < boxMaxPt[1] &&
		 point[2] > boxMinPt[2] && point[2] < boxMaxPt[2] )
	{
		return true;
	}
	else
	{
		return false;
	}
}

//Separating Axis Test for Triangle and Box
bool BoxTriangleSATest(Vector3f boxMinPt, Vector3f boxMaxPt, Vector3f triVertices[3])
{
	vector<Vector3f> boxpoints;
	boxpoints.push_back( Vector3f(boxMaxPt[0], boxMaxPt[1], boxMaxPt[2]) );
	boxpoints.push_back( Vector3f(boxMaxPt[0], boxMaxPt[1], boxMinPt[2]) );
	boxpoints.push_back( Vector3f(boxMaxPt[0], boxMinPt[1], boxMaxPt[2]) );
	boxpoints.push_back( Vector3f(boxMaxPt[0], boxMinPt[1], boxMinPt[2]) );
	boxpoints.push_back( Vector3f(boxMinPt[0], boxMaxPt[1], boxMaxPt[2]) );
	boxpoints.push_back( Vector3f(boxMinPt[0], boxMaxPt[1], boxMinPt[2]) );
	boxpoints.push_back( Vector3f(boxMinPt[0], boxMinPt[1], boxMaxPt[2]) );
	boxpoints.push_back( Vector3f(boxMinPt[0], boxMinPt[1], boxMinPt[2]) );

	vector<Vector3f> tripoints;
	tripoints.push_back( triVertices[0] );
	tripoints.push_back( triVertices[1] );
	tripoints.push_back( triVertices[2] );

	// test the x, y, and z axes
	if (!isect(boxpoints, tripoints, Vector3f(1, 0, 0))) return false;
	if (!isect(boxpoints, tripoints, Vector3f(0, 1, 0))) return false;
	if (!isect(boxpoints, tripoints, Vector3f(0, 0, 1))) return false;

	// test the triangle normal
	Vector3f triedge1 = tripoints[1] - tripoints[0];
	Vector3f triedge2 = tripoints[2] - tripoints[1];
	Vector3f trinormal = triedge1 CROSS triedge2;
	if (!isect(boxpoints, tripoints, trinormal)) return false;

	// test the 9 edge cross products
	Vector3f triedge3 = tripoints[0] - tripoints[2];

	Vector3f boxedge1 = Vector3f(1, 0, 0);
	Vector3f boxedge2 = Vector3f(0, 1, 0);
	Vector3f boxedge3 = Vector3f(0, 0, 1);

	if (!isect(boxpoints, tripoints, boxedge1 CROSS triedge1)) return false;
	if (!isect(boxpoints, tripoints, boxedge1 CROSS triedge2)) return false;
	if (!isect(boxpoints, tripoints, boxedge1 CROSS triedge3)) return false;

	if (!isect(boxpoints, tripoints, boxedge2 CROSS triedge1)) return false;
	if (!isect(boxpoints, tripoints, boxedge2 CROSS triedge2)) return false;
	if (!isect(boxpoints, tripoints, boxedge2 CROSS triedge3)) return false;

	if (!isect(boxpoints, tripoints, boxedge3 CROSS triedge1)) return false;
	if (!isect(boxpoints, tripoints, boxedge3 CROSS triedge2)) return false;
	if (!isect(boxpoints, tripoints, boxedge3 CROSS triedge3)) return false;
}


float getmin(const vector<Vector3f> &points, Vector3f axis)
{
	//float min = std::numeric_limits<float>::max(); 
	const float FLOAT_MAX = 100000000;
	float min = FLOAT_MAX;

	for (int ctr = 0; ctr < points.size(); ctr++)
	{
		float dotprod = points[ctr] DOT axis;
		if (dotprod < min) min = dotprod;
	}
	return min;
}

float getmax(const vector<Vector3f> &points, Vector3f axis)
{
	//float max = -std::numeric_limits<float>::max(); 
	const float FLOAT_MIN = -100000000;
	float max = FLOAT_MIN;

	for (int ctr = 0; ctr < points.size(); ctr++)
	{
		float dotprod = points[ctr] DOT axis;
		if (dotprod > max) max = dotprod;
	}
	return max;
}

bool isect(const vector<Vector3f> &points1, const vector<Vector3f> &points2, Vector3f axis)
{
	if (getmin(points1, axis) > getmax(points2, axis)) return false;
	if (getmax(points1, axis) < getmin(points2, axis)) return false;
	return true;     
}

