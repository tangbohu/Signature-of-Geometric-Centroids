///////////////////////////////////////////////////////////////
//
// RayMeshTest.cpp
//
//   Ray-Mesh Intersection Test
//
// by Song Peng ( songpeng@ustc.edu.cn )
// 
// 02/Apr/2015
//
///////////////////////////////////////////////////////////////

#include "RayMeshTest.h"
#include "HelpDefines.h"
#include "HelpStructs.h"
#include "Helper.h"


//**************************************************************************************//
//                        Intersection Test for General Rays 
//**************************************************************************************//

bool LineMeshIntersect(Vector3f rayOrg, Vector3f rayDir, vector<Triangle*> _triList, HitPoint &hitPoint)
{
	if ( _triList.size() == 0 )
		return false;

	vec rayPt1 = rayOrg + 1.0f*rayDir;
	float minHitDist = MAX_INT;

	// Note: current exhaustive search will be slow for hig-resolution meshes
	for (int i=0; i<_triList.size(); i++) 
	{
		// Perform ray-triangle intersection test
		Triangle *tri = _triList[i];
		Vector3f tempHitPos;
		float currDist = LineHitTriangle(rayOrg, rayDir, *tri, tempHitPos, false);

		// Find the hit point with smallest hit distance
		if( currDist > 0 )
		{	
			if ( currDist < minHitDist )
			{
				minHitDist = currDist;
				hitPoint.dist = currDist;
				hitPoint.pos  = tempHitPos;
				hitPoint.nor  = tri->normal;
			}
		}
	}

	if ( minHitDist != MAX_INT )
		return true;
	else
		return false;
}




//**************************************************************************************//
//                       Intersection Test for Axis-aligned Rays
//**************************************************************************************//

vector<HitPoint> RayMeshIntersect(Vector3f rayOrg, Vector3f rayDir, vector<Triangle*> _triList)
{
	vector<HitPoint> hitPoints;

	//////////////////////////////////////////////////////
	// 1. Compute all the intersection points between the ray and the mesh

	for (int i=0; i<_triList.size(); i++)
	{
		Vector3f triMinPt; 
		Vector3f triMaxPt;
		GetTriangleBBox(*_triList[i], triMinPt, triMaxPt);

		if ( RayBoxIntersect(rayOrg, rayDir, triMinPt, triMaxPt) == false )
		{
			continue;
		}

		Vector3f hitPos;
		double hitDist;
		hitDist = LineHitTriangle(rayOrg, rayDir, *_triList[i], hitPos, false);

		if ( hitDist > 0.0 )
		{
			HitPoint hitPt;
			hitPt.dist = hitDist;
			hitPt.pos  = hitPos;
			hitPt.nor  = _triList[i]->normal;

			hitPoints.push_back( hitPt );
		}
	}
	//printf(" %d ", x);


	//////////////////////////////////////////////////////
	// 2. Remove duplicated intersection points

	RemoveDuplicatedPoints( hitPoints );


	//////////////////////////////////////////////////////
	// 3. Sort the intersection points according to their distance to rayOrg

	vector<double> hitDepths;
	for (int i=0; i<hitPoints.size(); i++)
	{
		hitDepths.push_back( hitPoints[i].dist );
	}

	vector<int> hitPtIndices = BubbleSort(hitDepths, true);
	vector<HitPoint> tempPoints = hitPoints;
	hitPoints.clear();
	for (int i=0; i<hitPtIndices.size(); i++)
	{
		hitPoints.push_back( tempPoints[hitPtIndices[i]]);
	}

	//if ( hitPoints.size() == 3 )
	//{
	//	printf("%.12f  %.12f  %.12f \n", hitDepths[0], hitDepths[1], hitDepths[2]);
	//
	//}

	return hitPoints;
}

void RemoveDuplicatedPoints( vector<HitPoint> &hitPtList)
{
	vector<double> depthList;
	vector<Vector3f> pointList;
	vector<Vector3f> normalList;

	for (int i=0; i<hitPtList.size(); i++)
	{
		depthList.push_back( hitPtList[i].dist );
		pointList.push_back( hitPtList[i].pos  );
		normalList.push_back( hitPtList[i].nor );
	}

	vector<double> depthListCopy = depthList;
	vector<Vector3f> pointListCopy = pointList;
	vector<Vector3f> normalListCopy = normalList;


	//printf("Before: ");
	//for (int i=0; i<depthList.size(); i++)
	//{
	//	printf("[i=%d %f] ", i, depthList[i]);
	//}
	//printf("\n");

	depthList.clear();
	pointList.clear();
	normalList.clear();
	for (int i=0; i<depthListCopy.size(); i++)
	{
		if ( IsDepthInList( depthListCopy[i], depthList) == false )
		{
			depthList.push_back( depthListCopy[i] );
			pointList.push_back( pointListCopy[i] );
			normalList.push_back( normalListCopy[i] );
		}
	}

	hitPtList.clear();

	for (int i=0; i<depthList.size(); i++)
	{
		HitPoint hitPt;
		hitPt.dist = depthList[i];
		hitPt.pos  = pointList[i];
		hitPt.nor  = normalList[i];

		hitPtList.push_back( hitPt );
	}

	//printf("After: ");
	//for (int i=0; i<depthList.size(); i++)
	//{
	//	printf("[i=%d %f] ", i, depthList[i]);
	//}
	//printf("\n\n");
}

bool IsDepthInList(double testDepth, vector<double> depthList)
{
	const double distThres = 0.00001;// Note: this value depends on the scale and shape of the 3D model 

	for (int i=0; i<depthList.size(); i++)
	{
		double depthDist = fabs( testDepth- depthList[i] );

		if ( depthDist < distThres  )
		{
			return true;
		}
	}

	return false;
}

void GetTriangleBBox(Triangle tri, Vector3f &boxMinPt, Vector3f &boxMaxPt)
{
	Vector3f v0 = tri.v[0];
	Vector3f v1 = tri.v[1];
	Vector3f v2 = tri.v[2];

	boxMinPt[0] = _MIN(v0[0], _MIN(v1[0], v2[0]));
	boxMinPt[1] = _MIN(v0[1], _MIN(v1[1], v2[1]));
	boxMinPt[2] = _MIN(v0[2], _MIN(v1[2], v2[2]));
	boxMaxPt[0] = _MAX(v0[0], _MAX(v1[0], v2[0]));
	boxMaxPt[1] = _MAX(v0[1], _MAX(v1[1], v2[1]));
	boxMaxPt[2] = _MAX(v0[2], _MAX(v1[2], v2[2]));
}

bool RayBoxIntersect(Vector3f rayOrg, Vector3f rayDir, Vector3f boxMinPt, Vector3f boxMaxPt)
{
	// X-axis aligned ray
	if ( rayDir == Vector3f(1,0,0) || rayDir == Vector3f(-1,0,0) )
	{
		if ( rayOrg[1] < boxMinPt[1] || rayOrg[1] > boxMaxPt[1] ||
			 rayOrg[2] < boxMinPt[2] || rayOrg[2] > boxMaxPt[2] )
		{
			return false;
		}
		else
		{
			return true;
		}
	}

	// Y-axis aligned ray
	if ( rayDir == Vector3f(0,1,0) || rayDir == Vector3f(0,-1,0) )
	{
		if ( rayOrg[0] < boxMinPt[0] || rayOrg[0] > boxMaxPt[0] ||
			 rayOrg[2] < boxMinPt[2] || rayOrg[2] > boxMaxPt[2] )
		{
			return false;
		}
		else
		{
			return true;
		}
	}

	// Z-axis aligned ray
	if ( rayDir == Vector3f(0,0,1) || rayDir == Vector3f(0,0,-1) )
	{
		if ( rayOrg[0] < boxMinPt[0] || rayOrg[0] > boxMaxPt[0] ||
			 rayOrg[1] < boxMinPt[1] || rayOrg[1] > boxMaxPt[1] )
		{
			return false;
		}
		else
		{
			return true;
		}
	}

	printf("Warning: The ray direction is not axis-aligned. \n");
	return false;
}