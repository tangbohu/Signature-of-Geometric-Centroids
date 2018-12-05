///////////////////////////////////////////////////////////////
//
// ClipMeshInBox.cpp
//
//   Clip Mesh that is Inside a Box
//
// by Song Peng ( songpeng@ustc.edu.cn )
// 
// 07/May/2015
//
//
///////////////////////////////////////////////////////////////

#include "ClipMeshInBox.h"


//**************************************************************************************//
//                                Clip Mesh within a Box
//**************************************************************************************//

vector<Triangle*> ClipMeshInBox(vector<Triangle*> origTris, Vector3f boxCen, Vector3f boxSize)
{
	// some need to be deleted here
	vector<Triangle*> clipTris;
	clipTris = CopyTriangles( origTris );

	vector<Plane> clipPlanes = CreateClipPlanes(boxCen, boxSize);

	//int i=0;
	for (int i=0; i<clipPlanes.size(); i++)
	{
		clipTris = ClipTrisWithPlane(clipTris, clipPlanes[i].normal, clipPlanes[i].point);
	}

	return clipTris;
}

vector<Triangle*> CopyTriangles(vector<Triangle*> triList)
{
	vector<Triangle*> newTris;

	for (int i=0; i<triList.size(); i++)
	{
		Triangle *tri = new Triangle();
		tri->v[0] = triList[i]->v[0];
		tri->v[1] = triList[i]->v[1];
		tri->v[2] = triList[i]->v[2];
		tri->ComputeNormal();
		tri->ComputeArea();

		newTris.push_back( tri );
	}

	return newTris;
}

vector<Plane> CreateClipPlanes(Vector3f center, Vector3f size)
{
	vector<Plane> clipPlanes;
	Plane planePX, planePY, planePZ;
	planePX.normal = Vector3f( 1, 0, 0);   planePX.point  = center + 0.5f*size;
	planePY.normal = Vector3f( 0, 1, 0);   planePY.point  = center + 0.5f*size;
	planePZ.normal = Vector3f( 0, 0, 1);   planePZ.point  = center + 0.5f*size;

	Plane planeNX, planeNY, planeNZ;
	planeNX.normal = Vector3f(-1, 0, 0);   planeNX.point  = center - 0.5f*size;
	planeNY.normal = Vector3f( 0,-1, 0);   planeNY.point  = center - 0.5f*size;
	planeNZ.normal = Vector3f( 0, 0,-1);   planeNZ.point  = center - 0.5f*size;

	clipPlanes.push_back( planePX );
	clipPlanes.push_back( planePY );
	clipPlanes.push_back( planePZ );
	clipPlanes.push_back( planeNX );
	clipPlanes.push_back( planeNY );
	clipPlanes.push_back( planeNZ );

	return clipPlanes;
}




//**************************************************************************************//
//                                Clip Mesh using a Plane
//**************************************************************************************//

vector<Triangle*> ClipTrisWithPlane(vector<Triangle*> triList, Vector3f planeNor, Vector3f planePt)
{
	vector<Triangle*> outTriList;

	for (int i=0; i<triList.size(); i++)
	{
		vector<Triangle*> clippedTris = PlaneTriTest(planeNor, planePt, triList[i]);

		for (int j=0; j<clippedTris.size(); j++)
		{
			outTriList.push_back( clippedTris[j] );
		}
	}

	return outTriList;
}


vector<Triangle*> PlaneTriTest(Vector3f planeNor, Vector3f planePt, Triangle *triangle)
{
	// Compute signed distance between each vertex and the plane
	float pointDists[3];
	for (int i=0; i<3; i++)
	{
		float dotR0 = planeNor DOT (triangle->v[i]-planePt);
		pointDists[i] = dotR0;
	}

	// Save triangle vertices that are upper or lower than the plane (along plane normal direction)
	vector<TempPoint> upperPts;
	vector<TempPoint> lowerPts;
	for (int i=0; i<3; i++)
	{
		TempPoint tempPoint;
		tempPoint.point = triangle->v[i];
		tempPoint.dist  = fabs(pointDists[i]);

		if ( pointDists[i] > 0 )
			upperPts.push_back( tempPoint );	
		else
			lowerPts.push_back( tempPoint );
	}

	// Compute the two cross point (if exist) between the triangle and the plane
	Vector3f crossPtPair[2];
	int ptIndex = 0;
	for (int i=0; i<upperPts.size(); i++)
	{
		for (int j=0; j<lowerPts.size(); j++)
		{
			Vector3f crossPt = LerpPoint(upperPts[i].point, upperPts[i].dist, lowerPts[j].point, lowerPts[j].dist);
			crossPtPair[ptIndex] = crossPt;
			ptIndex++;
		}
	}


	vector<Triangle*> outTris;

	///////////////////////////////////////////////////////////////////////////
	// Case 1: Ignore the triangle if it is on top of the plane
	if ( upperPts.size() == 3 )
	{
		return outTris;
	}

	///////////////////////////////////////////////////////////////////////////
	// Case 2: Save the triangle if it is on top of the plane
	if ( lowerPts.size() == 3 )
	{
		outTris.push_back( triangle);
		return outTris;
	}

	///////////////////////////////////////////////////////////////////////////
	// Case 3: Save one clipped triangle if two vertices of the input triangle
	//         are on top of the plane (along the normal direction)
	if ( upperPts.size() == 2 )
	{
		Triangle *clipTri = new Triangle();
		clipTri->v[0] = crossPtPair[0];
		clipTri->v[1] = crossPtPair[1];
		clipTri->v[2] = lowerPts[0].point;
		clipTri->CorrectNormal( triangle->normal );
		clipTri->ComputeArea();

		outTris.push_back( clipTri );

		delete triangle;
	}

	///////////////////////////////////////////////////////////////////////////
	// Case 4: Save one clipped triangle if one vertex of the input triangle
	//         is on top of the plane (along the normal direction)
	else
	{
		Triangle *clipTri1 = new Triangle();
		clipTri1->v[0] = lowerPts[0].point;
		clipTri1->v[1] = lowerPts[1].point;
		clipTri1->v[2] = crossPtPair[0];
		clipTri1->CorrectNormal( triangle->normal );
		clipTri1->ComputeArea();

		Triangle *clipTri2 = new Triangle();
		clipTri2->v[0] = crossPtPair[0];
		clipTri2->v[1] = crossPtPair[1];
		clipTri2->v[2] = lowerPts[1].point;
		clipTri2->CorrectNormal( triangle->normal );
		clipTri2->ComputeArea();

		outTris.push_back( clipTri1 );
		outTris.push_back( clipTri2 );

		delete triangle;
	}

	return outTris;
}

Vector3f LerpPoint(Vector3f pointA, float distA, Vector3f pointB, float distB)
{
	float ratio = distA / (distA+distB);
	Vector3f newPt = pointA + ratio*(pointB-pointA);
	return newPt;
}
