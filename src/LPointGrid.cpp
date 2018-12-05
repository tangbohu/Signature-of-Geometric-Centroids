///////////////////////////////////////////////////////////////
//
// LPointGrid.cpp
//
//   3D Point Grid Class
//
// by Song Peng ( songpeng@ustc.edu.cn )
// 
// 01/Apr/2015
//
///////////////////////////////////////////////////////////////

#include "Controls.h"
#include <GL/glut.h>
#include "RayMeshTest.h"
#include "LVoxel.h"
#include "LPointGrid.h"


//**************************************************************************************//
//                                   Initialization
//**************************************************************************************//

LPointGrid::LPointGrid()
{
	pointNum  = 0;
	pointList = NULL;
}

LPointGrid::~LPointGrid()
{
	if( pointList != NULL )	   
	{
		delete pointList;
		pointList = NULL;
	}
}

void LPointGrid::InitLPointGrid(Vector3f _bboxMinPt, Vector3f voxelSize, Vector3i volumeSize)
{
	bboxMinPt = _bboxMinPt;
	gridSize  = voxelSize / (float)(VOXEL_EDGE_SAMPLE_NUM-1);

	gridDimen[0] = (VOXEL_EDGE_SAMPLE_NUM-1) * volumeSize[0] + 1;
	gridDimen[1] = (VOXEL_EDGE_SAMPLE_NUM-1) * volumeSize[1] + 1;
	gridDimen[2] = (VOXEL_EDGE_SAMPLE_NUM-1) * volumeSize[2] + 1;

	if( pointList != NULL )	
	{
		delete pointList;
		pointList = NULL;
	}
	pointNum  = gridDimen[0] * gridDimen[1] * gridDimen[2];
	pointList = new SampPoint[pointNum];
}

void LPointGrid::SetXAxisHitPoints(vector<vector<HitPoint>> _XRayHitPts)
{
	XRayHitPts = _XRayHitPts;
}




//**************************************************************************************//
//                                 Compute Point Grid
//**************************************************************************************//

void LPointGrid::ComputePointGrid()
{
	// Note: the order of FOR loops is important
	for( int k=0; k< gridDimen[2]; k++ )
	for( int j=0; j< gridDimen[1]; j++ )
	for( int i=0; i< gridDimen[0]; i++ )
	{
		Vector3i ptCoord = Vector3i(i, j, k);
		int index = GetPointIndex( ptCoord );

		pointList[index].state = POINT_UNKNOWN;
		pointList[index].coord = ptCoord;
		pointList[index].pos = GetPointPosition( ptCoord );
	}
}

void LPointGrid::ComputePointState()
{
	for( int j=0; j<gridDimen[1]; j++ )
	for( int k=0; k<gridDimen[2]; k++ )
	{
		vector<HitPoint> hitPoints = GetRayHitPoints(Vector3i(0,j,k), Vector3i(1,0,0));

		// Case 1: if there is no hit point, all sampling point are outside of the mesh
		if ( hitPoints.size() == 0  )
		{
			for (int i=0; i<gridDimen[0]; i++)
			{
				int index = GetPointIndex(Vector3i(i,j,k));
				pointList[index].state = POINT_OUT_MESH;
			}
		}


		// Case 2: if there is more than 1 hit point, only consider the sampling point from the origin to
		//         the 1st hit point (note: here exists simplifications for cases with 2 or more hit points)
		else
		{
			for (int i=0; i<gridDimen[0]; i++)
			{
				int index = GetPointIndex(Vector3i(i,j,k));

				if ( pointList[index].pos[0] < hitPoints[0].pos[0] )
				{
					pointList[index].state = POINT_IN_MESH;
				}
				else
				{
					pointList[index].state = POINT_OUT_MESH;
				}
			}

			//if ( hitPoints.size() > 1 )
			//{
			//	printf("Waring: there are %d hit points. \n", hitPoints.size());
			//}
		}

		//else if ( hitPoints.size() == 2 )
		//{
		//	for (int i=0; i<gridDimen[0]; i++)
		//	{
		//		int index = GetPointIndex(Vector3i(i,j,k));

		//		//float angle = acos(hitPoints[0].nor DOT Vector3f(1,0,0)) * 180/M_PI; // 0.86 => 36 degree

		//		if ( pointList[index].pos[0] < hitPoints[0].pos[0] ||
		//			 pointList[index].pos[0] > hitPoints[1].pos[0] )
		//		{
		//			pointList[index].state = POINT_OUT_MESH;
		//		}

		//		else
		//		{
		//			pointList[index].state = POINT_IN_MESH;
		//		}
		//	}

		//	printf("Warning: there are 2 intersecting points with the mesh. \n");
		//}

		//else
		//{
		//	printf("Warning: there are more than 2 intersecting points with the mesh. \n");
		//}

		//int index = j*(volumeSize.z+1)+k;
		//printf("[%d %d] ", index, hitPtLists.size());
	}
}

vector<SampPoint> LPointGrid::GetPointsInsideVoxel(Vector3i voxelPos)
{
	vector<SampPoint> insidePoints; 

	for (int k=0; k<VOXEL_EDGE_SAMPLE_NUM; k++)
	for (int j=0; j<VOXEL_EDGE_SAMPLE_NUM; j++)
	for (int i=0; i<VOXEL_EDGE_SAMPLE_NUM; i++)
	{
		Vector3i pointCoord = voxelPos*(VOXEL_EDGE_SAMPLE_NUM-1) + Vector3i(i,j,k);
		int pointIndex = GetPointIndex( pointCoord );

		insidePoints.push_back( pointList[pointIndex] );
	}

	return insidePoints;
}




//**************************************************************************************//
//                                 Utility Functions
//**************************************************************************************//

int LPointGrid::GetPointIndex(Vector3i pointCoord)
{
	int index = pointCoord[2]*gridDimen[1]*gridDimen[0] + pointCoord[1]*gridDimen[0] + pointCoord[0];
	return index;	
}

Vector3f LPointGrid::GetPointPosition(Vector3i pointCoord)
{
	Vector3f pickPoint;

	pickPoint[0] = bboxMinPt[0] + pointCoord[0]*gridSize[0]; 
	pickPoint[1] = bboxMinPt[1] + pointCoord[1]*gridSize[1]; 
	pickPoint[2] = bboxMinPt[2] + pointCoord[2]*gridSize[2];

	//printf("coord [%d %d %d]  pos [%.3f %.3f %.3f] \n", pointCoord.x, pointCoord.y, pointCoord.z, pickPoint.x, pickPoint.y, pickPoint.z);
	return pickPoint;
}

vector<HitPoint> LPointGrid::GetRayHitPoints(Vector3i pointCoord, Vector3i rayDir)
{
	vector<HitPoint> hitPoints;

	if ( rayDir != Vector3i(1,0,0) )
	{
		printf("Warning: The ray direction supposes to be along x-axis. \n");
		return hitPoints;
	}

	int index = pointCoord[1] * gridDimen[2] + pointCoord[2];
	hitPoints = XRayHitPts[index];

	return hitPoints;
}




//**************************************************************************************//
//                                 Draw Point Grid
//**************************************************************************************//

void LPointGrid::DrawPointGrid(float pointSize, Vector3f pointColor)
{
	if ( pointNum == 0 )
		return;

	glDisable( GL_LIGHTING );
	glPointSize( pointSize );
	glColor3f( pointColor[0], pointColor[1], pointColor[2] );

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();

	glBegin(GL_POINTS);
	for (int i=0; i<pointNum; i++ )
	{
		if ( pointList[i].state == POINT_UNKNOWN )
			glColor3f( 0.1, 0.1, 0.8 );
		else if( pointList[i].state == POINT_OUT_MESH )
			glColor3f( 0.8, 0.8, 0.1 );
		else
			glColor3f( pointColor[0], pointColor[1], pointColor[2] );

		if ( pointList[i].state == POINT_IN_MESH ||
			 pointList[i].state == POINT_UNKNOWN )
		{
			glVertex3f(pointList[i].pos[0], pointList[i].pos[1], pointList[i].pos[2]);
		}
	}
	glEnd();

	glPopMatrix();
}

void LPointGrid::DrawCrossPoints()
{
	if ( XRayHitPts.size() == 0 )
		return;

	glDisable(GL_LIGHTING);
	glPointSize(4.0);
	glLineWidth(2.0);

	for( int j=0; j<gridDimen[1]; j++ )
	for( int k=0; k<gridDimen[2]; k++ )
	{
		Vector3f point = GetPointPosition(Vector3i(0,              j, k));
		Vector3f endPt = GetPointPosition(Vector3i(gridDimen[0]-1, j, k));

		int index = j*gridDimen[2] + k;
		vector<HitPoint> hitPoints = XRayHitPts[index];

		if ( hitPoints.size() == 0 )
			continue;

		// Draw the X-axis aligned ray
		//glColor3f(0.1,0.1,0.9);
		//glBegin(GL_LINES);
		//glVertex3f(point.x, point.y, point.z);
		//glVertex3f(endPt.x, endPt.y, endPt.z);
		//glEnd();

		// Draw the hit points along the ray
		glPointSize(8.0);
		glBegin(GL_POINTS);
		for (int i=0; i<hitPoints.size(); i++)
		{
			if ( i %2 == 0 )  	    glColor3f(0.9,0.1,0.1);
			else					glColor3f(0.2,0.2,0.9);

			//printf("i=%d [%.3f %.3f %.3f] \n", i, hitPoints[i].x, hitPoints[i].y, hitPoints[i].z);
			glVertex3f(hitPoints[i].pos[0], hitPoints[i].pos[1], hitPoints[i].pos[2]);
		}
		glEnd();
	}

	glPointSize(1.0);
	glLineWidth(1.0);
	glEnable(GL_LIGHTING);
}


