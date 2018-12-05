///////////////////////////////////////////////////////////////
//
// ScanSampler.cpp
//
//   Sampling Scan Mesh Vertices Class
//
// by Song Peng ( songpeng@ustc.edu.cn )
// 
// 11/May/2015
//
///////////////////////////////////////////////////////////////

#include "time.h"
#include "Helper.h"
#include "Controls.h"
#include "Scan.h"
#include "ScanSampler.h"


//**************************************************************************************//
//                                   Initialization
//**************************************************************************************//

ScanSampler::ScanSampler()
{

}

ScanSampler::~ScanSampler()
{

}

vector<Point> ScanSampler::GetSamplePoints()
{
	return meshSamplePts;
}

vector<Point> ScanSampler::GetFeaturePoints()
{
	return meshFeaturePts;
}

vector<int> ScanSampler::GetFeaturePointIndices()
{
	return featPtIndices;
}




//**************************************************************************************//
//                        Mesh Vertex Sampling (uniform sampling)
//**************************************************************************************//

void ScanSampler::SampleMesh_Feature(vector<Point> verList)
{
	vector<int> randPtIndices;
	meshFeaturePts = GetRandomPointList(meshSamplePts, FEATURE_POINT_NUM, FEATURE_POINT_DIST, randPtIndices);
	featPtIndices  = randPtIndices;
}

void ScanSampler::SampleMesh_Uniform(vector<Point> verList)
{
	vector<int> randPtIndices;

	meshSamplePts = GetRandomPointList(verList, SAMPLE_POINT_NUM, SAMPLE_POINT_DIST, randPtIndices);
}




//**************************************************************************************//
//                             Sample Point From a List
//**************************************************************************************//

// TODO: use a more general way for randomly sampling the mesh points
vector<Point> ScanSampler::GetRandomPointList(vector<Point> pointList, int sampleNum, const float distThres, vector<int> &randPtIndices)
{
	vector<Point> samplePtList;

	if ( sampleNum > pointList.size() || sampleNum <= 0 )
	{
		printf("Warning: The sampleNum %d should not be larger than pointNum %d\n\n", sampleNum, pointList.size());
		return samplePtList;
	}

	const int maxIterNum = 1000;   // Maximum iteration number to find a new sample point

	int iterNum = 0;
	while( samplePtList.size() < sampleNum )
	{
		int tempIndex = GetRandomObjIndex( pointList.size() );
		Point currPt = pointList[tempIndex];

		if ( IsValidSample(samplePtList, currPt, distThres) )
		{
			samplePtList.push_back( currPt );
			randPtIndices.push_back( tempIndex );  // Save the index of the sampled point

			iterNum = 0;
		}
		else
		{
			iterNum++;
		}

		if ( iterNum > maxIterNum )
		{
			printf("Warning: cannot find a new sample point and exit the loop with %d samplePt. \n", samplePtList.size());
			break;
		}
	}

	return samplePtList;
}

int ScanSampler::GetRandomObjIndex(int objNum)
{
	float randValue = objNum * ( rand()/(RAND_MAX+1.0) );

	if ( randValue == objNum )
		randValue = objNum - 1;
		

	return randValue;
}

bool ScanSampler::IsValidSample(vector<Point> samplePts, Point currPt, const float distThres)
{
	for (int i=0; i<samplePts.size(); i++)
	{
		float ptDist = len(samplePts[i].pos-currPt.pos);

		if ( ptDist < distThres)
		{
			return false;
		}
	}

	return true;
}




//**************************************************************************************//
//                                 Draw Mesh Samples 
//**************************************************************************************//

void ScanSampler::DrawFeaturePoints()
{
	DrawPointList( meshFeaturePts, 8.0, Vector3f(0.9,0.1,0.1) );
}

void ScanSampler::DrawSamplePoints()
{
	DrawPointList( meshSamplePts, 6.0, Vector3f(0.2,0.2,0.9) );
}

void ScanSampler::DrawPointList(vector<Point> pintList, int pointSize, Vector3f pointColor)
{
	if ( pintList.size() == 0 )
		return;

	glDisable( GL_LIGHTING );
	glColor3f(pointColor[0], pointColor[1], pointColor[2]);
	glLineWidth( 2.0 );
	glPointSize( pointSize );

	// Set the offset value
	const float epsilon = 1e-3;
	float projMat[16];
	glGetFloatv(GL_PROJECTION_MATRIX, projMat);

	// Offset the drawing by modifying the projection matrix
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glTranslatef(0.0f, 0.0f, -epsilon);
	glMultMatrixf(projMat);

	// Draw mesh sample points
	for (int i=0; i<pintList.size(); i++)
	{
		glBegin(GL_POINTS);
		glVertex3f(pintList[i].pos[0], pintList[i].pos[1], pintList[i].pos[2]);
		glEnd();
	}

	glPopMatrix();

	glMatrixMode(GL_MODELVIEW);
	glLineWidth( 1.0 );
	glPointSize( 1.0 );
	glEnable(GL_LIGHTING);
}