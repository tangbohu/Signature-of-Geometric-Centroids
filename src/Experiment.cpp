///////////////////////////////////////////////////////////////
//
// Experiment.cpp
//
//   Experiment to Find Best Parameters for Local Voxelizer
//
// by Song Peng ( songpeng@ustc.edu.cn )
// 
// 02/May/2015
//
///////////////////////////////////////////////////////////////

#include "time.h"
//#include <GL/glut.h>
#include "Controls.h"
#include "Helper.h"
#include "HelpDefines.h"
#include "Scan.h"
#include "math3D.h"
#include "LVoxelizer.h"
#include "Register.h"
#include "Experiment.h"


//**************************************************************************************//
//                                   Initialization
//**************************************************************************************//

Experiment::Experiment()
{

}

Experiment::~Experiment()
{

}

void Experiment::InitScanPair(Scan *_sceneScan, Scan *_modelScan)
{
	sceneScan = _sceneScan;
	modelScan = _modelScan;
}




//**************************************************************************************//
//                               World and Scan Matrix
//**************************************************************************************//

//void Experiment::PerformExperiment()
//{
//	
//}

void Experiment::MatchScanPair()
{
	printf("computing descriptors for scene scan... \n");

	sceneScan->SampleScan();
	sceneScan->DescribeScan();
	
	printf("computing descriptors for model scan... \n");

	modelScan->SampleScan();
	modelScan->DescribeScan();

	printf("performing experiment to create RP-curve... \n");

	vector<float> ratioThresList;

	// For experiment: local features and combinations
	ratioThresList.push_back( 0.98 );
	ratioThresList.push_back( 0.94 );
	ratioThresList.push_back( 0.90 );
	ratioThresList.push_back( 0.85 );
	ratioThresList.push_back( 0.80 );
	ratioThresList.push_back( 0.75 );
	ratioThresList.push_back( 0.70 );
	ratioThresList.push_back( 0.60 );
	ratioThresList.push_back( 0.50 );
	ratioThresList.push_back( 0.40 );

	// For experiment: voxel grid resolution
	//ratioThresList.push_back( 0.98 );
	//ratioThresList.push_back( 0.90 );
	//ratioThresList.push_back( 0.85 );
	//ratioThresList.push_back( 0.80 );
	//ratioThresList.push_back( 0.75 );
	//ratioThresList.push_back( 0.70 );
	//ratioThresList.push_back( 0.60 );
	//ratioThresList.push_back( 0.50 );
	//ratioThresList.push_back( 0.40 );
	//ratioThresList.push_back( 0.30 );


	RPDotList.clear();

	for (int i=0; i<ratioThresList.size(); i++)
	{
		float distRatioThres = ratioThresList[i];
		printf("i=%d  ratioThres: %.2f   ", i, distRatioThres);

		RPDot RPPoint;
		FindHighMaches( distRatioThres, RPPoint );

		RPDotList.push_back( RPPoint );
	}

#ifdef USE_DESCRIPTOR_OURS
	char *fileName = "E:\\Project_Local\\ShapeMatching\\SourceCode\\test_ours.dat";
#endif 
#ifdef USE_DESCRIPTOR_3DSC
	char *fileName = "E:\\Project_Local\\ShapeMatching\\SourceCode\\test_3dsc.dat";
#endif 
#ifdef USE_DESCRIPTOR_SHOT
	char *fileName = "E:\\Project_Local\\ShapeMatching\\SourceCode\\test_shot.dat";
#endif 

#ifdef USE_DESCRIPTOR_SGC
	char *fileName = "E:\\Project_Local\\ShapeMatching\\SourceCode\\test_sgc.dat";
#endif 

	WriteRPCurveFile( fileName );
}	



void Experiment::FindHighMaches(float distRatioThres, RPDot &RPPoint)
{
	const float pointDistThres = 0.10; // TODO: need to be tuned
	//const float pointDistThres = 0.15; // TODO: need to be tuned

	int truePositiveNum = 0;
	int falsePositiveNum = 0;
	int totalPositiveNum = 0;
	int totalMatchNum = 0;

	//for (int i=0; i<sceneScan->sampDescriptors.size(); i++)
	for (int i=0; i<sceneScan->sampDescriptors.size(); i=i+10)
	{
		//totalPositiveNum++;

		int modelSampPtID;
		bool isMatch = FindModelHighMaches( i, distRatioThres, modelSampPtID );

		if ( isMatch )
		{
			totalMatchNum++;

			bool isRealMatch = ValidateMatch(i, modelSampPtID, pointDistThres);

			if ( isRealMatch )    truePositiveNum++;
			else                  falsePositiveNum++;

			if ( isRealMatch )    totalPositiveNum++;

			//if ( isRealMatch )
			//{
			//	vector<Point> sceneSampPts = sceneScan->sampler.GetSamplePoints();
			//	sceneScan->BuildLocalVoxelizer(sceneSampPts[i], LOCAL_SHAPE_RADIUS, true);

			//	vector<Point> modelSampPts = modelScan->sampler.GetSamplePoints();
			//	modelScan->BuildLocalVoxelizer(modelSampPts[modelSampPtID], LOCAL_SHAPE_RADIUS, true);
			//}
		}

		else
		{
			bool isRealMatch = ValidateMatch(i, modelSampPtID, pointDistThres);

			//if ( !isRealMatch )    totalPositiveNum++;
			if ( isRealMatch )    totalPositiveNum++;
		}
	}

	RPPoint.one_precision = falsePositiveNum / (float)totalMatchNum;
	RPPoint.recall = truePositiveNum / (float)totalPositiveNum;

	printf("1-precision %.2f  recall: %.2f   \n", RPPoint.one_precision, RPPoint.recall);
	printf("total: %d  match: %d true: %d  false: %d \n", totalPositiveNum, totalMatchNum, truePositiveNum, falsePositiveNum);

}

void Experiment::WriteRPCurveFile(const char *fileName)
{
	FILE *fp;
	if ( (fp=fopen(fileName, "w+")) == NULL )
	{
		printf("Error: file not exists! \n");	
	}
	else
	{
		// Scan number
		//fprintf( fp, "Scan Num: %d \n\n", RPDotList.size());

		// Align local matrix of each Scan 
		for (int i=0; i<RPDotList.size(); i++)
		{
			fprintf(fp, "%.3f %.3f \n", RPDotList[i].one_precision, RPDotList[i].recall);
		}

		fprintf(fp, "\n");
	}

	fclose(fp); 
}


bool Experiment::FindModelHighMaches(int sceneSampPtID, float distRatioThres, int &modelSampPtID)
{
	vector<float> sceneSampDesc = sceneScan->sampDescriptors[sceneSampPtID]->descVector;

	vector<double> distList;
	float minDist = MAX_FLOAT;
	for (int i=0; i<modelScan->sampDescriptors.size(); i++)
	{
		vector<float> modelSampDesc = modelScan->sampDescriptors[i]->descVector;

		float descDist = helpRegister.CompareDescriptors( sceneSampDesc, modelSampDesc );
		distList.push_back( descDist );
	}

	sortedIndices = BubbleSort(distList, true);
	//int modelSampPtID = sortedIndices[0];
	modelSampPtID = sortedIndices[0];

	//vector<Point> modelSampPts = modelScan->sampler.GetSamplePoints();
	//modelScan->BuildLocalVoxelizer(modelSampPts[modelSampPtID], LOCAL_SHAPE_RADIUS, true);

	float distRatio = distList[0] / distList[1];
	//printf("distRatio: %.2f \n", distRatio);

	//printf(" distList[0]: %.2f ", distList[0]);

	//if ( distList[0] < distRatioThres )    return true;
	//else                                   return false;

	if ( distRatio < distRatioThres )    return true;
	else								 return false;
}

bool Experiment::ValidateMatch(int sceneSampPtID, int modelSampPtID, float pointDistThres)
{
	Vector3f modelSampPt = modelScan->sampDescriptors[modelSampPtID]->keyPtLRF.keyPoint.pos;
	Vector3f sceneSampPt = sceneScan->sampDescriptors[sceneSampPtID]->keyPtLRF.keyPoint.pos;

	Vector3f newModelSampPt;
	MultiplyVector(modelSampPt, modelScan->alignLocalMat, newModelSampPt);

	Vector3f newSceneSampPt;
	MultiplyVector(sceneSampPt, sceneScan->alignLocalMat, newSceneSampPt);

	float pointDist = len(newSceneSampPt-newModelSampPt);

	//printf("pointDist: %.2f \n", pointDist);

	if ( pointDist < pointDistThres )	return true;
	else                                return false;
}

/*float Experiment::CompareDescriptors(int featureType, vector<float> modelDesc, vector<float> dataDesc)
{
	// Some descriptors could be empty
	if ( modelDesc.size() == 0 || dataDesc.size() == 0 )
	{
		return MAX_FLOAT;
	}

	/////////////////////////////////////////////////////////////
	// 1. Compute distance vector for all voxels

	vector<float> diffDescVec;

	if ( featureType == FEATURE_SURFACE_NORMAL )
	{
		for (int i=0; i<modelDesc.size(); i=i+3)
		{
			Vector3f modelNor;
			modelNor[0] = modelDesc[i];
			modelNor[1] = modelDesc[i+1];
			modelNor[2] = modelDesc[i+2];

			Vector3f dataNor;
			dataNor[0] = dataDesc[i];
			dataNor[1] = dataDesc[i+1];
			dataNor[2] = dataDesc[i+2];

			float value = modelNor DOT dataNor;
			diffDescVec.push_back( value );
		}
	}
	else
	{
		for (int i=0; i<modelDesc.size(); i++)
		{
			float value = modelDesc[i] - dataDesc[i];
			diffDescVec.push_back( value );
		}
	}


	/////////////////////////////////////////////////////////////
	// 2. Compute Euclidean distance for the two descriptors 

	float tempSum = 0;
	for (int i=0; i<diffDescVec.size(); i++)
	{
		tempSum += diffDescVec[i] * diffDescVec[i]; 
	}

	float descDist = sqrt( tempSum );
	//printf(" descDist: %.2f \n", descDist);

	return descDist;
}*/