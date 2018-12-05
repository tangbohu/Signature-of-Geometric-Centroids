///////////////////////////////////////////////////////////////
//
// Register.cpp
//
//   Class for Registering Scans 
//
// by Song Peng ( songpeng@ustc.edu.cn )
// 
// 11/May/2015
//
// TODOs:
// 1. Speed up the scan desciption and matching
// 2. Handle scans with a small overlapping
// 3. Retain the debugging capabilities
//
///////////////////////////////////////////////////////////////

#include "time.h"
#include <GL/glut.h>
#include "KDtree.h"
#include "ICP.h"
#include "XForm.h"
#include "math3D.h"
#include "Helper.h"
#include "HelpDefines.h"
#include "Controls.h"
#include "Scan.h"
#include "ScanAlign.h"
#include "LVoxelizer.h"
#include "Register.h"

// For Debug
#define HIGH_MATCH_NUM                 20


//**************************************************************************************//
//                                   Initialization
//**************************************************************************************//

Register::Register()
{

}

Register::~Register()
{

}

void Register::InitRegister(vector<Scan*> _scanList)
{
	scanList = _scanList;
}




//**************************************************************************************//
//                          Experiment: Scan Pair Overlap
//**************************************************************************************//

void Register::Experiment_Overlap(int startScanID, int endScanID)
{
	if ( scanList.size() == 0 )
		return;

	startScanID = 0;
	endScanID   = scanList.size()-1;

	/////////////////////////////////////////////////////////
	// 1. Describe each scan using local voxelizer

	DescribeScans( startScanID, endScanID );
	SetFixedScan( scanList[startScanID] );


	/////////////////////////////////////////////////////////
	// 2. Loop through all the scan pairs

	const float realOverlapThres  = 0.02;
	const float rotErrorThres     = 10;
	const float tranErrorThres    = 0.06;

	vector<ScanOverlap> overlapList;

	for (int i=startScanID; i<endScanID; i++)
	{
		for (int j=i+1; j<=endScanID; j++)
		{
			int tagtScanID = i;
			int dataScanID = j;

			// 2.1 Keep alignLocalMat of the two scans  
			double tagtRealMat[16];
			double dataRealMat[16];
			memcpy(tagtRealMat, scanList[tagtScanID]->alignLocalMat, sizeof(double)*16);
			memcpy(dataRealMat, scanList[dataScanID]->alignLocalMat, sizeof(double)*16);

			// 2.2 Computed the actual overlap for the scan pair
			float realOverlap = EvaluateAlign( scanList[tagtScanID],  scanList[dataScanID] );

			// 2.3 Align the scan pair if they have a certain amount of overlap
			AlignError alignError;
			float compOverlap = 0;
			float rotError    = 0;
			float tranError   = 0;
			int isScanAligned = 0;
			if ( realOverlap > realOverlapThres )
			{
				alignError = AlignScanWithTarget(scanList[dataScanID], scanList[tagtScanID], false);
				compOverlap = EvaluateAlign( scanList[tagtScanID],  scanList[dataScanID] );

				// 2.4 Compare alignment matrix with the ground truth
				ScanAlign scanAlign;
				scanAlign.CompareAlignMatrix(scanList[dataScanID]->alignLocalMat, dataRealMat, rotError, tranError);

				// 2.5 Determine if the alignment is counted as correct
				if ( rotError < rotErrorThres && tranError < tranErrorThres )
				{
					isScanAligned = true;
				}
			}

			// 2.6 Reset alignLocalMat of the two scans  
			memcpy(scanList[tagtScanID]->alignLocalMat, tagtRealMat, sizeof(double)*16);
			memcpy(scanList[dataScanID]->alignLocalMat, dataRealMat, sizeof(double)*16);

			printf("[%d %d]  real: %.3f  comp: %.3f  rot: %6.3f  tran: %.3f  aligned: %d \n", 
				i, j, realOverlap, compOverlap, rotError, tranError, isScanAligned);

			// 2.7 Save the overlap info
			ScanOverlap overlap;
			overlap.tagtScanID  = tagtScanID;
			overlap.dataScanID  = dataScanID;
			overlap.isAligned   = isScanAligned;
			overlap.realOverlap = realOverlap;
			overlap.compOverlap = compOverlap;
			overlap.rotError    = rotError;
			overlap.tranError   = tranError;

			overlapList.push_back( overlap );
		}
	}


	/////////////////////////////////////////////////////////
	// 3. Write the pairwise alignment result

	char *rawFileName  = "E:\\Project_Local\\ShapeMatching\\SourceCode\\ScanOverlaps_Chicken_S8.dat";
	char *simpFileName = "E:\\Project_Local\\ShapeMatching\\SourceCode\\ScanOverlaps_Chicken_S8_simple.dat";

	WriteScanOverlapFile( rawFileName, overlapList );
	WriteScanOverlapFile_Simple( simpFileName, overlapList );
}

void Register::WriteScanOverlapFile(const char *fileName, vector<ScanOverlap> overlapList)
{
	FILE *fp;
	if ( (fp=fopen(fileName, "w+")) == NULL )
	{
		printf("Error: file not exists! \n");	
	}
	else
	{
		// Error list size
		fprintf( fp, "Overlap List Size: %d \n\n", overlapList.size());

		// Align local matrix of each Scan 
		for (int i=0; i<overlapList.size(); i++)
		{
			int tagtScanID    = overlapList[i].tagtScanID;
			int dataScanID    = overlapList[i].dataScanID;
			float realOverlap = overlapList[i].realOverlap;
			float compOverlap = overlapList[i].compOverlap;
			float rotError    = overlapList[i].rotError;
			float tranError   = overlapList[i].tranError;
			int isAligned     = overlapList[i].isAligned;

			fprintf(fp, "[%2d %2d]  real: %.4f  comp: %.4f  rot: %8.4f  tran:%.4f  aligned: %d  \n",
				    tagtScanID, dataScanID, realOverlap, compOverlap, rotError, tranError, isAligned);
		}

		fprintf(fp, "\n");
	}

	fclose(fp); 
}

void Register::WriteScanOverlapFile_Simple(const char *fileName, vector<ScanOverlap> overlapList)
{
	FILE *fp;
	if ( (fp=fopen(fileName, "w+")) == NULL )
	{
		printf("Error: file not exists! \n");	
	}
	else
	{
		// Error list size
		fprintf( fp, "Overlap List Size: %d \n\n", overlapList.size());

		// Align local matrix of each Scan 
		for (int i=0; i<overlapList.size(); i++)
		{
			float realOverlap = overlapList[i].realOverlap;
			int isAligned     = overlapList[i].isAligned;

			fprintf(fp, "%.4f %d \n", realOverlap, isAligned);
		}

		fprintf(fp, "\n");
	}

	fclose(fp); 
}





//**************************************************************************************//
//                          Experiment: Alignment Accuracy
//**************************************************************************************//

void Register::Experiment_Accuracy(int startScanID, int endScanID)
{
	if ( scanList.size() == 0 )
		return;

	// 1. Describe each scan using local voxelizer
	DescribeScans( startScanID, endScanID );
	SetFixedScan( scanList[startScanID] );

	// 2. Pairwise match neighboring scans
	vector<AlignError> errorList;
	for (int i=startScanID; i<_MIN((int)scanList.size()-1, endScanID); i++)
	{
		printf("Aligning scan %2d with scan %d...\n", i, i+1);
		AlignError alignError = AlignScanWithTarget(scanList[i+1], scanList[i], false);

		printf("i=%d align rotError: %.3f  tranError: %.3f \n", i, alignError.rotError, alignError.tranError);

		errorList.push_back( alignError );
	}

	// 3. Report scan alignment errors
	printf("\n==========Align Error List============\n");
	for (int i=0; i<errorList.size(); i++)
	{
		Match match = errorList[i].match;

		printf("i=%2d  [%2d %2d]  rot: %.3f  tran: %.3f  dist: %.2f  init: %.2f  align: %.2f \n", 
			i, match.tagtScanID, match.dataScanID, errorList[i].rotError, errorList[i].tranError, match.descDist, match.initScore, match.alignScore);
	}

	// 4. Write scan alignment error file
	char *fileName = "E:\\Project_Local\\ShapeMatching\\SourceCode\\AlignErrors_.dat";
	WriteAlignErrorFile( fileName, errorList );
}


void Register::WriteAlignErrorFile(const char *fileName, vector<AlignError> errorList)
{
	FILE *fp;
	if ( (fp=fopen(fileName, "w+")) == NULL )
	{
		printf("Error: file not exists! \n");	
	}
	else
	{
		// Error list size
		fprintf( fp, "Error List Size: %d \n\n", errorList.size());

		// Align local matrix of each Scan 
		for (int i=0; i<errorList.size(); i++)
		{
			Match match     = errorList[i].match;
			float rotError  = errorList[i].rotError;
			float tranError = errorList[i].tranError;

			fprintf(fp, "[%2d %2d]  rot: %.4f  tran:%.4f  dist: %.3f  initScore: %.3f  alignScore: %.3f \n",
				match.tagtScanID, match.dataScanID, rotError, tranError, 
				match.descDist,   match.initScore,  match.alignScore);
		}

		fprintf(fp, "\n");
	}

	fclose(fp); 
}




//**************************************************************************************//
//                          Align Data Scan to a Target Scan
//**************************************************************************************//

void Register::AlignScanPair()
{
	if ( scanList.size() != 2 )
		return;

	int tagtScanID = 0;
	int dataScanID = 1;

	DescribeScans( tagtScanID, dataScanID );
	SetFixedScan( scanList[tagtScanID] );

	AlignError alignError = AlignScanWithTarget(scanList[dataScanID], scanList[tagtScanID], false);
	printf("align rotError: %.3f  tranError: %.3f \n", alignError.rotError, alignError.tranError);
}

void Register::UpdateLVoxelizers(Match match)
{
	int tagtScanID   = match.tagtScanID;
	int tagtSampPtID = match.tagtSampPtID;
	int dataScanID   = match.dataScanID;
	int dataFeatPtID = match.dataFeatPtID;

	Descriptor *tagtDesc = scanList[tagtScanID]->sampDescriptors[tagtSampPtID];
	Descriptor *dataDesc = scanList[dataScanID]->featDescriptors[dataFeatPtID];

	scanList[tagtScanID]->UpdateLVoxelizer( tagtDesc );
	scanList[dataScanID]->UpdateLVoxelizer( dataDesc );
}


AlignError Register::AlignScanWithTarget(Scan *dataScan, Scan *tagtScan, bool useICP)
{
	clock_t beginTime = clock();

	AlignError alignError;


	//////////////////////////////////////////////////////////////
	// 1. Set target scan as the only aligned scan (reuse some functions)

	alignedScans.clear();
	alignedScans.push_back( tagtScan );

	alignedSampPts.clear();
	vector<Point> scanSampPts = tagtScan->sampler.GetSamplePoints();
	for ( int i=0; i<scanSampPts.size(); i++ )
	{
		Vector3f outPtPos;
		MultiplyVector(scanSampPts[i].pos, tagtScan->alignLocalMat, outPtPos);

		alignedSampPts.push_back( outPtPos );
	}


	//////////////////////////////////////////////////////////////
	// 2. Compute local voxelizers for feature points on the data and match them with that of the model

	float maxAlignScore = MIN_FLOAT;
	Match outMatch;

	for (int i=0; i<dataScan->featDescriptors.size(); i++)
	{
		// Find matched descriptor candidates
		vector<Match> matchList = GetMatchedDescriptors(dataScan->scanID, i);
		if ( matchList.size() == 0 )
			continue;

		// Select best matched descriptor
		Match bestMatch = SelectBestMatch(matchList, dataScan, i, true);

		// Update best alignment if current match is better 
		if ( bestMatch.initScore > maxAlignScore )
		{
			maxAlignScore = bestMatch.initScore;
			outMatch      = bestMatch;

			if ( bestMatch.initScore > 0.1f )
			{
				printf("i=%d align: %.3f  desc: %.3f [initScore: %.3f]\n", i, bestMatch.alignScore, bestMatch.descDist, bestMatch.initScore);
			}
		}
	}

	if ( maxAlignScore == MIN_FLOAT )
	{
		return alignError;
	}


	//////////////////////////////////////////////////////////////
	// 3. Align data to the model using the best match

	printf("start aligning the scan...\n");

	alignError = PerformAlignUsingMatch( outMatch, useICP ); // Note: here only refine alignment with the target scan using ICP

#ifndef USE_DESCRIPTOR_OURS
	UpdateLVoxelizers( outMatch );
#endif

	clock_t endTime = clock();
	float elapsed = ((float) (endTime - beginTime)) / CLOCKS_PER_SEC; // Time unit: second
	printf("=> Matching scan %2d takes: %.2f secs, the score is: %.3f [initScore: %.3f]\n", 
		dataScan->scanID, elapsed, outMatch.alignScore, outMatch.initScore);

	return alignError;
}




//**************************************************************************************//
//                             Register All Specified Scans 
//**************************************************************************************//

void Register::AlignAllScans(int startScanID, int endScanID)
{
	clock_t beginTime = clock();

	/////////////////////////////////////////////////////////////
	// 1. Describe each scan using local voxelizer

	printf("==============Describing Scans==============\n");

	DescribeScans(startScanID, endScanID);
	SetFixedScan( scanList[startScanID] );

	clock_t endTime1 = clock();


	/////////////////////////////////////////////////////////////
	// 2. Register each scan (relative to the reference scan)

	printf("==============Matching Scans==============\n");

	// TODO: may need back-tracking here
	for (int i=startScanID+1; i<=_MIN((int)scanList.size()-1, endScanID); i++)
	{
		float alignScore = AlignScan( scanList[i] );

		UpdateAlignedScans( scanList[i] );

		printf("i=%d  score=%.2f \n", i, alignScore);
	}

	clock_t endTime2 = clock();


	/////////////////////////////////////////////////////////////
	// 3. Report the timing statistics

	float descTime  = ((float) (endTime1 - beginTime)) / CLOCKS_PER_SEC; // Time unit: second
	float matchTime = ((float) (endTime2 - endTime1))  / CLOCKS_PER_SEC; // Time unit: second
	float totalTime = ((float) (endTime2 - beginTime)) / CLOCKS_PER_SEC; // Time unit: second
	printf("The total time for registration: %.2f secs. \n", totalTime);
	printf("=> Scan description time: %.2f secs. \n", descTime);
	printf("=> Scan matching time:    %.2f secs. \n", matchTime);
}


void Register::DescribeScans(int startScanID, int endScanID)
{
	for (int i=startScanID; i<=_MIN((int)scanList.size()-1, endScanID); i++)
	{
		clock_t beginTime = clock();
		printf("computing descriptors for scan %2d ... ", i);

		scanList[i]->SampleScan();
		scanList[i]->DescribeScan();

		clock_t endTime = clock();
		float elapsed = ((float) (endTime - beginTime)) / CLOCKS_PER_SEC; // Time unit: second
		printf("taken %.2f secs \n", elapsed);
	}
}

void Register::SetFixedScan(Scan *fixedScan)
{
	kdTreeMaxDist = KDTREE_MAX_DIST_TIMES * fixedScan->bbox.GetDiagonalDist();

	UpdateAlignedScans( fixedScan );
}

void Register::UpdateAlignedScans(Scan *dataScan)
{
	//dataScan->isAligned = true;

	// TODO: remove sample points that are too close to each other (with similar descriptors)
	alignedScans.push_back( dataScan );

	vector<Point> scanSampPts = dataScan->sampler.GetSamplePoints();
	for ( int i=0; i<scanSampPts.size(); i++ )
	{
		Vector3f outPtPos;
		MultiplyVector(scanSampPts[i].pos, dataScan->alignLocalMat, outPtPos);

		alignedSampPts.push_back( outPtPos );
	}
}




//**************************************************************************************//
//                        Align Data Scan to All Registered Scans
//**************************************************************************************//

float Register::AlignScan(Scan *dataScan)
{
	clock_t beginTime = clock();

	//////////////////////////////////////////////////////////////
	// 1. Compute local voxelizers for feature points on the data and match them with that of the model

	float maxAlignScore = MIN_FLOAT;
	Match outMatch;

#ifdef USE_DESCRIPTOR_SGC
	for (int i = 0; i<dataScan->featSGCDescriptors.size(); i++)
#else
	for (int i=0; i<dataScan->featDescriptors.size(); i++)
#endif
	{
		// Find matched descriptor candidates
		vector<Match> matchList = GetMatchedDescriptors(dataScan->scanID, i);

		if (matchList.size() == 0)
		{
			continue;
		}
		// Select best matched descriptor
		Match bestMatch = SelectBestMatch(matchList, dataScan, i, false);


		// Update best alignment if current match is better 
		if ( bestMatch.alignScore > maxAlignScore )
		{
			maxAlignScore = bestMatch.alignScore;
			outMatch      = bestMatch;

			if ( bestMatch.alignScore > 0.2f )
			{
				printf("i=%d align: %.3f  desc: %.3f [initScore: %.3f]\n", i, bestMatch.alignScore, bestMatch.descDist, bestMatch.initScore);
			}

			// Exit if the alignment score is high enough
			if ( maxAlignScore > SCAN_ALIGN_SCORE_THRES )
			{
				cout << "engough high" << endl;
				break;
			}
		}


	}


	if ( maxAlignScore == MIN_FLOAT )
	{
		return maxAlignScore;
	}


	//////////////////////////////////////////////////////////////
	// 2. Align data to the model using the best match

	printf("start aligning the scan...\n");

	//EvaluateMatch( outMatch, false ); // Note: here only refined alignment with the target scan using ICP
	PerformAlignUsingMatch( outMatch, true ); // Note: here only refine alignment with the target scan using ICP
	//RefineAlign_ICP( dataScan );      // Note: refine data scan alignment with all the aligned scans using ICP

	clock_t endTime = clock();
	float elapsed = ((float) (endTime - beginTime)) / CLOCKS_PER_SEC; // Time unit: second
	printf("=> Matching scan %2d takes: %.2f secs, the score is: %.3f [initScore: %.3f]\n", 
		    dataScan->scanID, elapsed, outMatch.alignScore, outMatch.initScore);

	return outMatch.alignScore;
}

vector<Match> Register::GetMatchedDescriptors(int dataScanID, int dataFeatPtID)
{
	vector<Match> matchList;
	float minDist = MAX_FLOAT;

#ifdef USE_DESCRIPTOR_SGC

	SGCentroid* dataFeatDesc = scanList[dataScanID]->featSGCDescriptors[dataFeatPtID];

	// TODO: there are too many aligned descriptors for comparion which could increase the
	//       the chance of matching error; need to resolve this issue later
	for (int i = alignedScans.size() - 1; i >= 0; i--)
	{
		for (int j = 0; j<alignedScans[i]->sampSGCDescriptors.size(); j++)
		{
			SGCentroid* tagtSampDesc = alignedScans[i]->sampSGCDescriptors[j];
			float descDist = CompareSGCDescriptors(dataFeatDesc, tagtSampDesc);

		//	cout << endl<<"our dist is " << descDist << " ";
			// Note: need to normalize descriptor values to make this threshold consistent
			if (descDist > DESCRIPTOR_DIST_THRES)
				continue;

			Match match;
			match.dataScanID = dataScanID;
			match.dataFeatPtID = dataFeatPtID;

			match.tagtScanID = alignedScans[i]->scanID;
			match.tagtSampPtID = j;
			match.descDist = descDist;

			ReplaceSmallerMatchInList(match, matchList, DESCRIPTOR_TOP_MATCH_NUM);
		}
	}

#else

	vector<float> dataFeatDesc = scanList[dataScanID]->featDescriptors[dataFeatPtID]->descVector;

	// TODO: there are too many aligned descriptors for comparion which could increase the
	//       the chance of matching error; need to resolve this issue later
	for (int i=alignedScans.size()-1; i>=0; i--)
	{
		for (int j=0; j<alignedScans[i]->sampDescriptors.size(); j++)
		{
			vector<float> tagtSampDesc = alignedScans[i]->sampDescriptors[j]->descVector;
			float descDist = CompareDescriptors( dataFeatDesc, tagtSampDesc );

			// Note: need to normalize descriptor values to make this threshold consistent
			if ( descDist > DESCRIPTOR_DIST_THRES )
				continue;

			Match match;
			match.dataScanID   = dataScanID;
			match.dataFeatPtID = dataFeatPtID;

			match.tagtScanID   = alignedScans[i]->scanID;
			match.tagtSampPtID = j;
			match.descDist     = descDist;

			ReplaceSmallerMatchInList(match, matchList, DESCRIPTOR_TOP_MATCH_NUM);
		}
	}

#endif
	if (matchList.size() == 0)
	{
		return matchList;
	}

	SortDescriptorMatches( matchList );
	return matchList;
}

float Register::CompareSGCDescriptors(SGCentroid* modelDesc, SGCentroid* dataDesc)
{
	// Some descriptors could be empty
	if (modelDesc->FeatureDetails.size() == 0 || dataDesc->FeatureDetails.size() == 0)
	{
		return MAX_FLOAT;
	}

	return modelDesc->getDistanceWithOthers(*dataDesc);

}


float Register::CompareDescriptors(vector<float> modelDesc, vector<float> dataDesc)
{
	// Some descriptors could be empty
	if ( modelDesc.size() == 0 || dataDesc.size() == 0 )
	{
		return MAX_FLOAT;
	}

	/////////////////////////////////////////////////////////////
	// 1. Compute distance vector for all voxels

	vector<float> diffDescVec;
	for (int i=0; i<modelDesc.size(); i++)
	{
		float value = modelDesc[i] - dataDesc[i];
		diffDescVec.push_back( value );
	}


	/////////////////////////////////////////////////////////////
	// 2. Compute Euclidean distance for the two descriptors 

	float tempSum = 0;
	for (int i=0; i<diffDescVec.size(); i++)
	{
		tempSum = tempSum + (diffDescVec[i]*diffDescVec[i]); 
	}

	float descDist = sqrt( tempSum );
	//printf(" descDist: %.2f \n", descDist);

	return descDist;
}

void Register::ReplaceSmallerMatchInList(Match currMatch, vector<Match> &matchList, const int maxValueNum)
{
	if ( matchList.size() < maxValueNum )
	{
		matchList.push_back( currMatch );
	}

	else
	{
		for (int i=0; i<matchList.size(); i++)
		{
			if ( currMatch.descDist < matchList[i].descDist )
			{
				matchList[i] = currMatch;
				return;
			}
		}
	}
}

void Register::SortDescriptorMatches(vector<Match> &matchList)
{
	if ( matchList.size() == 0 )
		return;

	// Sort the distances of the matched descriptors 
	vector<float> descDistList;
	for (int i=0; i<matchList.size(); i++)
	{
		float descDist = matchList[i].descDist;
		descDistList.push_back( descDist ); 
	}
	vector<int> sortedIndices = BubbleSort(descDistList, true);		

	// Sort the match orders based on the descriptor distances
	vector<Match> origMatches = matchList;
	matchList.clear();
	for (int i=0; i<sortedIndices.size(); i++)
	{
		int index = sortedIndices[i];
		matchList.push_back( origMatches[index] );
	}
}

Match Register::SelectBestMatch(vector<Match> matchList, Scan *dataScan, int dataFeatPtID, bool useInitScore)
{
	Match bestMatch;
	if ( matchList.size() == 0 )
	{
		printf("Warning: There is no descriptor match candidate. \n");
		return bestMatch;
	}
	else if ( matchList.size() == 1 )
	{
		EvaluateMatch( matchList[0], true );
		bestMatch = matchList[0];
		return bestMatch;
	}
	else
	{
		EvaluateMatch( matchList[0], true );
		bestMatch = matchList[0];

		for (int i=1; i<matchList.size(); i++)
		{
			EvaluateMatch( matchList[i], true );

			// Select best match based on the alignment score using LRF only
			if ( useInitScore )
			{
				if ( matchList[i].initScore > bestMatch.initScore)
				{
					bestMatch = matchList[i];
				}
			}

			// Select best match based on the alignment score after using ICP 
			else
			{
				if ( matchList[i].alignScore > bestMatch.alignScore)
				{
					bestMatch = matchList[i];
				}
			}

		}
			
		return bestMatch;
	}
}

AlignError Register::PerformAlignUsingMatch(Match match, bool useICP)
{
	AlignError alignError;
	
	int dataScanID   = match.dataScanID;
	int dataFeatPtID = match.dataFeatPtID;
	int tagtScanID   = match.tagtScanID;
	int tagtSampPtID = match.tagtSampPtID;
	float descDist   = match.descDist;

	Scan *dataScan = scanList[dataScanID];

	////////////////////////////////////////////////////////////////
	// 1. Compute the scan alignment matrix using two LRFs

	ComputeScanAlignMatrix(scanList[tagtScanID], tagtSampPtID, dataScan, dataFeatPtID);

	double beforeICPAlignMat[16];
	memcpy(beforeICPAlignMat, dataScan->alignLocalMat, sizeof(double)*16);


	////////////////////////////////////////////////////////////////
	// 2. Refine the scan alignment matrix using ICP

	RefineAlign_ICP(scanList[tagtScanID], dataScan);

	double afterICPAlignMat[16];
	memcpy(afterICPAlignMat, dataScan->alignLocalMat, sizeof(double)*16);


	////////////////////////////////////////////////////////////////
	// 3. Compare the scan alignment matrices before and after using ICP

	ScanAlign scanAlign;
	float rotError, tranError;
	scanAlign.CompareAlignMatrix(beforeICPAlignMat, afterICPAlignMat, rotError, tranError);


	////////////////////////////////////////////////////////////////
	// 4. Use the alignment matrix with applying ICP refinement

	if ( useICP == false )
	{
		memcpy(dataScan->alignLocalMat, beforeICPAlignMat, sizeof(double)*16);
	}


	////////////////////////////////////////////////////////////////
	// 5. Return the error of alignment matrix using LRFs

	alignError.match     = match;
	alignError.rotError  = rotError;
	alignError.tranError = tranError;

	return alignError;
}




//**************************************************************************************//
//                             Evaluate Descriptor Match
//**************************************************************************************//

void Register::EvaluateMatch(Match &match, bool isResetAlignMat)
{
	int dataScanID   = match.dataScanID;
	int dataFeatPtID = match.dataFeatPtID;
	int tagtScanID   = match.tagtScanID;
	int tagtSampPtID = match.tagtSampPtID;
	float descDist   = match.descDist;

	Scan *dataScan = scanList[dataScanID];

	/////////////////////////////////////////////////////
	// 1. Keep the align matrix of the data scan

	double origAlignMat[16];
	memcpy(origAlignMat, dataScan->alignLocalMat, sizeof(double)*16);


	/////////////////////////////////////////////////////
	// 2. Evaluate the align score before and after ICP refinment
	ComputeScanAlignMatrix(scanList[tagtScanID], tagtSampPtID, dataScan, dataFeatPtID);
	float initAlignScore = EvaluateAlign( dataScan, false );
	float newAlignScore = 0.0;

	if ( initAlignScore > INIT_ALIGN_SCORE_THRES )
	{

		// Refine the align matrix of the data scan and evaluate it
		RefineAlign_ICP(scanList[tagtScanID], dataScan);
		//RefineAlign_ICP( dataScan ); // Note: this function is very slow
		newAlignScore = EvaluateAlign( dataScan, false );

	}
	else
	{

		newAlignScore = initAlignScore;
	}

	match.initScore  = initAlignScore;
	match.alignScore = newAlignScore;

	/////////////////////////////////////////////////////
	// 3. Reset the align matrix of the data scan

	memcpy(dataScan->alignLocalMat, origAlignMat, sizeof(double)*16);
}

void Register::ComputeScanAlignMatrix(Scan *tagtScan, int tagtSampPtID, Scan *dataScan, int dataFeatPtID)
{
	// Compute align matrix based on the original local mesh
	ScanAlign alignScan;
	alignScan.InitScanAlign(tagtScan->sampSGCDescriptors[tagtSampPtID]->keyPtLRF,
		dataScan->featSGCDescriptors[dataFeatPtID]->keyPtLRF);
	alignScan.ComputeAlignMatrix();

	// Obtain the align matrix that aligns data scan to the transformed model scan
	double modelAlignMat[16];
	memcpy(modelAlignMat, tagtScan->alignLocalMat, sizeof(double)*16);
	modelAlignMat[12] = tagtScan->alignLocalMat[12] - ALIGNED_SCAN_SHIFT_VECTOR[0];
	modelAlignMat[13] = tagtScan->alignLocalMat[13] - ALIGNED_SCAN_SHIFT_VECTOR[1];
	modelAlignMat[14] = tagtScan->alignLocalMat[14] - ALIGNED_SCAN_SHIFT_VECTOR[2];

	MultiplyMatrix(modelAlignMat, alignScan.alignMatrix, alignScan.alignMatrix);

	// Shift the aligned scene for rendering
	memcpy(dataScan->alignLocalMat,  alignScan.alignMatrix, sizeof(double)*16);
	dataScan->alignLocalMat[12] = alignScan.alignMatrix[12] + ALIGNED_SCAN_SHIFT_VECTOR[0];
	dataScan->alignLocalMat[13] = alignScan.alignMatrix[13] + ALIGNED_SCAN_SHIFT_VECTOR[1];
	dataScan->alignLocalMat[14] = alignScan.alignMatrix[14] + ALIGNED_SCAN_SHIFT_VECTOR[2];
}

float Register::EvaluateAlign(Scan *dataScan, bool isPrint)
{
	/////////////////////////////////////////////////////////////////
	// 1. Get uniform sample points on model + data

	vector<Point> origDataSampPts  = dataScan->sampler.GetSamplePoints();

	if ( alignedSampPts.size() == 0 || origDataSampPts.size() == 0 )
	{
		printf("Warning: the uniform sample points on model or data is empty. \n");
		return MAX_FLOAT;
	}


	/////////////////////////////////////////////////////////////////
	// 2. Save sample points on the data scan

	currDataSamplePts.clear();
	for (int i=0; i<origDataSampPts.size(); i++)
	{
		Vector3f outPt;
		MultiplyVector(origDataSampPts[i].pos, dataScan->alignLocalMat, outPt);
		currDataSamplePts.push_back( outPt );
	}


	/////////////////////////////////////////////////////////////////
	// 3. For each data sample points, find the closest sample point on the model

	float alignScore = ComputeAverageDistance(alignedSampPts, kdTreeMaxDist, currDataSamplePts, currModelClosePts, isPrint);

	return alignScore;
}

float Register::EvaluateAlign(Scan *tagtScan, Scan *dataScan)
{
	kdTreeMaxDist = KDTREE_MAX_DIST_TIMES * tagtScan->bbox.GetDiagonalDist();

	/////////////////////////////////////////////////////////////////
	// 1. Get uniform sample points on model + data

	vector<Point> tagtSampPts = tagtScan->sampler.GetSamplePoints();
	vector<Point> dataSampPts = dataScan->sampler.GetSamplePoints();

	if ( tagtSampPts.size() == 0 || dataSampPts.size() == 0 )
	{
		printf("Warning: the uniform sample points on model or data is empty. \n");
		return MAX_FLOAT;
	}


	/////////////////////////////////////////////////////////////////
	// 2. Save sample points on the model

	vector<Vector3f> currTagtSamplePts;
	for (int i=0; i<tagtSampPts.size(); i++)
	{
		Vector3f outPt;
		MultiplyVector(tagtSampPts[i].pos, tagtScan->alignLocalMat, outPt);

		currTagtSamplePts.push_back( outPt );
	}


	/////////////////////////////////////////////////////////////////
	// 3. Save sample points on the data

	currDataSamplePts.clear();
	for (int i=0; i<dataSampPts.size(); i++)
	{
		Vector3f outPt;
		MultiplyVector(dataSampPts[i].pos, dataScan->alignLocalMat, outPt);

		currDataSamplePts.push_back( outPt );
	}


	/////////////////////////////////////////////////////////////////
	// 4. For each data sample points, find the closest sample point on the model

	float alignScore = ComputeAverageDistance(currTagtSamplePts, kdTreeMaxDist, currDataSamplePts, currModelClosePts, false);

	return alignScore;
}

float Register::ComputeAverageDistance(vector<Vector3f> kdTreePointSet, float kdTreeMaxDist, vector<Vector3f> dataPointSet, vector<Vector3f> &kdTreeClosePts, bool isPrint)
{
	kdTreeClosePts.clear();

	/////////////////////////////////////////////////////////////////
	// 1. For each data sample points, find the closest sample point on the model

	KDtree kdTree( kdTreePointSet );
	//const float kdTreeMaxDist = 3.0 * scanList[modelID]->bbox.GetDiagonalDist();

	vector<float> distList;
	for (int i=0; i<dataPointSet.size(); i++)
	{
		const float *tempPt = kdTree.closest_to_pt(dataPointSet[i], kdTreeMaxDist, NULL);

		Vector3f closestPt;
		closestPt[0] = tempPt[0];
		closestPt[1] = tempPt[1];
		closestPt[2] = tempPt[2];

		kdTreeClosePts.push_back( closestPt );
		float dist = len(closestPt - dataPointSet[i]);
		distList.push_back( dist );
	}

	//printf("[%d %d] \n", currDataSamplePts.size(), currModelClosePts.size());


	/////////////////////////////////////////////////////////////////
	// 2. Compute the corresponding point between two point sets (using a closest distance threshold)

	float avgPtDist = 0;
	for (int i=0; i<distList.size(); i++)
	{
		avgPtDist += distList[i];
	}
	avgPtDist = avgPtDist / (float)distList.size();

	float closeDistThres = SAMPLE_POINT_DIST;
	int closePtNum = 0;
	for (int i=0; i<distList.size(); i++)
	{
		if ( distList[i] < closeDistThres )
		{
			closePtNum++;
		}
	}
	

	/////////////////////////////////////////////////////////////////
	// 3. Compute the score (i.e., overlap) as the ratio between correponding
	//    point and the minimum number of scan sample points

	int scanMinPtNum = _MIN(kdTreePointSet.size(), dataPointSet.size());
	float closePtRatio = closePtNum / (float)scanMinPtNum;
	//float closePtRatio = closePtNum/(float)distList.size();

	if ( isPrint )
	{
		printf("avgDist: %.3f  ratio: %.3f\n", avgPtDist, closePtRatio );
	}

	//printf("avgDist: %.3f  Input:[%d %d] Output:[%d %d] ratio: %.3f\n", 
	//	   avgPtDist, kdTreePointSet.size(), dataPointSet.size(), distList.size(), closePtNum, closePtRatio );

	return closePtRatio;
}




//**************************************************************************************//
//                               Refine Scan Alignment
//**************************************************************************************//

void Register::RefineAlign_ICP(Scan *tagtScan, Scan *dataScan)
{
	xform modelXForm;
	xform dataXForm;

	modelXForm = xform( tagtScan->alignLocalMat  );
	dataXForm  = xform( dataScan->alignLocalMat );

	ICP(tagtScan->GetTriMesh(), dataScan->GetTriMesh(), modelXForm, dataXForm, 0, false, false);

	memcpy(dataScan->alignLocalMat, dataXForm, sizeof(double)*16);
}

void Register::RefineAlign_ICP(Scan *dataScan)
{
	if ( alignedScans.size() == 0 )
		return;

	xform modelXForm;
	xform dataXForm;
	modelXForm = xform( alignedScans[0]->alignLocalMat  ); 
	dataXForm  = xform( dataScan->alignLocalMat );

	TriMesh *modelMesh = new TriMesh();
	for (int i=0; i<alignedScans.size(); i++)
	{
		// Avoid translating the aligned scan twice
		double alignLocalMat[16];
		memcpy(alignLocalMat, alignedScans[i]->alignLocalMat, sizeof(double)*16);
		alignLocalMat[12] = alignedScans[i]->alignLocalMat[12] - ALIGNED_SCAN_SHIFT_VECTOR[0];
		alignLocalMat[13] = alignedScans[i]->alignLocalMat[13] - ALIGNED_SCAN_SHIFT_VECTOR[1];
		alignLocalMat[14] = alignedScans[i]->alignLocalMat[14] - ALIGNED_SCAN_SHIFT_VECTOR[2];

		vector<Vector3f> tempVertices = alignedScans[i]->GetTriMesh()->vertices;

		for (int j=0; j<tempVertices.size(); j++)
		{
			Vector3f outPtPos;
			//MultiplyVector(tempVertices[j], alignedScans[i]->alignLocalMat, outPtPos);
			MultiplyVector(tempVertices[j], alignLocalMat, outPtPos);

			//outPtPos = outPtPos - ALIGNED_SCAN_SHIFT_VECTOR;
			modelMesh->vertices.push_back( outPtPos );
		}
	}

	ICP(modelMesh, dataScan->GetTriMesh(), modelXForm, dataXForm, 0, false, false);

	memcpy(dataScan->alignLocalMat, dataXForm, sizeof(double)*16);

	delete modelMesh;
}




//**************************************************************************************//
//                              Functions for Debug
//**************************************************************************************//

void Register::Function_HighMatch(Scan *dataScan, Scan *tagtScan)
{
	// specify a key point, find highest matches on the other scans
	tagtScan->SampleScan();
	tagtScan->DescribeScan();

	FindHighMaches(dataScan, tagtScan);
}

void Register::FindHighMaches(Scan *dataScan, Scan *tagtScan)
{
	vector<float> dataFeatDesc = dataScan->voxelizer.descVector;
	if ( dataFeatDesc.size() == 0 )
		return;

	vector<double> distList;
	float minDist = MAX_FLOAT;
	for (int i=0; i<tagtScan->sampDescriptors.size(); i++)
	{
		vector<float> tagtSampDesc = tagtScan->sampDescriptors[i]->descVector;

		float descDist = CompareDescriptors( dataFeatDesc, tagtSampDesc );
		distList.push_back( descDist );
	}

	sortedIndices = BubbleSort(distList, true);
	int tagtSampPtID = sortedIndices[0];

	vector<Point> tagtSampPts = tagtScan->sampler.GetSamplePoints();
	tagtScan->BuildLocalVoxelizer(tagtSampPts[tagtSampPtID], LOCAL_SHAPE_RADIUS, true);

	ComputeScanAlignMatrix(tagtScan, dataScan);
}

void Register::ComputeScanAlignMatrix(Scan *tagtScan, Scan *dataScan)
{
	// Compute align matrix based on the original local mesh
	ScanAlign alignScan;
	alignScan.InitScanAlign(tagtScan->voxelizer.keyPtLRF, dataScan->voxelizer.keyPtLRF);
	alignScan.ComputeAlignMatrix();

	// Obtain the align matrix that aligns data scan to the transformed model scan
	double modelAlignMat[16];
	memcpy(modelAlignMat, tagtScan->alignLocalMat, sizeof(double)*16);
	modelAlignMat[12] = tagtScan->alignLocalMat[12] - ALIGNED_SCAN_SHIFT_VECTOR[0];
	modelAlignMat[13] = tagtScan->alignLocalMat[13] - ALIGNED_SCAN_SHIFT_VECTOR[1];
	modelAlignMat[14] = tagtScan->alignLocalMat[14] - ALIGNED_SCAN_SHIFT_VECTOR[2];

	MultiplyMatrix(modelAlignMat, alignScan.alignMatrix, alignScan.alignMatrix);

	// Shift the aligned scene for rendering
	memcpy(dataScan->alignLocalMat,  alignScan.alignMatrix, sizeof(double)*16);
	dataScan->alignLocalMat[12] = alignScan.alignMatrix[12] + ALIGNED_SCAN_SHIFT_VECTOR[0];
	dataScan->alignLocalMat[13] = alignScan.alignMatrix[13] + ALIGNED_SCAN_SHIFT_VECTOR[1];
	dataScan->alignLocalMat[14] = alignScan.alignMatrix[14] + ALIGNED_SCAN_SHIFT_VECTOR[2];
}




//**************************************************************************************//
//                                  Drawing Stuff
//**************************************************************************************//

void Register::DrawClosestPoint()
{
	glDisable( GL_LIGHTING );
	glColor3f(0.3, 0.1, 0.9);
	glPointSize( 4.0 );

	glColor3f(0.3, 0.1, 0.9);
	for (int i=0; i<currDataSamplePts.size(); i++)
	{
		glBegin(GL_POINTS);
		glVertex3f(currDataSamplePts[i][0], currDataSamplePts[i][1], currDataSamplePts[i][2]);
		glEnd();
	}

	glColor3f(0.8, 0.6, 0.1);
	for (int i=0; i<currModelClosePts.size(); i++)
	{
		glBegin(GL_POINTS);
		glVertex3f(currModelClosePts[i][0], currModelClosePts[i][1], currModelClosePts[i][2]);
		glEnd();
	}

	glColor3f(0.8, 0.1, 0.1);
	for (int i=0; i<currModelClosePts.size(); i++)
	{
		glBegin(GL_LINES);
		glVertex3f(currModelClosePts[i][0], currModelClosePts[i][1], currModelClosePts[i][2]);
		glVertex3f(currDataSamplePts[i][0], currDataSamplePts[i][1], currDataSamplePts[i][2]);
		glEnd();
	}

	glPointSize( 1.0 );
	glEnable(GL_LIGHTING);
}

void Register::DrawHighMatches(Scan *tagtScan)
{
	vector<Point> scanSamplePts = tagtScan->sampler.GetSamplePoints();

	if ( sortedIndices.size() == 0 || scanSamplePts.size() == 0 )
		return;

	glDisable( GL_LIGHTING );
	glColor3f(0.1, 0.1, 0.1);
	glPointSize( 8.0 );

	// Draw mesh sample points with high match value
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadMatrixd( tagtScan->globalMatrix );
	for (int i=0; i<HIGH_MATCH_NUM; i++)
	{
		if ( i == 0 )         glColor3f(0.9, 0.1, 0.1);
		else if ( i == 1 )    glColor3f(0.1, 0.9, 0.1);
		else if ( i == 2 )    glColor3f(0.1, 0.1, 0.9);
		else if ( i == 3 )    glColor3f(0.9, 0.9, 0.1);
		else if ( i == 4 )    glColor3f(0.9, 0.1, 0.9);
		else                  glColor3f(0.1, 0.1, 0.1);         

		int index = sortedIndices[i];

		glBegin(GL_POINTS);
		glVertex3f(scanSamplePts[index].pos[0], scanSamplePts[index].pos[1], scanSamplePts[index].pos[2]);
		glEnd();
	}
	glPopMatrix();

	glPointSize( 1.0 );
	glEnable(GL_LIGHTING);
}