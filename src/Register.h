///////////////////////////////////////////////////////////////
//
// Register.h
//
//   Class for Registering Scans 
//
// by Song Peng ( songpeng@ustc.edu.cn )
// 
// 11/May/2015
//
///////////////////////////////////////////////////////////////

#ifndef _REGISTER_H
#define _REGISTER_H

#include "HelpStructs.h"
#include <vector>

using namespace std;

struct Match
{
	int dataScanID;        // ScanID of the data scan
	int dataFeatPtID;      // Descriptor ID of feature point of the data scan

	int tagtScanID;        // ScanID that data scan is to be aligned
	int tagtSampPtID;      // Descriptor ID that current descriptor is matched
	float descDist;        // Distance between two descriptor vectors   

	float alignScore;      // Alignment score after ICP refinement
	float initScore;       // Alignment score before ICP refinement

	Match & operator=(const Match &match)
	{
		if( this == &match )
			return *this;

		this->dataScanID   = match.dataScanID;
		this->dataFeatPtID = match.dataFeatPtID;

		this->tagtScanID   = match.tagtScanID;
		this->tagtSampPtID = match.tagtSampPtID;
		this->descDist     = match.descDist;

		this->alignScore   = match.alignScore;
		this->initScore    = match.initScore;

		return *this;
	};
};

struct ScanOverlap
{
	int tagtScanID;
	int dataScanID;
	int isAligned;

	float realOverlap;
	float compOverlap; 

	float rotError;        // Rotational error of alignment matrix using LRFs
	float tranError;       // Translational error of alignment matrix using LRFs
};

struct AlignError 
{
	Match match;           // Scan alignment corresponding to a match
	float rotError;        // Rotational error of alignment matrix using LRFs
	float tranError;       // Translational error of alignment matrix using LRFs
};


class Register 
{
public:  
	vector<Scan*>  scanList;                // Triangular mesh of each scan

	vector<Scan*> alignedScans;             // All the scans aligned (relative to a reference scan)
	vector<Vector3f> alignedSampPts;        // All the aligned scan sample points
    float kdTreeMaxDist;                    // Max distance for building kd-tree 

	vector<Vector3f> currDataSamplePts;     // Sample points on the data scan (after multiplying the scan local matrix)
	vector<Vector3f> currModelClosePts;     // Closest point on the model scan to corresponding data point
	vector<int> sortedIndices;              // Temp: mark sample points with highest scores


public:
	Register();
	~Register();
	void InitRegister(vector<Scan*> _scanList);

	// Experiment: Scan Pair Overlap
	void Experiment_Overlap(int startScanID, int endScanID);
	void WriteScanOverlapFile(const char *fileName, vector<ScanOverlap> overlapList);
	void WriteScanOverlapFile_Simple(const char *fileName, vector<ScanOverlap> overlapList);

	// Experiment: Alignment Accuracy
	void Experiment_Accuracy(int startScanID, int endScanID);
	void WriteAlignErrorFile(const char *fileName, vector<AlignError> errorList);

	// Align Data Scan to a Target Scan
	void AlignScanPair();
	AlignError AlignScanWithTarget(Scan *dataScan, Scan *tagtScan, bool useICP);
	void UpdateLVoxelizers(Match match);

	// Register All Specified Scans 
	void AlignAllScans(int startScanID, int endScanID);
	void DescribeScans(int startScanID, int endScanID);
	void SetFixedScan(Scan *fixedScan);
	void UpdateAlignedScans(Scan *dataScan);

	// Align Data Scan to All Registered Scans
	float AlignScan(Scan *dataScan);
	vector<Match> GetMatchedDescriptors(int dataScanID, int dataFeatPtID);
	float CompareDescriptors(vector<float> modelDesc, vector<float> dataDesc);
	float CompareSGCDescriptors(SGCentroid* modelDesc, SGCentroid* dataDesc);

	void ReplaceSmallerMatchInList(Match currMatch, vector<Match> &matchList, const int maxValueNum);
	void SortDescriptorMatches(vector<Match> &matchList);		
	Match SelectBestMatch(vector<Match> matchList, Scan *dataScan, int dataFeatPtID, bool useInitScore);
	AlignError PerformAlignUsingMatch(Match match, bool useICP);

	// Evaluate Descriptor Match
	void EvaluateMatch(Match &match, bool isResetAlignMat);
	void ComputeScanAlignMatrix(Scan *tagtScan, int tagtSampPtID, Scan *dataScan, int dataFeatPtID);
	float EvaluateAlign(Scan *dataScan, bool isPrint);
	float EvaluateAlign(Scan *tagtScan, Scan *dataScan);
	float ComputeAverageDistance(vector<Vector3f> kdTreePointSet, float kdTreeMaxDist, vector<Vector3f> dataPointSet, vector<Vector3f> &kdTreeClosePts, bool isPrint);

	// Refine Scan Alignment
	void RefineAlign_ICP(Scan *tagtScan, Scan *dataScan);
	void RefineAlign_ICP(Scan *dataScan);

	// Function for Debug
	void Function_HighMatch(Scan *dataScan, Scan *tagtScan);
	void FindHighMaches(Scan *dataScan, Scan *tagtScan);
	void ComputeScanAlignMatrix(Scan *tagtScan, Scan *dataScan);

	// Drawing Stuff
	void DrawClosestPoint();
	void DrawHighMatches(Scan *tagtScan);
};

#endif
