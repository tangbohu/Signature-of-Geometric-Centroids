///////////////////////////////////////////////////////////////
//
// Experiment.h
//
//   Experiment to Find Best Parameters for Local Voxelizer
//
// by Song Peng ( songpeng@ustc.edu.cn )
// 
// 02/May/2015
//
///////////////////////////////////////////////////////////////

#ifndef _EXPERIMENT_H
#define _EXPERIMENT_H

#include "HelpStructs.h"
#include <vector>

using namespace std;

struct RPDot
{
	float one_precision;
	float recall;
};


class Experiment 
{
public:
	Scan* sceneScan;
	Scan* modelScan;

	vector<int> sortedIndices;              // Temp: mark sample points with highest scores

	vector<RPDot> RPDotList;

	Register helpRegister;                  // For utilizing some functions in Register class


public:
	Experiment();
	~Experiment();
	void InitScanPair(Scan *_sceneScan, Scan *_modelScan);

	void MatchScanPair();
	void FindHighMaches(float distRatioThres, RPDot &RPPoint);
	bool FindModelHighMaches(int sceneSampPtID, float distRatioThres, int &modelSampPtID);
	bool ValidateMatch(int sceneSampPtID, int modelSampPtID, float pointDistThres);
	void GetMatchedDescriptor(vector<float> dataFeatDesc, bool isSpeedUp, int &tagtScanID, int &tagtSampPtID);
	
	void WriteRPCurveFile(const char *fileName);
};

#endif
