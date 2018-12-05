///////////////////////////////////////////////////////////////
//
// ScanSampler.h
//
//   Sampling Scan Mesh Vertices Class
//
// by Song Peng ( songpeng@ustc.edu.cn )
// 
// 11/May/2015
//
///////////////////////////////////////////////////////////////

#ifndef _SCAN_SAMPLER_H
#define _SCAN_SAMPLER_H

#include <vector>
#include "HelpStructs.h"

using namespace std;

class ScanSampler 
{
private:
	vector<Point> meshFeaturePts;      // Sample pints on a mesh for matching scans
	vector<Point> meshSamplePts;       // Sample pints on a mesh for describing a scan
	vector<int> featPtIndices;         // The original point index of each feature point


public:
	ScanSampler();
	~ScanSampler();
	vector<Point> GetFeaturePoints();
	vector<Point> GetSamplePoints();
	vector<int> GetFeaturePointIndices();

	// Mesh Vertex Sampling (uniform sampling)
	void SampleMesh_Feature(vector<Point> verList);
	void SampleMesh_Uniform(vector<Point> verList);

	// Sample Point From a List
	vector<Point> GetRandomPointList(vector<Point> pointList, int sampleNum, const float distThres, vector<int> &randPtIndices);
	int GetRandomObjIndex(int objNum);
	bool IsValidSample(vector<Point> currDataSamplePts, Point currPt, const float distThres);

	// Draw Mesh Sampling
	void DrawSamplePoints();
	void DrawFeaturePoints();
	void DrawPointList(vector<Point> pintList, int pointSize, Vector3f pointColor);
};

#endif