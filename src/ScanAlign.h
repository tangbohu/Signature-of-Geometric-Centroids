///////////////////////////////////////////////////////////////
//
// ScanAlign.h
//
//   Aligning Two Scan Meshes Class
//
// by Song Peng ( songpeng@ustc.edu.cn )
// 
// 11/May/2015
//
///////////////////////////////////////////////////////////////

#ifndef _SCAN_ALIGN_H
#define _SCAN_ALIGN_H

#include "HelpStructs.h"

using namespace std;

class LRF;


class ScanAlign
{
public:
	LRF modelLRF;                   // LRF of model scan
	LRF dataLRF;                    // LRF of data scan
	double alignMatrix[16];         // Align matrix for aligning data to model

public:
	ScanAlign();
	~ScanAlign();
	void InitScanAlign(LRF _modelLRF, LRF _dataLRF);

	// Compute Alignment Matrix
	void ComputeAlignMatrix();

	// Compare Alignment Matrix
	void CompareAlignMatrix(double alignMat[16], double realMat[16], float &rotError, float &tranError);
};

#endif