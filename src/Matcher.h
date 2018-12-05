///////////////////////////////////////////////////////////////
//
// Matcher.h
//
//   Match Model and Data Meshes
//
// by Song Peng ( songpeng@ustc.edu.cn )
// 
// 11/May/2015
//
///////////////////////////////////////////////////////////////

#ifndef _MATCHER_H
#define _MATCHER_H

#include "HelpStructs.h"
#include <vector>
#include "Register.h"
#include "Model.h"

using namespace std;


class Matcher 
{
public:
	vector<Scan*>  scanList;         // Triangular mesh of each scan
	float scaleFactor;               // Scale factor for all scanList (i.e., normalize 1st scan with in [-1.0, 1.0])
	double worldMatrix[16];          // World matrix
	double worldAxesMatrix[16];      // World axes matrix

	Register myRegister;             // Reconstruct a 3D model from the input scans
	vector<Point> modelPoints;       // Summarize vertices (with normals) for aligned scans
	Model model;                     // Load and draw the reconstructed 3D model

	int tagtScanID;                  // Debug variable: for testing pairwise scan alignment


public:
	Matcher();
	~Matcher();
	void ResetScene();
	void ClearScans();

	// World and Scan Matrix
	void InitWorldMatrix(Vector3f position, Vector3f rotAxis, float rotAngle);
	void InitWorldAxesMatrix(Vector3f position, Vector3f rotAxis, float rotAngle);
	void UpdateGlobalMatrix();
	void GetGlobalMatrix(double *localMat, double *globalMat);

	// Load Scan Meshes
	void LoadScans(char *scanFilePath);
	void SetScanModelsMTL();
	Vector3f GenerateRandomMTL(Vector3f existColors[10], int colorNum);
	float GetMinColorDistance(Vector3f existColors[10], int colorNum, Vector3f randColor);

	// Get Scan File Names
	vector<string> GetScanFileNames(char *scanFilePath);
	int GetScanFileNumber(string scanFileName);
	void SortScanFileNames(vector<string> &scanFileNames, bool isPrint);
	string GetFirstFileName(const char scanFilePath[]);

	// Reconstruct 3D Model from Scans
	void AlignAllScans(int startScanID, int endScanID);
	Point PickScanSurfPoint(int scanID, int winX, int winY);
	void ComputeScanDescriptor(int scanID, Point keyPoint, float radius);

	// Poisson Reconstruction
	void PoissonReconstruction(char *modelFileName, int PSRLevel);
	void SummarizeScanPoints();
	string GetNPTSFileName(const char modelFileName[]);

	// Functions for Debug
	void Function_Sample(int startScanID, int endScanID);
	void Function_Evaluate(int startScanID, int endScanID);
	void Function_Test(int startScanID, int endScanID);
	void Function_Match(int dataScanID, int endScanID);
	void Function_Experiment(int startScanID, int endScanID);
	void Function_MeshClipping();
	void Function_MeshTriNum();
	void Function_PairwiseAlign();
	void Function_RefineAlign();

	// Read and Write Files
	void WriteAlignMatrices(const char *fileName);
	void ReadAlignMatrices(const char *fileName);
	void WriteModelPoints(const char *fileName);
	void ReadModel(const char *fileName);

	// Draw Scans
	void DrawScans(int startScanID, int endScanID, int mode);
	void DrawScansWire(int startScanID, int endScanID);
	void DrawScansBBox(int startScanID, int endScanID);
	void DrawScansFeaturePoints(int startScanID, int endScanID);
	void DrawScansSamplePoints(int startScanID, int endScanID);
	void DrawAlignedScans(int startScanID, int endScanID);

	// Draw Voxelizer 
	void DrawScansRay(int startScanID, int endScanID);
	void DrawScansKeyPoint(int startScanID, int endScanID);
	void DrawScansLRFShape(int startScanID, int endScanID);
	void DrawScansLRFSphere(int startScanID, int endScanID);
	void DrawScansLRF(int startScanID, int endScanID);
	void DrawScansLocalShape(int startScanID, int endScanID);
	void DrawScansLocalShapeWire(int startScanID, int endScanID);
	void DrawScansLocalBBox(int startScanID, int endScanID);
	void DrawScansLocalGrid(int startScanID, int endScanID);
	void DrawTargetScanHighMatches();

	void DrawScansLocalSphere(int startScanID, int endScanID);


	// Draw Scan Alignment
	void DrawClosestPoint();
	void ShowPickedScan(int pickScanID);
	void DrawModelPoints();
	
	// Draw 3D Model
	void DrawModel();
	void DrawModelWire();
	void DrawModelBBox();

	// Draw World Axes and Planes
	void DrawWorldAxes(int winW, int winH, float currFovy);
	void RotateWorldAxes(vec rotAxis, float rotAngle);

	// Transform Scan
	void TranslateScan(int scanID, Vector3f transVec);
	void RotateScan(int scanID, Vector3f rotAxis, float rotAngle);
	void ScaleWorld_AroundScan(int scanID, Vector3f scaleVec);

	// Transform World
	void TranslateWorld(Vector3f transVec);
	void RotateWorld(Vector3f rotAxis, float rotAngle);
	void ScaleWorld(Vector3f scaleVec);
};

#endif
