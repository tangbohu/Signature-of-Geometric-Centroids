#ifndef _SGC_DES_H
#define _SGC_DES_H

#include <vector>
#include "HelpStructs.h"
#include "SGCBin.h"
#include "LRF.h"
using namespace std;

class SGCBin;

class SGCentroid
{
public:

	vector<SGCBin> FeatureDetails;
	LRF keyPtLRF;                  // LRF defined by the support of the key point


public:
	SGCentroid();
	~SGCentroid();

	void ClearSGCentroid();
	void InsertBin(SGCBin bin);
	double getDistanceWithOthers(SGCentroid& that);
	//void InitSGCentroid(Vector3i _gridSize, float _radius, vector<Triangle*> _LRFTris, vector<Point> _shapePts);

	//// Preprocess Local Shape
	//vector<Point> GetNeighborPoints(Vector3f center, float radius, vector<Point> pointList);

	//// Compute Descriptor 
	//void ComputeDescriptor_SGC(vector<float> &descVec);

	//// Draw volume
	//void DrawLocalSphere(float lineWidth, Vector3f color);
	//void DrawSpheBinGrid(float lineWidth, vec color);
	//void DrawLocalShapeWire(vec color, float width);
};

#endif
