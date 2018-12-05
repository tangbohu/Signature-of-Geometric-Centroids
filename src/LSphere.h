///////////////////////////////////////////////////////////////
//
// LSphere.h
//
//   Build Local Sphere for Describing Local Shape
//
// by Song Peng ( songpeng@ustc.edu.cn )
// 
// 14/May/2015
//
///////////////////////////////////////////////////////////////

#ifndef _L_SPHERE_H
#define _L_SPHERE_H

#include "vec.h"
#include <vector>
#include "HelpStructs.h"

using namespace std;

class LSBin;

class LSphere 
{
public:
	Vector3i gridSize;
	float radius;
	vector<LSBin> binGrid;

	vector<Triangle*> shapeTris;   // Transformed mesh triangles in the sphere (radius = R)
	vector<Point> shapePts;     // Transformed mesh vertices in the sphere (radius = R)


public:
	LSphere();
	~LSphere();
	void ClearLSphere();
	void InitLSphere(Vector3i _gridSize, float _radius, vector<Triangle*> _LRFTris, vector<Point> _shapePts);

	// Preprocess Local Shape
	void ComputPointDensity(float localRadius);
	vector<Point> GetNeighborPoints(Vector3f center, float radius, vector<Point> pointList);

	// Compute Bin Grid
	void ComputeBinGrid();	
	void ComputBinState();

	// Compute Descriptor (3DSC and SHOT)
	void ComputeDescriptor_3DSC(vector<float> &descVec);
	void ComputeDescriptor_SHOT(vector<float> &descVec, Vector3f keyPtNormal);

	// Draw volume
	void DrawLocalSphere(float lineWidth, Vector3f color);
	void DrawSpheBinGrid(float lineWidth, vec color);
	void DrawLocalShapeWire(vec color, float width);
};

#endif
