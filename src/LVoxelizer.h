///////////////////////////////////////////////////////////////
//
// LVoxelizer.h
//
//   Local Voxelizer Class
//
// by Song Peng ( songpeng@ustc.edu.cn )
// 
// 13/May/2015
//
///////////////////////////////////////////////////////////////

#ifndef _L_VOXELIZER_H
#define _L_VOXELIZER_H

#include "HelpStructs.h"
#include "LRF.h"
#include "LVoxel.h"
#include "LVolume.h"
#include "LSBin.h"
#include "LSphere.h"
#include "SGCentroid.h"

using namespace std;

struct Descriptor
{
	LRF keyPtLRF;                  // LRF defined by the support of the key point
	vector<float> descVector;      // Local shape descriptor vector (could be empty, e.g., 
	                               // when there is no local triangles around an isolated sample point)
};

class LVoxelizer 
{
public:
	Point keyPoint;                // Key point to define local shape 
	float radius;                  // Radius to define local shape
	vector<Triangle*> LRFTris;     // Local triangles that define LRF (using radius R; no transformation)
	LRF keyPtLRF;                  // LRF defined by the support of the key point

	vector<Triangle*> shapeTris;   // Local triangles that define volume (using radius sqrt(2)*R; no transformation)
	LVolume localVolume;           // Local shape volume
	vector<float> descVector;      // Local shape descriptor vector

	LSphere localSphere;           // Local shape sphere
	vector<Point> shapePts;        // TODO: the names of shapePts and shapeTris are confusing; change it later
	vector<Point> LRFPts;        // TODO: the names of shapePts and shapeTris are confusing; change it later


public:
	LVoxelizer();
	~LVoxelizer();
	void InitLVoxelizer(Point _keyPoint, float _radius, vector<Triangle*> _LRFTris, vector<Triangle*> _shapeTris, vector<Point> _LRFPts, vector<Point> _shapePts);
	
	// Compute Descriptor using LSphere
	Descriptor* ComputeDescriptor_Spherical(int featureTypeID, bool isDebug);
	vector<Point> TransformPoints(vector<Point> _shapePoints);

	// Compute Descriptor using LVolume
	Descriptor* ComputeDescriptor(int featureTypeID, bool isDebug);
	SGCentroid * ComputeSGC();

	vector<Triangle*> TransformTriangles(vector<Triangle*> _LRFTris);

};

#endif