///////////////////////////////////////////////////////////////
//
// LVoxel.h
//
//   Local Volume Voxel Class
//
// by Song Peng ( songpeng@ustc.edu.cn )
// 
// 11/May/2015
//
///////////////////////////////////////////////////////////////

#ifndef _L_VOXEL_H
#define _L_VOXEL_H

#include "vec.h"
#include <vector>
#include "HelpStructs.h"

using namespace std;

class LVoxel 
{
public:
	Vector3f center;               // Voxel center position
	Vector3f size;                 // Voxel size
	Vector3i pos;                  // Discrete position in the volume
	int state;                     // Voxel [IN]/[ON]/[OUT] model 

private:
	vector<Triangle*> origTris;   // Triangles inside/cross the voxel (coordiantes in LRF)
	vector<Triangle*> clipTris;    // Triangles inside the voxel after box-mesh clipping (coordiantes in LRF)
	vector<SampPoint> sampPts;     // Sampled points inside partial voxel


public:
	LVoxel();
	~LVoxel();
	void ClearClipTris();
	vector<Triangle*>  GetVoxelOrigTris();

	// Compute Voxel Triangles
	void AddTriangle(Triangle* _tri);
	void ClipLocalMesh();
	void SetSamplePoints(vector<SampPoint> _sampPts);

	// Compute Voxel Attributes
	float ComputeSurfaceArea();
	Vector3f ComputeSurfaceCenter();
	Vector3f ComputeSurfaceNormal();
	float ComputeSurfaceCurvature();
	float ComputeShapeVolume();
	Vector3f ComputeShapeCenter();
	Vector3f ComputeSurfaceColor();

	// Draw Voxels
	void DrawVoxel(float lineWidth, Vector3f color);
	void DrawSamplePoints(float pointSize, Vector3f color);
	void DrawVoxelTris();
	void DrawVoxelTrisWire(Vector3f color, float width);
	void DrawTris(vector<Triangle*> triList);
	void DrawTrisWire(vector<Triangle*> triList, Vector3f color, float width);
	void DrawDebug();

	int getSampledPointsNum();
	bool isPointIn(Vector3f Point);
};

#endif
