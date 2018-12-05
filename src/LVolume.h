///////////////////////////////////////////////////////////////
//
// LVolume.h
//
//   Build Local Volume for Describing Local Shape
//
// by Song Peng ( songpeng@ustc.edu.cn )
// 
// 11/May/2015
//
///////////////////////////////////////////////////////////////

#ifndef _L_VOLUME_H
#define _L_VOLUME_H

#include "vec.h"
#include <vector>
#include "HelpStructs.h"
#include "SGCentroid.h"

using namespace std;

class LVoxel;
class LPointGrid;
struct HitPoint;

class LVolume 
{
public:
	Vector3i volumeSize;            // Volume size [W H D]   
	Vector3f volumeDimen;           // Volume dimension [X Y Z]

	Vector3f volumeMinPt;           // Volume minimum point (w.r.t LRF)
	Vector3f volumeMaxPt;           // Volume maximum point (w.r.t LRF)
	Vector3f volumeCenPt;           // Volume center point  (w.r.t LRF)
	Vector3f voxelSize;             // Voxel size [x y z] = [X/W Y/H Z/D]
	vector<LVoxel> voxelGrid;       // Voxel array (size: W*H*D)

	vector<Triangle*> sphereTris;   // Transformed mesh triangles in the sphere (radius = sqrt(2)*R)
	vector<Triangle*> volumeTris;   // Transformed mesh triangles in the cubical volume (edge = 2*R)
	vector<Point> volumePoints;

	LPointGrid *pointGrid;          // Sample point grid for the volume


public:
	LVolume();
	~LVolume();
	void ClearLVolume();
	void InitLVolume(Vector3i _volumeSize, Vector3f _volumeDimen, vector<Triangle*> _sphereTris);
	void InitLVolumeWithTransformPoints(Vector3i _volumeSize, Vector3f _volumeDimen, vector<Point> _volumePoints);


	// Compute Voxel Grid
	void ComputeVoxelGrid();
	void ComputeVoxelState();
	int MapToGrid(float value, float min, float max, int N);
	int GetVoxelIndex(Vector3i voxelPos);

	// Voxel Local Mesh Clipping 
	void ClipMeshInVoxels();

	// Compute Point Grid
	void ComputePointGrid();
	void ComputeCrossPoints(vector<vector<HitPoint>> &XRayHitPts);
	void UpdateTestTriangles(Vector3i voxelPos, vector<Triangle*> &testTris);


	void ComputeSGC(SGCentroid& feature);

	// Compute Descriptor (using Single Feature)
	void ComputeDescriptor_SurfArea(vector<float> &descVec);
	void ComputeDescriptor_SurfPos(vector<float> &descVec);
	void ComputeDescriptor_Normal(vector<float> &descVec, Vector3f keyPtNormal);
	void ComputeDescriptor_Curvature(vector<float> &descVec);
	void ComputeDescriptor_ShapeVolume(vector<float> &descVec);
	void ComputeDescriptor_ShapePos(vector<float> &descVec);
	void ComputeDescriptor_Color(vector<float> &descVec);

	// Compute Descriptor (using Feature Combinations)
	void ComputeDescriptor_Area_SurfCen(vector<float> &descVec);
	void ComputeDescriptor_Area_Normal(vector<float> &descVec, Vector3f keyPtNormal);
	void ComputeDescriptor_Volume_VolCen(vector<float> &descVec);
	void ComputeDescriptor_SurfCen_Normal(vector<float> &descVec, Vector3f keyPtNormal);
	void ComputeDescriptor_Area_SurfCen_Normal(vector<float> &descVec, Vector3f keyPtNormal);

	// Draw volume
	void DrawLocalBBox(float lineWidth, Vector3f color);
	void DrawVoxelGrid(float lineWidth, Vector3f color);
	void DrawLocalShape();
	void DrawLocalShapeWire(vec color, float width);
	
	// Draw Point Grid
	void DrawPointGrid(float pointSize, Vector3f pointColor);
	void DrawCrossPoints();
	void DrawVolumeSamplePts(int pointSize, Vector3f color);
};

#endif
