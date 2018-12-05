///////////////////////////////////////////////////////////////
//
// LPointGrid.h
//
//   3D Local Point Grid Class
//
// by Song Peng ( songpeng@ustc.edu.cn )
// 
// 01/Apr/2015
//
///////////////////////////////////////////////////////////////

#ifndef _L_POINT_GRID_H
#define _L_POINT_GRID_H

#include "HelpStructs.h"
#include <vector>

struct Point;
struct HitPoint;

using namespace std;

class LPointGrid 
{
public:
	Vector3f bboxMinPt;                   // Origin of the point grid
	Vector3f gridSize;                    // Distance between neighboring points
	Vector3i gridDimen;                   // [(n-1)*W+1  (n-1)*H+1  (n-1)*D+1]

private:
	int pointNum;                         // Sampling point number: (W+1)*(H+1)*(D+1)
	SampPoint *pointList;                 // Sampling point array (size: (W+1)*(H+1)*(D+1))  
	vector<vector<HitPoint>> XRayHitPts;  // HitPoint array for Y-axis aligned rays


public:
	LPointGrid();
	~LPointGrid();
	void InitLPointGrid(Vector3f _bboxMinPt, Vector3f voxelSize, Vector3i volumeSize);
	void SetXAxisHitPoints(vector<vector<HitPoint>> _XRayHitPts);

	// Compute Point Grid
	void ComputePointGrid();
	void ComputePointState();
	vector<SampPoint> GetPointsInsideVoxel(Vector3i voxelPos);

	// Utility Functions
	Vector3f GetPointPosition(Vector3i pointCoord);
	int GetPointIndex(Vector3i pointCoord);
	vector<HitPoint> GetRayHitPoints(Vector3i pointCoord, Vector3i rayDir);

	// Draw Point Grid
	void DrawPointGrid(float pointSize, Vector3f pointColor);
	void DrawCrossPoints();
};

#endif