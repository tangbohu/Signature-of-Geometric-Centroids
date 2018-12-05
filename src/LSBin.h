///////////////////////////////////////////////////////////////
//
// LSBin.h
//
//   Local Spherical Bin Class
//
// by Song Peng ( songpeng@ustc.edu.cn )
// 
// 14/May/2015
//
///////////////////////////////////////////////////////////////

#ifndef _LS_BIN_H
#define _LS_BIN_H

#include "vec.h"
#include <vector>
#include "HelpStructs.h"

using namespace std;

class LSBin 
{
public:
	float minRadius;             // Minimum radial distance   ( Vector3f[0] )
	float maxRadius;             // Maximum radial distance
	float minTheta;              // Minimum polar angle       ( Vector3f[1]; range: [0, 180] )
	float maxTheta;              // Maximum polar angle
	float minAlpha;              // Minimum azimuthal angle   ( Vector3f[2]; range: [0, 360) )
	float maxAlpha;              // Maximum azimuthal angle

	Vector3i pos;                // Discrete position in the sphere
	int state;                   // Bin [IN]/[ON]/[OUT] model 
	Vector3f cen;                // Average position of bin's 8 corner points
	float volume;                // Bin volume size (relative to the sphere volume)

	vector<Triangle*> binTris;   // Triangles inside the bin
	vector<Point> binPoints;     // Vertices inside the bin


public:
	LSBin();
	~LSBin();
	void InitBin(float _minRadius, float _maxRadius, float _minTheta, float _maxTheta, float _minAlpha, float _maxAlpha);
	void PrintBin();

	// Compute Bin Variables
	void ComputeBinCenter();
	float ComputeBin3DSCValue();
	vector<int> ComputeBinSHOTValues(Vector3f keyPtNormal);

	// Construct Spherical Bin
	void GetBinTriangles(vector<Triangle*> localTris);
	void GetBinPoints(vector<Point> localPts);
	//void GetBinPoints(vector<Triangle*> localTris);
	//void UpdateBinPoints(vector<Point> &binPtList);

	// Utility Functions
	bool IsTriangleInBin(Vector3f vertices[3]);
	bool IsPointInBin(Vector3f point);
	Vector3f GetCartesianCoord(float radius, float theta, float alpha);
	Vector3f GetSphericalCoord(Vector3f point);
	int GetPointIndexInList(Vector3f tagtPt, vector<Vector3f> ptList);

	// Draw Shape in Bin
	void DrawBinTriangles(Vector3f color);
	void DrawBinPoints(Vector3f color);

	// Draw Bin
	void DrawSphericalBin(int drawStep, Vector3f color);
	void DrawBinCenter(Vector3f color, int pointSize);
	void DrawSphericalCurve_Alpha(float radius, float theta, float alphaStart,  float alphaEnd, float alphaStep);
	void DrawSphericalCurve_Theta(float radius, float alpha, float thetaStart,  float thetaEnd, float thetaStep);
	void DrawSphericalLine_Radius(float theta,  float alpha, float radiusStart, float radiusEnd);


};

#endif
