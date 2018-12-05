///////////////////////////////////////////////////////////////
//
// LRF.h
//
//   Local Reference Frame Class
//
// by Song Peng ( songpeng@ustc.edu.cn )
// 
// 27/Mar/2015
//
///////////////////////////////////////////////////////////////

#ifndef _LRF_H
#define _LRF_H

#include <vector>
#include "HelpStructs.h"

using namespace std;


class LRF 
{
public:
	Point keyPoint;               // Key point to define local shape 
	float radius;                 // Radius to define local shape
	vector<Triangle*> LRFTris;    // Local triangles that define LRF (using radius R)
	vector<Point> LRFPoints;    // Local triangles that define LRF (using radius R)

	vec localAxes[3];             // LRF axes       
	double LRFMatrix[16];         // LRF matrix (from LRF to mesh local space; saved in opengl-style)

public:
	LRF();
	~LRF();
	LRF & operator=(const LRF &lrf);
	void InitLRF(Point _keyPoint, float _radius, vector<Triangle*> _LRFTris, vector<Point> _LRFpoints);

	// Build Local Reference Frame (LRF)
	void BuildLocalFrame(bool onlyPoint= false);
	void ComputeLRFMatrix();


	// Compute Scatter Matrix for the point
	void ComputeScatterMatrix(vector<Point> LRFpoints, Vector3f keyPointPos, float radius, double matrix[9]);

	// Compute Scatter Matrix for the Triangles
	void ComputeScatterMatrix(vector<Triangle*> LRFTris, Vector3f keyPointPos, float radius, double matrix[9]);
	float ComputeLocalShapeArea(vector<Triangle*> LRFTris);
	float GetAreaWeight(Triangle tri, float totalArea);
	float GetDistWeight(Triangle tri, Vector3f keyPointPos, float radius);
	void GetScatterMatrix(Triangle tri, Vector3f keyPointPos, double matrix[9]);
	void ComputeMatrix(vec P, vec Q, double mat[9]);

	// Compute Local Frame Axes
	void ComputeLocalFrameAxes(double matrix[9], float radius);
	void GetLocalAxisSign(vector<Triangle*> LRFTris, vector<Point> LRFPoints, Vector3f keyPointPos, float radius, vec &localAxis);
	void GetLocalAxisSign_Normal(vector<Triangle*> LRFTris, vector<Point> LRFPoints, Point keyPoint, float radius, vec &localAxis);

	// Draw LRF
	void DrawKeyPoint(Vector3f color);
	void DrawLRFSphere(Vector3f color);
	void DrawLRF();
	void DrawLRFShape();
};

#endif