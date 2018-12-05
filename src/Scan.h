///////////////////////////////////////////////////////////////
//
// Scan.h
//
//   General Triangular Scan Mesh Class
//
// by Song Peng ( songpeng@ustc.edu.cn )
// 
// 11/May/2015
//
///////////////////////////////////////////////////////////////

#ifndef _SCAN_H
#define _SCAN_H

#include <vector>
#include "TriMesh.h"
#include "LVoxelizer.h"
#include "SGCentroid.h"
#include "ScanSampler.h"

using namespace std;
using namespace trimesh;

class Scan 
{
private:
	TriMesh *trimesh;

public:
	bool isColored;                        // If the mesh has vertex color
	vector<Point> verList;                 // Vertex positions (with normals) 
	vector<Triangle*> triList;             // Triangles (vertex positions)
	_Box bbox;                             // Bounding box of the model
	int scanID;                            // Scan ID according to the saving index in scanList

	float mtlAmbient[4];                   // Model rendering material
	float mtlDiffuse[4];
	float mtlSpecular[4];
	float mtlEmission[4];

	double localMatrix[16];                // Initial mesh local matrix  (OpenGL style matrix)
	double globalMatrix[16];               // Initial mesh global matrix (OpenGL style matrix)
	double alignLocalMat[16];              // Aligned mesh local matrix  (OpenGL style matrix)
	double alignGlobalMat[16];             // Aligned mesh global matrix (OpenGL style matrix)

	ScanSampler sampler;                   // Sample mesh vertices
	vector<Descriptor*> featDescriptors;   // Descriptors of sampled feature points
	vector<Descriptor*> sampDescriptors;   // Descriptors of sampled sample points 
	vector<SGCentroid*> sampSGCDescriptors;
	vector<SGCentroid*> featSGCDescriptors;   // Descriptors of sampled feature points

	Vector3f pickPoint;                    // 3D picking point on the target scan mesh       (for debug)
	Vector3f origPoint;                    // 3D OpenGL origin point in the world coordinate (for debug)
	Point intersetPt;                      // 3D intersected point on the target scan mesh   (for debug)
	LVoxelizer voxelizer;                  // Local voxelizer around the intersected point   (for debug)
	SGCentroid sgc;
	//bool isAligned;                        // If the scan aligned successfully


public:
	Scan();
	~Scan();
	void ClearScan();
	void InitPose();
	void ResetPose();
	TriMesh* GetTriMesh();

	// Load Mesh Model
	void LoadMesh(const char *fileName);
	void ProcessMesh(const Vector3f diffuseMTL, const float scaleFactor);
	void ComputeBoundingBox();
	void SetDiffuseMTL(Vector3f diffuseMTL);
	float GetScaleFactor();
	void NormalizeModel(const float scale);
	void CheckCurvatures();
	bool GetModelVertices(vector<Point> &vertices);
	bool GetModelTriangles(vector<Triangle*> &triangles);

	// Describe Scan with Descriptors
	void SampleScan();
	void DescribeScan();
	void ComputeSampleDescriptors();
	void ComputeFeatureDescriptors();
	Descriptor* BuildLocalVoxelizer(Point keyPoint, float radius, bool isDebug);
	SGCentroid* BuildSGC(Point keyPoint, float radius, bool isDebug);

	void UpdateLVoxelizer(Descriptor* descriptor);  // For debug
	Descriptor* BuildLocalSphere(Point keyPoint, float radius, bool isDebug); // For comparison

	// Get Mesh Local Shape
	void GetLocalShape_Tri(Point keyPoint, float radius, vector<Triangle*> &LRFTris, vector<Triangle*> &shapeTris);
	void GetLocalShape_Ver(Point keyPoint, float radius, vector<Point> &LRFPts, vector<Point> &shapePts);
	bool IsTriangleInsideSphere(Triangle tri, Vector3f center, float radius);
	bool IsVertexInsideSphere(Point ver, Vector3f center, float radius);

	// Pick Scan Mesh Point
	Point PickMeshSurfPoint(int winX, int winY);
	Vector3f GetOGLPos(int winX, int winY);

	// Draw Mesh Model
	void DrawScan(bool isAlign);
	void DrawTStrips(const TriMesh *themesh);
	void DrawPoints(const  TriMesh *themesh);

	void DrawScanWire(Vector3f color, float width);
	void DrawScanBBox();
	void DrawFeaturePoints(bool isAlign);
	void DrawSamplePoints(bool isAlign);
	void DrawEvaluatePoints(bool isAlign);

	// Draw Local Shape, LRF, Volume etc.
	void DrawRay();
	void DrawKeyPoint();
	void DrawLRFShape();
	void DrawLRFSphere();
	void DrawLRF();
	void DrawLocalShape();
	void DrawLocalShapeWire();
	void DrawLocalBBox();
	void DrawLocalGrid();

	void DrawLocalSphere();
};

#endif