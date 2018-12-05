///////////////////////////////////////////////////////////////
//
// Model.h
//
//   Reconstructed 3D Model Class
//
// by Song Peng ( songpeng@ustc.edu.cn )
// 
// 11/May/2015
//
///////////////////////////////////////////////////////////////

#ifndef _MODEL_H
#define _MODEL_H

#include <vector>
#include "TriMesh.h"

using namespace std;
using namespace trimesh;


class Model 
{
private:
	TriMesh *trimesh;

public:
	bool isColored;                        // If the mesh has vertex color
	vector<Point> verList;                 // Vertex positions (with normals) 
	vector<Triangle*> triList;             // Triangles (vertex positions)
	_Box bbox;                             // Bounding box of the model

	float mtlAmbient[4];                   // Model rendering material
	float mtlDiffuse[4];
	float mtlSpecular[4];
	float mtlEmission[4];

	//double localMatrix[16];                // Initial mesh local matrix  (OpenGL style matrix) 
	//double globalMatrix[16];               // Initial mesh global matrix (OpenGL style matrix)


public:
	Model();
	~Model();
	void ClearScan();
	//void InitPose();
	//void ResetPose();
	TriMesh* GetTriMesh();

	// Load Mesh Model
	void LoadMesh(const char *fileName);
	void ProcessMesh(const Vector3f diffuseMTL, const float scaleFactor);
	void ComputeBoundingBox();
	void SetDiffuseMTL(Vector3f diffuseMTL);
	void NormalizeModel(const float scale);
	bool GetModelVertices(vector<Point> &vertices);
	bool GetModelTriangles(vector<Triangle*> &triangles);

	// Draw Mesh Model
	void DrawModel(bool isAlign);
	void DrawTStrips(const TriMesh *themesh);
	void DrawModelWire(Vector3f color, float width);
	void DrawModelBBox(Vector3f color, float width);
};

#endif