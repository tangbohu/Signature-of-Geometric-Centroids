///////////////////////////////////////////////////////////////
//
// HelpStructs.h
//
//   Common Structures
//
// by Song Peng ( songpeng@ustc.edu.cn )
// 
// 01/Apr/2015
//
///////////////////////////////////////////////////////////////


#ifndef _HELP_STRUCTS_H
#define _HELP_STRUCTS_H

#include "vec.h"
using namespace trimesh;


///////////////////////////////////////////////////////////////
// Define Structs
///////////////////////////////////////////////////////////////

// 3D Sampling Point
struct SampPoint
{
	Vector3f pos;               // Point position 
	Vector3i coord;			    // Discrete coordinates
	int state;                  // Point [IN]/[OUT] model
};


// 3D Point (with normal) 
struct Point
{
	Vector3f pos;     // Point position
	Vector3f nor;     // Point normal

	Vector3f color;   // Point color (if available)
	float curv;       // Point curvature

	Point & operator=(const Point &pt)
	{
		if( this == &pt )
			return *this;

		this->pos   = pt.pos;
		this->nor   = pt.nor;

		this->color = pt.color;
		this->curv  = pt.curv;

		return *this;
	};
};


// 3D Triangle
struct Triangle
{
	//int id;              // For mesh deformation
	Vector3f v[3];       // Vertex positions 

	Vector3f color[3];   // Vertex colors
	float curv[3];       // Vertex mean curvatures

	Vector3f normal;     // Face normal
	float area;          // Face area

	bool isCross;        // Is triangle cross the cube (i.e., voxel)

	//bool isMark;         // Mark a triangle for debug
	//bool isTwoUpper;     // If there is two upper vertices along the plane normal
	//Vector3f tempPts[2]; // Two intersecting points between the triangle and a plane

	void ComputeNormal()
	{
		Vector3f tempNor = (v[1] - v[0]) CROSS (v[2] - v[0]);  // Assume the vertices are saved in counter-clockwise
		float tempNorLen = len(tempNor);

		if ( tempNorLen != 0 )    normal = tempNor / tempNorLen;
		else                      normal = Vector3f(1,0,0);     // Note: this default vector also can be others
	};

	void ComputeArea()
	{
		Vector3f normal  = (v[1] - v[0]) CROSS (v[2] - v[0]);
		area  = len(normal);
	};

	Vector3f ComputeCenter()
	{
		Vector3f center = (v[0]+v[1]+v[2])/3.0f;
		return center;
	};


	void CorrectNormal(Vector3f tagtNormal)
	{
		// Compute initial normal
		ComputeNormal();

		// Rearrange vertex order if needed
		float dotp = normal DOT tagtNormal;
		if ( dotp < 0 )
		{
			Vector3f triVers[3];
			for (int i=0; i<3; i++)
			{
				triVers[i] = v[i];
			}

			v[0] = triVers[0];
			v[1] = triVers[2];
			v[2] = triVers[1];
		}

		// Recompute the normal
		ComputeNormal();
	}
};


// 3D Box (avoid ambiguity with Box in Trimesh )
struct _Box 
{
	Vector3f minPt;
	Vector3f maxPt;
	Vector3f cenPt;

	_Box() :
	minPt(0, 0, 0),
	maxPt(0, 0, 0),
	cenPt(0, 0, 0)
	{}

	void Transform(vec transVec, float scale)
	{
		minPt += transVec;
		maxPt += transVec;
		cenPt += transVec;

		minPt *= scale;
		maxPt *= scale;
		cenPt *= scale;
	};

	_Box & operator=(const _Box &box)
	{
		if( this == &box )
			return *this;

		this->minPt = box.minPt;
		this->maxPt = box.maxPt;
		this->cenPt = box.cenPt;

		return *this;
	};

	float GetDiagonalDist()
	{
		return len(maxPt - minPt);
	};

	void PrintBox()
	{
		printf("Box: [%7.3f  %7.3f  %7.3f]      [%7.3f  %7.3f  %7.3f] \n", minPt[0], minPt[1], minPt[2], maxPt[0], maxPt[1], maxPt[2]);
	};
};


// 3D Circle
struct Circle
{
	Vector3f center;
	Vector3f normal;
	float radius;
};


#endif