///////////////////////////////////////////////////////////////
//
// LVoxel.cpp
//
//   Local Volume Voxel Class
//
// by Song Peng ( songpeng@ustc.edu.cn )
// 
// 11/May/2015
//
///////////////////////////////////////////////////////////////

#include "Controls.h"
#include "Helper.h"
#include "LVoxel.h"
#include "ClipMeshInBox.h"
#include "HelpDefines.h"
#include "BoxTriTest.h"

//**************************************************************************************//
//                                   Initialization
//**************************************************************************************//

LVoxel::LVoxel()
{
	center    =  Vector3f(0, 0, 0);
	size      =  Vector3f(1, 1, 1);
	pos       =  Vector3i(0, 0, 0);
	state     =  -1;
}

LVoxel::~LVoxel()
{
	ClearClipTris();
}

void LVoxel::ClearClipTris()
{
	for (int i=0; i<clipTris.size(); i++)
	{
		delete clipTris[i];
	}
	clipTris.clear();
}

vector<Triangle*> LVoxel::GetVoxelOrigTris()
{
	return origTris;
}




//**************************************************************************************//
//                              Compute Voxel Triangles
//**************************************************************************************//

void LVoxel::AddTriangle(Triangle* _tri)
{
	origTris.push_back( _tri );
}

void LVoxel::ClipLocalMesh()
{
	if ( origTris.size() == 0 )
		return;

	ClearClipTris();

	clipTris = ClipMeshInBox(origTris, center, size);
}

void LVoxel::SetSamplePoints(vector<SampPoint> _sampPts)
{
	sampPts = _sampPts;
}




//**************************************************************************************//
//                              Compute Voxel Attributes
//**************************************************************************************//

float LVoxel::ComputeSurfaceArea()
{
	if ( clipTris.size() == 0 ) 
		return DEFAULT_VOXEL_VALUE;

	// Compute local surface area within the voxel
	// TODO: need actual area for CROSS triangles (esp. for large triangles)
	float surfArea = 0;
	for (int i=0; i<clipTris.size(); i++)
	{
		surfArea += clipTris[i]->area; 
	}

	 // Normalize the area such its value is roughly in [0,1] (1.5 is a fixed parameter)
	surfArea = surfArea / (1.5f*size[0]*size[0]); 

	return surfArea;
}

Vector3f LVoxel::ComputeSurfaceCenter()
{
	if ( clipTris.size() == 0 )
		return Vector3f(0,0,0);

	// Compute local surface center using all triangles in the voxel
	Vector3f surfCenter = Vector3f(0,0,0);
	for (int i=0; i<clipTris.size(); i++)
	{
		Vector3f triCenter = clipTris[i]->ComputeCenter();
		surfCenter += triCenter;
	}
	surfCenter = surfCenter / (float)origTris.size();

	// Normalize the position such that its norm is within [0,1]
	surfCenter = (surfCenter - (center-0.5f*size)) / (SQUARE_ROOT_THREE*size[0]); 

	return surfCenter;
}

Vector3f LVoxel::ComputeSurfaceNormal()
{
	if ( clipTris.size() == 0 )
		return Vector3f(0,0,0);

	// Compute local surface normal within the voxel
	Vector3f surfNormal = Vector3f(0,0,0);
	for (int i=0; i<clipTris.size(); i++)
	{
		surfNormal += clipTris[i]->normal;
	}

	if ( clipTris.size() > 0 )
		surfNormal = surfNormal / (float)clipTris.size();

	//printf("surfNor: [%.2f %.2f %.2f] \n", surfNormal[0], surfNormal[1], surfNormal[2]);

	return surfNormal;
}

float LVoxel::ComputeSurfaceCurvature()
{	
	if ( origTris.size() == 0 )
	{
		return DEFAULT_VOXEL_VALUE;
	}

	int curvNum = 0;
	float avgCurvature = 0;

	for (int i=0; i<origTris.size(); i++)
	{
		// Note: curvature value could be negative; plus may not be a good way for accumulation
		for (int j=0; j<3; j++)
		{
			if ( IsPointInBox( origTris[i]->v[j], center, size) )
			{
				avgCurvature += fabs(origTris[i]->curv[j]);
				curvNum++;			
			}
		}
	}

	// Note: this curvature needs to be weighted when combining with other normalized features 
	//       by defining CURVATURE_WEIGHT (e.g., 0.02)
	if ( curvNum > 0 )
		avgCurvature = avgCurvature / (float)curvNum;

	return avgCurvature;
}

float LVoxel::ComputeShapeVolume()
{
	if ( origTris.size() == 0 || sampPts.size() == 0 )
		return DEFAULT_VOXEL_VALUE;

	// Compute the number of sampled points inside the mesh
	int inMeshPtNum = 0;
	for (int i=0; i<sampPts.size(); i++)
	{
		if ( sampPts[i].state == POINT_IN_MESH )
		{
			inMeshPtNum++;
		}

		if ( sampPts[i].state == POINT_UNKNOWN )
		{
			printf("Warning: There should not exist point with unknown state. \n");
		}
	}

	// Estimate local shape volume using sampled points (size range: [0,1])
	int totalPtNum  = sampPts.size();
	float volume = inMeshPtNum / (float)totalPtNum;   // Normalize the volume
	//printf("voxel volume: %.2f \n", volume);

	return volume;
}

Vector3f LVoxel::ComputeShapeCenter()
{
	if ( origTris.size() == 0 || sampPts.size() == 0 )
		return Vector3f(0,0,0);

	// Compute the number of sampled points inside the mesh
	Vector3f shapeCenter = Vector3f(0,0,0);
	int inMeshPtNum = 0;
	for (int i=0; i<sampPts.size(); i++)
	{
		if ( sampPts[i].state == POINT_IN_MESH )
		{
			shapeCenter += sampPts[i].pos;
			inMeshPtNum++;
		}

		if ( sampPts[i].state == POINT_UNKNOWN )
		{
			printf("Warning: There should not exist point with unknown state. \n");
		}
	}

	if ( inMeshPtNum == 0 )
	{
		return Vector3f(0,0,0);
	}

	// Estimate local shape volume using sampled points
	shapeCenter = shapeCenter / (float)inMeshPtNum;
	//printf("shape center [%.2f %.2f %.2f] \n", shapeCenter[0], shapeCenter[1], shapeCenter[2]);

	// Normalize the position such that its norm is within [0,1]
	shapeCenter = (shapeCenter - (center-0.5f*size)) / (SQUARE_ROOT_THREE*size[0]); // Normalize the position 
	return shapeCenter;
}

Vector3f LVoxel::ComputeSurfaceColor()
{
	if ( origTris.size() == 0 )
		return Vector3f(0,0,0);

	// Compute local surface center using all triangles in the voxel
	Vector3f surfColor = Vector3f(0,0,0);
	for (int i=0; i<origTris.size(); i++)
	{
		Vector3f triColor = (origTris[i]->color[0]+origTris[i]->color[1]+origTris[i]->color[2])/3.0f;
		surfColor += triColor;
	}
	surfColor = surfColor / (float)origTris.size();

	// Convert the 3D location into a 1D index number (using the LRF)
	//Vector3f tempPos = surfColor - (center-0.5f*size); 
	//float posIndex = tempPos[0]*size[1]*size[2] + tempPos[1]*size[2] + tempPos[2];
	//float colorIndex = surfColor[0] + surfColor[1]/255.0 + surfColor[2]/(255.0*255.0);

	//printf("color [%.2f %.2f %.2f] \n", surfColor[0], surfColor[1], surfColor[2]);

	return surfColor;
}




//**************************************************************************************//
//                                   Draw Voxels 
//**************************************************************************************//

void LVoxel::DrawVoxel(float lineWidth, Vector3f color)
{
	vec minPt = center - 0.5f*size;
	vec maxPt = center + 0.5f*size;

	DrawWireCuboid(minPt, maxPt, lineWidth, color);
}

void LVoxel::DrawSamplePoints(float pointSize, Vector3f color)
{
	glDisable( GL_LIGHTING );
	glColor3f(color[0], color[1], color[2]);
	glPointSize( pointSize );

	glBegin(GL_POINTS);
	for (int i=0; i<sampPts.size(); i++)
	{
		if ( sampPts[i].state == POINT_IN_MESH )
			glColor3f(0.6, 0.1, 0.8);
		else if ( sampPts[i].state == POINT_OUT_MESH ) 
			glColor3f(0.8, 0.6, 0.1);
		else
			glColor3f(0.1, 0.1, 0.1);

		if ( sampPts[i].state == POINT_IN_MESH )
			glVertex3f(sampPts[i].pos[0], sampPts[i].pos[1], sampPts[i].pos[2]);
	}
	glEnd();

	glPointSize( 1.0 );
	glEnable(GL_LIGHTING);
}


void LVoxel::DrawVoxelTris()
{
	//DrawTris( origTris );
	DrawTris( clipTris );
}

void LVoxel::DrawVoxelTrisWire(Vector3f color, float width)
{
	//DrawTrisWire(origTris, color, width);
	DrawTrisWire( clipTris, color, width);
}

void LVoxel::DrawTris(vector<Triangle*> triList)
{
	if ( triList.size() == 0 )
		return;

	glEnable(GL_LIGHTING);

	for (int i=0; i<triList.size(); i++)
	{
		Triangle *tri = triList[i];

		glBegin(GL_TRIANGLES);
		glNormal3fv( tri->normal );
		glVertex3fv( tri->v[0] );
		glVertex3fv( tri->v[1] );
		glVertex3fv( tri->v[2] );
		glEnd();
	}
}

void LVoxel::DrawTrisWire(vector<Triangle*> triList, Vector3f color, float width)
{
	if ( triList.size() == 0 )
		return;

	glDisable( GL_LIGHTING );
	glColor3f(color[0], color[1], color[2]);
	glLineWidth( width );

	// Set the offset value
	const float epsilon = 1e-3;
	float projMat[16];
	glGetFloatv(GL_PROJECTION_MATRIX, projMat);

	// Offset the drawing by modifying the projection matrix
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glTranslatef(0.0f, 0.0f, -epsilon);
	glMultMatrixf(projMat);
	for (int i=0; i<triList.size(); i++)
	{
		Triangle *tri = triList[i];

		glBegin(GL_LINE_LOOP);
		glVertex3fv( tri->v[0] );
		glVertex3fv( tri->v[1] );
		glVertex3fv( tri->v[2] );
		glEnd();
	}
	glPopMatrix();

	glMatrixMode(GL_MODELVIEW);
	glLineWidth( 1.0 );
	glEnable(GL_LIGHTING);
}

void LVoxel::DrawDebug()
{
	if ( origTris.size() == 0 )
		return;

	glDisable(GL_LIGHTING);
	glPointSize( 10.0 );
	glLineWidth( 5.0 );
	glColor3f( 0.9, 0.1, 0.1 );

	Vector3f planeNor = Vector3f(0,1,0); 
	Vector3f planePt  = center + 0.5f*size;
	Vector3f endPt    = planePt + 0.05f*planeNor;

	glColor3f( 0.1, 0.9, 0.1 );	
	glBegin(GL_POINTS );
	glVertex3f(planePt[0], planePt[1], planePt[2]);
	glEnd();

	glColor3f( 0.6, 0.9, 0.1 );	
	glBegin(GL_LINES );
	glVertex3f(planePt[0], planePt[1], planePt[2]);
	glVertex3f(endPt[0],    endPt[1],   endPt[2]);
	glEnd();

	//glLineWidth( 3.0 );

	//for (int i=0; i<tris.size(); i++)
	//{
	//	Triangle *tri = tris[i];

	//	if ( tri->isTwoUpper )    	glColor3f( 0.9, 0.1, 0.1 );
	//	else                        glColor3f( 0.1, 0.1, 0.9 );

	//	if ( tri->isMark )
	//	{
	//		Vector3f triCen = (tri->v[0]+tri->v[1]+tri->v[2]) / 3.0f;
	//		glBegin(GL_POINTS );
	//		glVertex3f(triCen[0], triCen[1], triCen[2]);
	//		glEnd();

	//		glBegin(GL_LINES);
	//		glVertex3f(tri->tempPts[0][0], tri->tempPts[0][1], tri->tempPts[0][2]);
	//		glVertex3f(tri->tempPts[1][0], tri->tempPts[1][1], tri->tempPts[1][2]);
	//		glEnd();
	//	}
	//}

	glEnable(GL_LIGHTING);
	glPointSize( 1.0 );
	glLineWidth( 1.0 );

}


int LVoxel::getSampledPointsNum()
{
	return sampPts.size();
}

bool LVoxel::isPointIn(Vector3f thisPoint)
{
	Vector3f boxMinPt = center - 0.5f * size;
	Vector3f boxMaxPt = center + 0.5f * size;

	return IsPointInsideBox(boxMinPt, boxMaxPt, thisPoint);
}
