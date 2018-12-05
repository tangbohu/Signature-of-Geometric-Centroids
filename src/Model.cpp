///////////////////////////////////////////////////////////////
//
// Model.cpp
//
//   Reconstructed 3D Model Class
//
// by Song Peng ( songpeng@ustc.edu.cn )
// 
// 11/May/2015
//
///////////////////////////////////////////////////////////////

#include "time.h"
#include "Controls.h"
#include "HelpStructs.h"
#include "HelpDefines.h"
#include "Helper.h"
#include "Model.h"


//**************************************************************************************//
//                                   Initialization
//**************************************************************************************//

Model::Model()
{
	trimesh   = NULL;
	isColored = false;
}

Model::~Model()
{
	ClearScan();
}

void Model::ClearScan()
{
	if ( trimesh != NULL )
	{
		delete trimesh;
		trimesh = NULL;
	}

	verList.clear();

	for (int i=0; i<triList.size(); i++)
	{
		delete triList[i];
	}
	triList.clear();
}

//void Model::InitPose()
//{
//	identity( localMatrix );
//}
//
//void Model::ResetPose()
//{
//	identity( localMatrix );
//}

TriMesh* Model::GetTriMesh()
{
	return trimesh;
}




//**************************************************************************************//
//                                   Load Mesh Model
//**************************************************************************************//

void Model::LoadMesh(const char *fileName)
{
	if( trimesh != NULL )
	{
		delete trimesh;
		trimesh = NULL;
	}
	
	// Read 3D model file (e.g., PLY, OBJ format)
	trimesh = new TriMesh;
	trimesh = trimesh->read(fileName);

	// Calculate various attributes of the trimesh
	trimesh->need_faces();
	trimesh->need_normals();
	trimesh->need_bbox();
	trimesh->need_tstrips();	
	trimesh->need_neighbors();
	trimesh->need_curvatures();
}

void Model::ProcessMesh(const Vector3f diffuseMTL, const float scaleFactor)
{
	// Compute bbox and transform the model
	ComputeBoundingBox();
	SetDiffuseMTL( diffuseMTL );
	NormalizeModel( scaleFactor );

	// Get mesh vertices and triangles
	GetModelVertices( verList );
	GetModelTriangles( triList );
}

void Model::ComputeBoundingBox()
{
	bbox.minPt = trimesh->bbox.min;
	bbox.maxPt = trimesh->bbox.max;
	bbox.cenPt = 0.5f * (bbox.minPt+bbox.maxPt);
}

void Model::SetDiffuseMTL(Vector3f diffuseMTL)
{
	mtlDiffuse[0]  = diffuseMTL[0];    
	mtlDiffuse[1]  = diffuseMTL[1];   
	mtlDiffuse[2]  = diffuseMTL[2];   
	mtlDiffuse[3]  = 1.0;
	
	mtlSpecular[0] = 0.6;   
	mtlSpecular[1] = 0.6;   
	mtlSpecular[2] = 0.6;    
	mtlSpecular[3] = 1.0;

	mtlAmbient[0]  = 0.1;   
	mtlAmbient[1]  = 0.1;   
	mtlAmbient[2]  = 0.1;    
	mtlAmbient[3]  = 1.0;

	mtlEmission[0] = 0.2;   
	mtlEmission[1] = 0.2;   
	mtlEmission[2] = 0.2;    
	mtlEmission[3] = 1.0;
}

void Model::NormalizeModel(const float scale)
{
	// Translate and scale the model to make it within [-a, a] (a is close to 1.0)
	for (int i=0; i<trimesh->vertices.size(); i++) 
	{
		Vector3f origVertex = trimesh->vertices[i];
		Vector3f currVertex = (origVertex - bbox.cenPt) * scale;

		trimesh->vertices[i] = currVertex;
	}

	// Translate and scale the bbox accordingly
	bbox.Transform(-bbox.cenPt, scale);

	// Print the bbox after the normalization
	Vector3f lenVec = bbox.maxPt - bbox.minPt;
	float bboxMaxLen = _MAX(lenVec[0], _MAX(lenVec[1], lenVec[2]));
	printf("\nAfter Normalization:  lenVec [%.3f %.3f %.3f]  maxLen: %.3f \n", lenVec[0], lenVec[1], lenVec[2], bboxMaxLen);
}

bool Model::GetModelVertices(vector<Point> &vertices)
{
	if ( trimesh == NULL )
		return false;

	// Determine if the mesh has vertex color
	if ( trimesh->colors.size() > 0 )
	{
		isColored = true;
		printf("Warning: This is a colored mesh model. \n");
	}
	else
	{
		isColored = false;
	}

	vertices.clear();
	for (int i=0; i<trimesh->vertices.size(); i++) 
	{
		Point pt;
		pt.pos  = trimesh->vertices[i];
		pt.nor  = trimesh->normals[i];
		pt.curv = (trimesh->curv1[i] + trimesh->curv2[i]) / 2.0; // Vertex mean curvature
		if ( isColored )
		{
			pt.color = trimesh->colors[i];
			//printf("i=%d  color: [%.2f %.2f %.2f] \n", i, pt.color[0], pt.color[1], pt.color[2]);
		}

		vertices.push_back( pt );
	}

	return true;
}

bool Model::GetModelTriangles(vector<Triangle*> &triangles)
{
	if ( trimesh == NULL || verList.size() == 0 )
		return false;

	triangles.clear();
	for (int i=0; i<trimesh->faces.size(); i++)
	{
		Triangle *tri = new Triangle();
		//tri->id = i;

		// Get vertex positions of three corners
		for (int j=0; j<3; j++)
		{
			int verIndex  = trimesh->faces[i].v[j];

			tri->v[j]     = verList[verIndex].pos;
			tri->curv[j]  = verList[verIndex].curv;
			tri->color[j] = verList[verIndex].color;
		}

		// Compute the area and the normal
		vec normal  = (tri->v[1] - tri->v[0]) CROSS (tri->v[2] - tri->v[0]);
		tri->area   = len(normal);
		tri->normal = normal / len(normal); 

		triangles.push_back( tri );
	}

	return true;
}




//**************************************************************************************//
//                                   Draw Mesh Model 
//**************************************************************************************//

void Model::DrawModel(bool isAlign)
{
	if ( trimesh == NULL )
		return;

	glEnable(GL_LIGHTING);
	if ( isColored )
	{
		glEnable(GL_COLOR_MATERIAL); 
	}

	glPushAttrib( GL_LIGHTING_BIT );

	// Set Material
	GLfloat shininess[] = { 76.8f };
	glMaterialfv(GL_FRONT_AND_BACK,  GL_AMBIENT,   mtlAmbient);
	glMaterialfv(GL_FRONT_AND_BACK,  GL_DIFFUSE,   mtlDiffuse);
	glMaterialfv(GL_FRONT_AND_BACK,  GL_SPECULAR,  mtlSpecular);
	glMaterialfv(GL_FRONT_AND_BACK,  GL_EMISSION,  mtlEmission);
	glMaterialfv(GL_FRONT_AND_BACK,  GL_SHININESS, shininess);

	// Draw the mesh model
	//glMatrixMode(GL_MODELVIEW);
	//glPushMatrix();
	DrawTStrips( trimesh );
	//glPopMatrix();

	glPopAttrib();

	if ( isColored )
	{
		glDisable(GL_COLOR_MATERIAL); 
	}
}

// Draw triangle strips. They are stored as length followed by values.
void Model::DrawTStrips(const TriMesh *themesh)
{
	// Vertex Positions
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_FLOAT, sizeof(trimesh->vertices[0]), &trimesh->vertices[0][0]);

	// Vertex Colors
	if (!trimesh->colors.empty()) 
	{	
		glEnableClientState(GL_COLOR_ARRAY);
		glColorPointer(3, GL_FLOAT, sizeof(trimesh->colors[0]), &trimesh->colors[0][0]);
	} 
	else       
	{ 
		glDisableClientState(GL_COLOR_ARRAY);  
	}

	// Vertex Normals
	if (!trimesh->normals.empty()) 
	{	
		glEnableClientState(GL_NORMAL_ARRAY);
		glNormalPointer(GL_FLOAT, sizeof(trimesh->normals[0]), &trimesh->normals[0][0]);
	} 
	else       
	{ 
		glDisableClientState(GL_NORMAL_ARRAY);  
	}

	// Draw the triangles
	const int *t = &themesh->tstrips[0];
	const int *end = t + themesh->tstrips.size();
	while ( likely(t < end) )
	{
		int striplen = *t++;
		glDrawElements(GL_TRIANGLE_STRIP, striplen, GL_UNSIGNED_INT, t);
		t += striplen;
	}
}

void Model::DrawModelWire(Vector3f color, float width)
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
	// Draw mesh wire frame
	for (int i=0; i<triList.size(); i++)
	{
		Triangle *tri = triList[i];

		glBegin(GL_LINE_LOOP);
		glVertex3f(tri->v[0][0], tri->v[0][1], tri->v[0][2]);
		glVertex3f(tri->v[1][0], tri->v[1][1], tri->v[1][2]);
		glVertex3f(tri->v[2][0], tri->v[2][1], tri->v[2][2]);
		glEnd();
	}
	glPopMatrix();

	glMatrixMode(GL_MODELVIEW);
	glLineWidth( 1.0 );
	glEnable(GL_LIGHTING);
}

void Model::DrawModelBBox(Vector3f color, float width)
{
	if ( trimesh == NULL )
		return;

	DrawWireCuboid(bbox.minPt, bbox.maxPt, width, color);
}