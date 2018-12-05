///////////////////////////////////////////////////////////////
//
// Scan.cpp
//
//   General Triangular Scan Mesh Class
//
// by Song Peng ( songpeng@ustc.edu.cn )
// 
// 11/May/2015
//
///////////////////////////////////////////////////////////////

#include "lineqn.h"
#include "math3D.h"
#include "time.h"
#include "RayMeshTest.h"
#include "Controls.h"
#include "HelpDefines.h"
#include "Helper.h"
#include "Scan.h"


//**************************************************************************************//
//                                   Initialization
//**************************************************************************************//

Scan::Scan()
{
	trimesh   = NULL;
	isColored = false;
	scanID    = 0;

	::identity( localMatrix );
	::identity(globalMatrix);
	::identity(alignLocalMat);
	::identity(alignGlobalMat);

	//isAligned = false;
}

Scan::~Scan()
{
	ClearScan();
}

void Scan::ClearScan()
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

	//for (int i=0; i<featDescriptors.size(); i++)
	//{
	//	delete featDescriptors[i];
	//}
	//featDescriptors.clear();

	for (int i=0; i<sampDescriptors.size(); i++)
	{
		delete sampDescriptors[i];
	}
	sampDescriptors.clear();
}

void Scan::InitPose()
{
	::identity(localMatrix);

	localMatrix[12]   = INITIAL_SCAN_SHIFT_VECTOR[0];
	localMatrix[13]   = INITIAL_SCAN_SHIFT_VECTOR[1];
	localMatrix[14]   = INITIAL_SCAN_SHIFT_VECTOR[2];

	::identity(alignLocalMat);

	alignLocalMat[12] = ALIGNED_SCAN_SHIFT_VECTOR[0];
	alignLocalMat[13] = ALIGNED_SCAN_SHIFT_VECTOR[1];
	alignLocalMat[14] = ALIGNED_SCAN_SHIFT_VECTOR[2];
}

void Scan::ResetPose()
{
	::identity(localMatrix);

	localMatrix[12]   = INITIAL_SCAN_SHIFT_VECTOR[0];
	localMatrix[13]   = INITIAL_SCAN_SHIFT_VECTOR[1];
	localMatrix[14]   = INITIAL_SCAN_SHIFT_VECTOR[2];
}

TriMesh* Scan::GetTriMesh()
{
	return trimesh;
}




//**************************************************************************************//
//                                   Load Mesh Model
//**************************************************************************************//

void Scan::LoadMesh(const char *fileName)
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

void Scan::ProcessMesh(const Vector3f diffuseMTL, const float scaleFactor)
{
	// Compute bbox and transform the model
	ComputeBoundingBox();
	SetDiffuseMTL( diffuseMTL );
	NormalizeModel( scaleFactor );
	CheckCurvatures();

	// Get mesh vertices and triangles
	GetModelVertices( verList );
	GetModelTriangles( triList );
}

void Scan::ComputeBoundingBox()
{
	bbox.minPt = trimesh->bbox.min;
	bbox.maxPt = trimesh->bbox.max;
	bbox.cenPt = 0.5f * (bbox.minPt+bbox.maxPt);
}

void Scan::SetDiffuseMTL(Vector3f diffuseMTL)
{
	mtlDiffuse[0]  = diffuseMTL[0];    
	mtlDiffuse[1]  = diffuseMTL[1];   
	mtlDiffuse[2]  = diffuseMTL[2];   
	mtlDiffuse[3]  = 1.0;
	
	mtlSpecular[0] = 0.5;   
	mtlSpecular[1] = 0.5;   
	mtlSpecular[2] = 0.5;    
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

float Scan::GetScaleFactor()
{
	ComputeBoundingBox();

	vec lenVec = bbox.maxPt - bbox.minPt;
	float bboxMaxLen = _MAX(lenVec[0], _MAX(lenVec[1], lenVec[2]));
	printf("Before Normalization: lenVec [%.3f %.3f %.3f]  maxLen: %.3f \n", lenVec[0], lenVec[1], lenVec[2], bboxMaxLen);

	// Translate and scale the model to make it within [-a, a] (a is close to 1.0)
	float scaleFactor = 2.0f / bboxMaxLen;
	printf("scaleFactor: %.6f \n", scaleFactor);

	return scaleFactor;
}

void Scan::NormalizeModel(const float scale)
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

void Scan::CheckCurvatures()
{
	//////////////////////////////////////////////////////////////////////////
	// 1. Compute min and max curvature values for two principle directions 

	float minCurv1 = MAX_FLOAT; 
	float maxCurv1 = MIN_FLOAT;  
	float minCurv2 = MAX_FLOAT;
	float maxCurv2 = MIN_FLOAT;

	float minCurv  = MAX_FLOAT;
	float maxCurv  = MIN_FLOAT;

	for (int i=0; i<trimesh->vertices.size(); i++)
	{
		if ( minCurv1 > trimesh->curv1[i] )    minCurv1 = trimesh->curv1[i];
		if ( maxCurv1 < trimesh->curv1[i] )    maxCurv1 = trimesh->curv1[i];

		if ( minCurv2 > trimesh->curv2[i] )    minCurv2 = trimesh->curv2[i];
		if ( maxCurv2 < trimesh->curv2[i] )    maxCurv2 = trimesh->curv2[i];

		float verCurv = (trimesh->curv1[i]+trimesh->curv2[i]) / 2.0;

		if ( minCurv  > verCurv )     minCurv  = verCurv;
		if ( maxCurv  < verCurv )     maxCurv  = verCurv;
	}

	printf("Curv1 Range: [%.2f %.2f]  \n", minCurv1, maxCurv1);
	printf("Curv2 Range: [%.2f %.2f]  \n", minCurv2, maxCurv2);
	printf("Curv  Range: [%.2f %.2f]  \n", minCurv, maxCurv);


	//////////////////////////////////////////////////////////////////////////
	// 2. Remove the vertex curvature value that is too small or too large 

	//const float curvUpperThres =  100;
	//const float curvLowerThres = -100;

	//for (int i=0; i<trimesh->vertices.size(); i++)
	//{
	//	if ( trimesh->curv1[i] < curvLowerThres )    trimesh->curv1[i] = curvLowerThres;
	//	if ( trimesh->curv1[i] > curvUpperThres )    trimesh->curv1[i] = curvUpperThres;

	//	if ( trimesh->curv2[i] < curvLowerThres )    trimesh->curv2[i] = curvLowerThres;
	//	if ( trimesh->curv2[i] > curvUpperThres )    trimesh->curv2[i] = curvUpperThres;
	//}
}

bool Scan::GetModelVertices(vector<Point> &vertices)
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

bool Scan::GetModelTriangles(vector<Triangle*> &triangles)
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
//                              Describe Scan with Descriptors
//**************************************************************************************//

void Scan::SampleScan()
{
	sampler.SampleMesh_Uniform( verList );
	sampler.SampleMesh_Feature( verList );
}

void Scan::DescribeScan()
{
	ComputeSampleDescriptors();
	ComputeFeatureDescriptors();
}

void Scan::ComputeFeatureDescriptors()
{

	// Get feature descriptors from sample descriptors
	vector<int> featPtIndices = sampler.GetFeaturePointIndices();
	for (int i=0; i<featPtIndices.size(); i++)
	{
		int index = featPtIndices[i];
#ifdef USE_DESCRIPTOR_SGC
		featSGCDescriptors.push_back(sampSGCDescriptors[index]);
#else	
		featDescriptors.push_back(sampDescriptors[index]);
#endif

	}
}

//void Scan::ComputeFeatureDescriptors()
//{
//	for (int i=0; i<featDescriptors.size(); i++)
//	{
//		delete featDescriptors[i];
//	}
//	featDescriptors.clear();
//
//	vector<Point> scanFeatPts = sampler.GetFeaturePoints();
//	for (int i=0; i<scanFeatPts.size(); i++)
//	{
//		Descriptor* descriptor = BuildLocalVoxelizer(scanFeatPts[i], LOCAL_SHAPE_RADIUS, false);	
//		featDescriptors.push_back( descriptor );
//	}
//}

void Scan::ComputeSampleDescriptors()
{
	for (int i=0; i<sampDescriptors.size(); i++)
	{
		delete sampDescriptors[i];
	}
	sampDescriptors.clear();

	for (int i = 0; i<sampSGCDescriptors.size(); i++)
	{
		delete sampSGCDescriptors[i];
	}
	sampSGCDescriptors.clear();


#ifdef USE_DESCRIPTOR_SGC
	cout << "compute SGC" << endl;
#endif

	vector<Point> scanSampPts = sampler.GetSamplePoints();
	for (int i=0; i<scanSampPts.size(); i++)
	{
#ifdef USE_DESCRIPTOR_OURS
	Descriptor* descriptor = BuildLocalVoxelizer(scanSampPts[i], LOCAL_SHAPE_RADIUS, false);
	sampDescriptors.push_back( descriptor );

#elif USE_DESCRIPTOR_SGC
		SGCentroid* descriptor = BuildSGC(scanSampPts[i], LOCAL_SHAPE_RADIUS, false);
		sampSGCDescriptors.push_back(descriptor);
#endif

	}
}




Descriptor* Scan::BuildLocalVoxelizer(Point keyPoint, float radius, bool isDebug)
{
	// Get local shape triangles
	vector<Triangle*> LRFTris;     
	vector<Triangle*> shapeTris; 
	vector<Point> LRFPts;
	vector<Point> shapePts;
	GetLocalShape_Tri(keyPoint, radius, LRFTris, shapeTris);
	GetLocalShape_Ver(keyPoint, radius, LRFPts, shapePts);

	// Compute descriptor for the local shape
	Descriptor *descriptor; 
	if ( isDebug )
	{
		// Define voxelizer as a global variable for debug
		voxelizer.InitLVoxelizer(keyPoint, radius, LRFTris, shapeTris, LRFPts, shapePts);
		descriptor = voxelizer.ComputeDescriptor( FEATURE_TYPE_ID, true );
	}
	else
	{
		// Define voxelizer as a local variable for generating results
		LVoxelizer tempVoxelizer;
		tempVoxelizer.InitLVoxelizer(keyPoint, radius, LRFTris, shapeTris, LRFPts, shapePts);
		descriptor = tempVoxelizer.ComputeDescriptor( FEATURE_TYPE_ID, false );
	}

	return descriptor;
}

SGCentroid* Scan::BuildSGC(Point keyPoint, float radius, bool isDebug)
{
	// Get local shape triangles
	vector<Triangle*> LRFTris;
	vector<Triangle*> shapeTris;
	vector<Point> LRFPts;
	vector<Point> shapePts;
	GetLocalShape_Tri(keyPoint, radius, LRFTris, shapeTris);
	GetLocalShape_Ver(keyPoint, radius, LRFPts, shapePts);

	// Compute descriptor for the local shape
	SGCentroid *descriptor=NULL;

	if (!isDebug)
	{
		// Define voxelizer as a global variable for debug
		voxelizer.InitLVoxelizer(keyPoint, radius, LRFTris, shapeTris, LRFPts, shapePts);
		descriptor = voxelizer.ComputeSGC( );
	}
	else
	{
		// Define voxelizer as a local variable for generating results
		LVoxelizer tempVoxelizer;
		tempVoxelizer.InitLVoxelizer(keyPoint, radius, LRFTris, shapeTris, LRFPts, shapePts);
		descriptor = tempVoxelizer.ComputeSGC();
	}

	return descriptor;
}


void Scan::UpdateLVoxelizer(Descriptor* descriptor)
{
	LRF keyPtLRF = descriptor->keyPtLRF;
	float radius = descriptor->keyPtLRF.radius;
	vector<Triangle*> dumyTris;
	vector<Point> dumyPts;

	voxelizer.InitLVoxelizer(keyPtLRF.keyPoint, keyPtLRF.radius, dumyTris, dumyTris, dumyPts, dumyPts);
	voxelizer.localVolume.InitLVolume(LOCAL_VOLUME_DIMEN, vec(2 * radius, 2 * radius, 2 * radius), dumyTris);

	voxelizer.keyPtLRF = keyPtLRF;
}

Descriptor* Scan::BuildLocalSphere(Point keyPoint, float radius, bool isDebug)
{
	// Get local shape triangles
	vector<Triangle*> LRFTris;     
	vector<Triangle*> shapeTris; 
	vector<Point> shapePts;
	vector<Point> LRFPts;

	GetLocalShape_Tri(keyPoint, radius, LRFTris, shapeTris);
	GetLocalShape_Ver(keyPoint, radius, LRFPts, shapePts);

	// Compute descriptor for the local shape
	Descriptor *descriptor; 
	if ( isDebug )
	{
		// Define voxelizer as a global variable for debug
		voxelizer.InitLVoxelizer(keyPoint, radius, LRFTris, shapeTris, LRFPts, shapePts);
		descriptor = voxelizer.ComputeDescriptor_Spherical( FEATURE_TYPE_ID, true);
	}
	else
	{
		// Define voxelizer as a local variable for generating results
		LVoxelizer tempVoxelizer;
		tempVoxelizer.InitLVoxelizer(keyPoint, radius, LRFTris, shapeTris, LRFPts, shapePts);
		descriptor = tempVoxelizer.ComputeDescriptor_Spherical( FEATURE_TYPE_ID, false);
	}

	return descriptor;
}




//**************************************************************************************//
//                               Get Mesh Local Shape
//**************************************************************************************//

void Scan::GetLocalShape_Tri(Point keyPoint, float radius, vector<Triangle*> &LRFTris, vector<Triangle*> &shapeTris)
{
	shapeTris.clear();
	LRFTris.clear();

	for (int i=0; i<triList.size(); i++)
	{
		// Save all the triangles in local sphere (with radius sqrt(2)*R)
		if ( IsTriangleInsideSphere(*triList[i], keyPoint.pos, SQUARE_ROOT_TWO*radius) )
		{
			shapeTris.push_back( triList[i] );

			// Save all the triangles in local sphere (with radius R)
			if ( IsTriangleInsideSphere(*triList[i], keyPoint.pos, radius) )
			{
				LRFTris.push_back( triList[i] );
			}
		}
	}

	//printf("localTriNum: %d \n", LRFTris.size());
}

void Scan::GetLocalShape_Ver(Point keyPoint, float radius, vector<Point> &LRFPts, vector<Point> &shapePts)
{
	shapePts.clear();

	for (int i=0; i<verList.size(); i++)
	{
		if (IsVertexInsideSphere(verList[i], keyPoint.pos, SQUARE_ROOT_TWO*radius))
		{
			shapePts.push_back(verList[i]);

			if (IsVertexInsideSphere(verList[i], keyPoint.pos, radius))
			{
				LRFPts.push_back(verList[i]);
			}
		}
	}
}

bool Scan::IsTriangleInsideSphere(Triangle tri, Vector3f center, float radius)
{
	bool triInSphere = false;

	// The triangle is considered as inside the sphere if there exists at least one
	// vertex of it inside the sphere
	for (int i=0; i<3; i++)
	{
		float dist = trimesh::dist(tri.v[i], center);

		if ( dist < radius)
		{
			triInSphere = true;
			break;
		}
	}

	return triInSphere;
}

bool Scan::IsVertexInsideSphere(Point ver, Vector3f center, float radius)
{
	float dist = len(ver.pos -center);

	if ( dist < radius )
		return true;
	else
		return false;
}




//**************************************************************************************//
//                               Pick Scan Mesh Point
//**************************************************************************************//

Point Scan::PickMeshSurfPoint(int winX, int winY)
{
	////////////////////////////////////////////////////////////////
	// 1. Construct the ray from the origin to the picking point

	double inveGlobalMat[16], tranInveGlobalMat[16];
	memcpy(inveGlobalMat, globalMatrix, sizeof(double)*16);
	if(invertMat(inveGlobalMat, 4)==1)  printf("Inverse Matrix Error \n");
	memcpy(tranInveGlobalMat, inveGlobalMat, sizeof(double)*16);
	transposeMat(tranInveGlobalMat, 4);

	double tempPos[3] = {0, 0, 0};
	origPoint[0] = inveGlobalMat[12] + dot3D(tempPos, tranInveGlobalMat  );
	origPoint[1] = inveGlobalMat[13] + dot3D(tempPos, tranInveGlobalMat+4);
	origPoint[2] = inveGlobalMat[14] + dot3D(tempPos, tranInveGlobalMat+8);

	// TODO: need to check why need to add 20 pixels
	pickPoint = GetOGLPos(winX, winY+20);
	Vector3f rayDir = (pickPoint-origPoint) /len(pickPoint-origPoint);


	////////////////////////////////////////////////////////////////
	// 2. Compute picking point position by ray-mesh intersection

	HitPoint hitPoint;
	bool isIntersect = LineMeshIntersect(origPoint, rayDir, triList, hitPoint);

	if ( isIntersect )
	{
		intersetPt.pos = hitPoint.pos;
		intersetPt.nor = hitPoint.nor;
	}
	else
	{
		printf("Warning: cannot find an intersecting point between the line and the mesh. \n");
	}

	return intersetPt;
}

Vector3f Scan::GetOGLPos(int winX, int winY)
{
	GLint viewport[4];
	GLdouble projection[16];
	GLfloat winZ;
	GLdouble posX, posY, posZ;

	glGetDoublev( GL_PROJECTION_MATRIX, projection );
	glGetIntegerv( GL_VIEWPORT, viewport );
	glReadPixels( winX, winY, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ );
	gluUnProject( (float)winX, (float)winY, winZ, globalMatrix, projection, viewport, &posX, &posY, &posZ);

	return Vector3f(posX, posY, posZ);
}




//**************************************************************************************//
//                                   Draw Mesh Model 
//**************************************************************************************//

void Scan::DrawScan(bool isAlign)
{
	if ( trimesh == NULL )
		return;

	//isColored = false;

	glEnable(GL_LIGHTING);
	if ( isColored )
	{
		glEnable(GL_COLOR_MATERIAL); 
	}

	glPushAttrib( GL_LIGHTING_BIT );

	// Set Material
	//GLfloat shininess[] = { 76.8f };
	GLfloat shininess[] = { 200.0f };  // For Chicken model
	glMaterialfv(GL_FRONT_AND_BACK,  GL_AMBIENT,   mtlAmbient);
	glMaterialfv(GL_FRONT_AND_BACK,  GL_DIFFUSE,   mtlDiffuse);
	glMaterialfv(GL_FRONT_AND_BACK,  GL_SPECULAR,  mtlSpecular);
	glMaterialfv( GL_FRONT_AND_BACK, GL_EMISSION,  mtlEmission);
	glMaterialfv( GL_FRONT_AND_BACK, GL_SHININESS, shininess);

	// Draw the mesh model
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	if ( isAlign )  glLoadMatrixd( alignGlobalMat );
	else		    glLoadMatrixd( globalMatrix );
//	DrawTStrips( trimesh );
	DrawPoints(trimesh);
	glPopMatrix();

	glPopAttrib();

	if ( isColored )
	{
		glDisable(GL_COLOR_MATERIAL); 
	}
}

// Draw triangle strips. They are stored as length followed by values.
void Scan::DrawTStrips(const TriMesh *themesh)
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

void Scan::DrawPoints(const  TriMesh *ply)
{

	for (int i = 0; i<ply->vertices.size(); i++)
	{
	
		glBegin(GL_POINTS);
		if (!ply->colors.empty()) 
		{
			glColor3f(ply->colors[i][0], ply->colors[i][1], ply->colors[i][2]);
		}
		glNormal3f(ply->normals[i][0], ply->normals[i][1], ply->normals[i][2]);
		glVertex3f(ply->vertices[i][0], ply->vertices[i][1], ply->vertices[i][2]);
		glEnd();

	}
}

void Scan::DrawScanWire(Vector3f color, float width)
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
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadMatrixd( globalMatrix );
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

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glMatrixMode(GL_MODELVIEW);
	glLineWidth( 1.0 );
	glEnable(GL_LIGHTING);
}

void Scan::DrawScanBBox()
{
	if ( trimesh == NULL )
		return;

	Vector3f color;
	color[0] = mtlDiffuse[0];
	color[1] = mtlDiffuse[1];
	color[2] = mtlDiffuse[2];

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadMatrixd( globalMatrix );
	DrawWireCuboid(bbox.minPt, bbox.maxPt, 3.0, color);
	glPopMatrix();
}

void Scan::DrawFeaturePoints(bool isAlign)
{
	if ( trimesh == NULL )
		return;

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	//glLoadMatrixd( globalMatrix );
	if ( isAlign )  glLoadMatrixd( alignGlobalMat );
	else		    glLoadMatrixd( globalMatrix );
	sampler.DrawFeaturePoints();
	glPopMatrix();
}

void Scan::DrawSamplePoints(bool isAlign)
{
	if ( trimesh == NULL )
		return;

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	//glLoadMatrixd( globalMatrix );
	if ( isAlign )  glLoadMatrixd( alignGlobalMat );
	else		    glLoadMatrixd( globalMatrix );
	sampler.DrawSamplePoints();
	glPopMatrix();
}




//**************************************************************************************//
//                           Draw Local Shape, LRF, Volume etc.
//**************************************************************************************//

void Scan::DrawRay()
{
	glDisable(GL_LIGHTING);

	glLineWidth(3.0);
	glPointSize(12.0);
	glColor3f(0.8,0.2,0.2);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadMatrixd( globalMatrix );

	glBegin(GL_LINES);
	glVertex3f(origPoint[0], origPoint[1], origPoint[2]);
	glVertex3f(pickPoint[0], pickPoint[1], pickPoint[2]);
	glEnd(); 

	//glBegin(GL_POINTS);
	//glVertex3f(origPoint[0], origPoint[1], origPoint[2]);
	//glVertex3f(pickPoint[0], pickPoint[1], pickPoint[2]);
	//glEnd();

	glColor3f(0.2,0.2,0.8);
	glBegin(GL_POINTS);
	glVertex3f(intersetPt.pos[0],  intersetPt.pos[1],  intersetPt.pos[2]);
	glEnd();

	glLineWidth(1.0);
	glPointSize(1.0);

	glPopMatrix();
}

void Scan::DrawLRFShape()
{
	if ( trimesh == NULL )
		return;

	glPushAttrib( GL_LIGHTING_BIT );

	GLfloat shininess[] = { 76.8f };
	glMaterialfv(GL_FRONT_AND_BACK,  GL_AMBIENT,   mtlAmbient);
	glMaterialfv(GL_FRONT_AND_BACK,  GL_DIFFUSE,   mtlDiffuse);
	glMaterialfv(GL_FRONT_AND_BACK,  GL_SPECULAR,  mtlSpecular);
	glMaterialfv( GL_FRONT_AND_BACK, GL_EMISSION,  mtlEmission);
	glMaterialfv( GL_FRONT_AND_BACK, GL_SHININESS, shininess);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
#ifdef SHOW_LRF_ALIGNED_SCAN
	glLoadMatrixd( alignGlobalMat );
#else
	glLoadMatrixd( globalMatrix );
#endif
	//glMultMatrixd( voxelizer.keyPtLRF.LRFMatrix );
	voxelizer.keyPtLRF.DrawLRFShape();
	glPopMatrix();

	glPopAttrib();
}

void Scan::DrawKeyPoint()
{
	if ( trimesh == NULL )
		return;

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
#ifdef SHOW_LRF_ALIGNED_SCAN
	glLoadMatrixd( alignGlobalMat );
#else
	glLoadMatrixd( globalMatrix );
#endif
	
	voxelizer.keyPtLRF.DrawKeyPoint( Vector3f(0.9,0.2,0.8) );

	glPopMatrix();
}


void Scan::DrawLRFSphere()
{
	if ( trimesh == NULL )
		return;

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
#ifdef SHOW_LRF_ALIGNED_SCAN
	glLoadMatrixd( alignGlobalMat );
#else
	glLoadMatrixd( globalMatrix );
#endif
	voxelizer.keyPtLRF.DrawLRFSphere( Vector3f(0.2,0.2,0.2) );

	glPopMatrix();
}


void Scan::DrawLRF()
{
	if ( trimesh == NULL )
		return;

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
#ifdef SHOW_LRF_ALIGNED_SCAN
	glLoadMatrixd( alignGlobalMat );
#else
	glLoadMatrixd( globalMatrix );
#endif
	voxelizer.keyPtLRF.DrawLRF();

	glPopMatrix();
}

void Scan::DrawLocalShape()
{
	if ( trimesh == NULL )
		return;

	glPushAttrib( GL_LIGHTING_BIT );

	GLfloat shininess[] = { 76.8f };
	glMaterialfv(GL_FRONT_AND_BACK,  GL_AMBIENT,   mtlAmbient);
	glMaterialfv(GL_FRONT_AND_BACK,  GL_DIFFUSE,   mtlDiffuse);
	glMaterialfv(GL_FRONT_AND_BACK,  GL_SPECULAR,  mtlSpecular);
	glMaterialfv( GL_FRONT_AND_BACK, GL_EMISSION,  mtlEmission);
	glMaterialfv( GL_FRONT_AND_BACK, GL_SHININESS, shininess);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
#ifdef SHOW_LRF_ALIGNED_SCAN
	glLoadMatrixd( alignGlobalMat );
#else
	glLoadMatrixd( globalMatrix );
#endif
	glMultMatrixd( voxelizer.keyPtLRF.LRFMatrix );
	voxelizer.localVolume.DrawLocalShape();
	glPopMatrix();

	glPopAttrib();
}

void Scan::DrawLocalShapeWire()
{
	if ( trimesh == NULL )
		return;

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
#ifdef SHOW_LRF_ALIGNED_SCAN
	glLoadMatrixd( alignGlobalMat );
#else
	glLoadMatrixd( globalMatrix );
#endif
	glMultMatrixd( voxelizer.keyPtLRF.LRFMatrix );
	voxelizer.localVolume.DrawLocalShapeWire(vec(0.1,0.1,0.1), 1.0);
	glPopMatrix();
}

void Scan::DrawLocalBBox()
{
	if ( trimesh == NULL )
		return;

	glDisable(GL_LIGHTING);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
#ifdef SHOW_LRF_ALIGNED_SCAN
	glLoadMatrixd( alignGlobalMat );
#else
	glLoadMatrixd( globalMatrix );
#endif
	glMultMatrixd( voxelizer.keyPtLRF.LRFMatrix );
	voxelizer.localVolume.DrawLocalBBox(4.0, vec(0.8, 0.6, 0.1));
	glPopMatrix();

	glEnable(GL_LIGHTING);
	glLineWidth(1.0f);
}

void Scan::DrawLocalGrid()
{
	if ( trimesh == NULL )
		return;

	glColor3f(0.8,0.2,0.2);
	glDisable(GL_LIGHTING);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
#ifdef SHOW_LRF_ALIGNED_SCAN
	glLoadMatrixd( alignGlobalMat );
#else
	glLoadMatrixd( globalMatrix );
#endif
	glMultMatrixd( voxelizer.keyPtLRF.LRFMatrix );
	//DrawWireCuboid(-voxelizer.radius*vec(1,1,1), voxelizer.radius*vec(1,1,1), 2.0, vec(0.9, 0.9, 0.1) );
	voxelizer.localVolume.DrawVoxelGrid(2.0, vec(0.1, 0.1, 0.1));
	//localVolume.DrawCrossPoints();
	//localVolume.DrawPointGrid(3.0, Vector3f(0.9,0.2,0.2));
	voxelizer.localVolume.DrawVolumeSamplePts(3.0, Vector3f(0.8,0.5,0.1));
	glPopMatrix();

	glEnable(GL_LIGHTING);
	glLineWidth(1.0f);
}

void Scan::DrawLocalSphere()
{
	if ( trimesh == NULL )
		return;

	glDisable(GL_LIGHTING);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
#ifdef SHOW_LRF_ALIGNED_SCAN
	glLoadMatrixd( alignGlobalMat );
#else
	glLoadMatrixd( globalMatrix );
#endif
	glMultMatrixd( voxelizer.keyPtLRF.LRFMatrix );
	voxelizer.localSphere.DrawLocalSphere(3.0, Vector3f(0.2,0.2,0.2));
	glPopMatrix();

	glEnable(GL_LIGHTING);
	glLineWidth(1.0f);
}