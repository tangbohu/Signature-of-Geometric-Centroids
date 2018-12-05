///////////////////////////////////////////////////////////////
//
// LVolume.cpp
//
//   Build Local Volume for Describing Local Shape
//
// by Song Peng ( songpeng@ustc.edu.cn )
// 
// 11/May/2015
//
///////////////////////////////////////////////////////////////

#include "Controls.h"
#include "Helper.h"
#include "math3D.h"
#include "BoxTriTest.h"
#include "RayMeshTest.h"
#include "LVoxel.h"
#include "LPointGrid.h"
#include "LVolume.h"


//**************************************************************************************//
//                                   Initialization
//**************************************************************************************//

LVolume::LVolume()
{
	pointGrid = NULL;
}

LVolume::~LVolume()
{
	ClearLVolume();
}

void LVolume::ClearLVolume()
{
	for (int i=0; i<sphereTris.size(); i++)
	{
		delete sphereTris[i];
	}
	sphereTris.clear();

	if( pointGrid != NULL )	   
	{
		delete pointGrid;
		pointGrid = NULL;
	}

	voxelGrid.clear();
}

void LVolume::InitLVolumeWithTransformPoints(Vector3i _volumeSize, Vector3f _volumeDimen, vector<Point> _volumePoints)
{
	ClearLVolume();

	volumeSize = _volumeSize;
	volumeDimen = _volumeDimen;
	volumePoints = _volumePoints;

	// Volume bounding box
	volumeMinPt = -0.5f * volumeDimen;
	volumeMaxPt = 0.5f * volumeDimen;
	volumeCenPt = Vector3f(0, 0, 0);

	// Voxel size
	voxelSize[0] = volumeDimen[0] / (float)volumeSize[0];
	voxelSize[1] = volumeDimen[1] / (float)volumeSize[1];
	voxelSize[2] = volumeDimen[2] / (float)volumeSize[2];
}


void LVolume::InitLVolume(Vector3i _volumeSize, Vector3f _volumeDimen, vector<Triangle*> _sphereTris)
{
	ClearLVolume();

	volumeSize   = _volumeSize;
	volumeDimen  = _volumeDimen;
	sphereTris   = _sphereTris;

	// Volume bounding box
	volumeMinPt  = -0.5f * volumeDimen;
	volumeMaxPt  =  0.5f * volumeDimen;
	volumeCenPt  =  Vector3f(0,0,0);

	// Voxel size
	voxelSize[0] = volumeDimen[0]/(float)volumeSize[0];
	voxelSize[1] = volumeDimen[1]/(float)volumeSize[1];
	voxelSize[2] = volumeDimen[2]/(float)volumeSize[2];
}




//**************************************************************************************//
//                                Compute Voxel Grid
//**************************************************************************************//

void LVolume::ComputeVoxelGrid()
{
	voxelGrid.clear();

	for( int k=0; k< volumeSize[2]; k++ )
	for( int j=0; j< volumeSize[1]; j++ )
	for( int i=0; i< volumeSize[0]; i++ )
	{
		LVoxel voxel;
		voxel.state = VOXEL_OUT_MESH;

		voxel.pos[0] = i; 
		voxel.pos[1] = j; 
		voxel.pos[2] = k; 

		voxel.size = voxelSize; 
		voxel.center[0] = volumeMinPt[0] + (i+0.5)*voxelSize[0]; 
		voxel.center[1] = volumeMinPt[1] + (j+0.5)*voxelSize[1]; 
		voxel.center[2] = volumeMinPt[2] + (k+0.5)*voxelSize[2]; 

		voxelGrid.push_back( voxel );
	}
}

void LVolume::ComputeVoxelState()
{
	//if ( sphereTris.size() ==0 )
	//	return;

	volumeTris.clear();

	for (int i=0; i<sphereTris.size(); i++)
	{
		bool isTriInVolume = false;
		Triangle *tri = sphereTris[i];

		/////////////////////////////////////////////////////////////////////
		// 1. A quick check to prune non-intersecting small triangles

		Vector3i triVerGrids[3];
		triVerGrids[0][0] = MapToGrid( tri->v[0][0], volumeMinPt[0], volumeMaxPt[0], volumeSize[0] );
		triVerGrids[0][1] = MapToGrid( tri->v[0][1], volumeMinPt[1], volumeMaxPt[1], volumeSize[1] );
		triVerGrids[0][2] = MapToGrid( tri->v[0][2], volumeMinPt[2], volumeMaxPt[2], volumeSize[2] );

		triVerGrids[1][0] = MapToGrid( tri->v[1][0], volumeMinPt[0], volumeMaxPt[0], volumeSize[0] );
		triVerGrids[1][1] = MapToGrid( tri->v[1][1], volumeMinPt[1], volumeMaxPt[1], volumeSize[1] );
		triVerGrids[1][2] = MapToGrid( tri->v[1][2], volumeMinPt[2], volumeMaxPt[2], volumeSize[2] );

		triVerGrids[2][0] = MapToGrid( tri->v[2][0], volumeMinPt[0], volumeMaxPt[0], volumeSize[0] );
		triVerGrids[2][1] = MapToGrid( tri->v[2][1], volumeMinPt[1], volumeMaxPt[1], volumeSize[1] );
		triVerGrids[2][2] = MapToGrid( tri->v[2][2], volumeMinPt[2], volumeMaxPt[2], volumeSize[2] );


		int minX = _MIN( triVerGrids[0][0], _MIN(triVerGrids[1][0], triVerGrids[2][0]) );
		int maxX = _MAX( triVerGrids[0][0], _MAX(triVerGrids[1][0], triVerGrids[2][0]) );
		int minY = _MIN( triVerGrids[0][1], _MIN(triVerGrids[1][1], triVerGrids[2][1]) );
		int maxY = _MAX( triVerGrids[0][1], _MAX(triVerGrids[1][1], triVerGrids[2][1]) );
		int minZ = _MIN( triVerGrids[0][2], _MIN(triVerGrids[1][2], triVerGrids[2][2]) );
		int maxZ = _MAX( triVerGrids[0][2], _MAX(triVerGrids[1][2], triVerGrids[2][2]) );


		/////////////////////////////////////////////////////////////////////
		// 2. Perform intersection test between the triangle and related voxels

		for (int l=minX; l<=maxX; l++)
		for (int m=minY; m<=maxY; m++)
		for (int n=minZ; n<=maxZ; n++)
		{
			Vector3i voxelPos = Vector3i(l, m , n);
			int voxelIndex = GetVoxelIndex( voxelPos );

			// The triangle could be out of the volume (since the sphere is larger than the volume)
			if ( voxelIndex < 0 || voxelIndex > voxelGrid.size() )
			{
				continue;
			}

			// Perform intersection test and update voxel state
			int triState = BoxTriangleIntersect(voxelGrid[voxelIndex].center, voxelSize, tri->v);
			if ( triState == TRIANGLE_IN_BOX || triState == TRIANGLE_CROSS_BOX )
			{
				isTriInVolume = true;

				if ( triState == TRIANGLE_IN_BOX    )    tri->isCross = false;
				if ( triState == TRIANGLE_CROSS_BOX )    tri->isCross = true;

				voxelGrid[voxelIndex].state = VOXEL_CROSS_MESH;
				voxelGrid[voxelIndex].AddTriangle( tri );
			}
		}


		/////////////////////////////////////////////////////////////////////
		// 3. Update triangles inside the volume (for visualization)

		if ( isTriInVolume )
		{
			volumeTris.push_back( tri ); // Save the triangles in the local volume
		}
	}
}

int LVolume::MapToGrid(float value , float min , float max , int N)
{
	if (value < min)		   return -1;
	else if (value > max)      return -1;
	else if (value == max)     return N-1;
	else					   return int((value - min) / (max-min) * N);
}

int LVolume::GetVoxelIndex(Vector3i voxelPos)
{
	if( voxelPos[0] < 0 || voxelPos[0] > volumeSize[0]-1 ||
		voxelPos[1] < 0 || voxelPos[1] > volumeSize[1]-1 ||
		voxelPos[2] < 0 || voxelPos[2] > volumeSize[2]-1 )
	{
		return -1;
	}
	else
	{
		int index = voxelPos[2]*volumeSize[1]*volumeSize[0] + voxelPos[1]*volumeSize[0] + voxelPos[0];
		return index;
	}
}




//**************************************************************************************//
//                            Voxel Local Mesh Clipping 
//**************************************************************************************//

void LVolume::ClipMeshInVoxels()
{
	for (int i=0; i<voxelGrid.size(); i++)
	{
		voxelGrid[i].ClipLocalMesh();
	}
}




//**************************************************************************************//
//                                Compute Point Grid  
//**************************************************************************************//



void LVolume::ComputePointGrid()
{
	// Create Point grid (note: need to test this change)
	if ( pointGrid != NULL )
	{
		delete pointGrid;
	}
	pointGrid = new LPointGrid();
	pointGrid->InitLPointGrid(volumeMinPt, voxelSize, volumeSize);
	pointGrid->ComputePointGrid();

	// Compute ray-mesh hit points
	vector<vector<HitPoint>> XRayHitPts;
	ComputeCrossPoints( XRayHitPts );
	pointGrid->SetXAxisHitPoints( XRayHitPts );

	// Compute point state (i.e., in or out mesh)
	pointGrid->ComputePointState();
}


void LVolume::ComputeCrossPoints(vector<vector<HitPoint>> &XRayHitPts)
{
	XRayHitPts.clear();

	Vector3i gridDimen = pointGrid->gridDimen;

	///////////////////////////////////////////////////////////
	// Compute intersection points for x-axis aligned rays

	for( int j=0; j< gridDimen[1]; j++ )
	for( int k=0; k< gridDimen[2]; k++ )
	{
		// Construct each x-axis aligned ray
		Vector3f point = pointGrid->GetPointPosition(Vector3i(0,              j, k));
		Vector3f endPt = pointGrid->GetPointPosition(Vector3i(gridDimen[0]-1, j, k));

		// Find all the triangles that could intersect with the ray
		vector<Triangle*> testTris;
		for (int i=0; i<volumeSize[0]; i++)
		{
			int voxelY = j / (VOXEL_EDGE_SAMPLE_NUM-1); // Position of the intersection voxel
			int voxelZ = k / (VOXEL_EDGE_SAMPLE_NUM-1);
			int remaiY = j % (VOXEL_EDGE_SAMPLE_NUM-1); // Is the ray the shared edge between neighboring voxels
			int remaiZ = k % (VOXEL_EDGE_SAMPLE_NUM-1);

			// Case 1: The ray along x-axis is the shared edge between the neighboring voxels
			if ( remaiY ==0 || remaiZ == 0 )
			{
				UpdateTestTriangles(Vector3i(i, voxelY,   voxelZ  ), testTris);
				UpdateTestTriangles(Vector3i(i, voxelY-1, voxelZ  ), testTris);
				UpdateTestTriangles(Vector3i(i, voxelY,   voxelZ-1), testTris);
				UpdateTestTriangles(Vector3i(i, voxelY-1, voxelZ-1), testTris);
			}
			// Case 2: The ray along x-axis is within one voxel
			else
			{
				UpdateTestTriangles(Vector3i(i, voxelY,   voxelZ  ), testTris);
			}
		}

		// Perform intersection test between the ray and the triangeles
		vector<HitPoint> hitPoints = RayMeshIntersect(point, Vector3f(1,0,0), testTris);
		XRayHitPts.push_back( hitPoints );
	}
}

void LVolume::UpdateTestTriangles(Vector3i voxelPos, vector<Triangle*> &testTris)
{
	if ( voxelPos[1] < 0              )    voxelPos[1] = 0;
	if ( voxelPos[1] > volumeSize[1]-1)    voxelPos[1] = volumeSize[1]-1;
	if ( voxelPos[2] < 0              )    voxelPos[2] = 0;
	if ( voxelPos[2] > volumeSize[2]-1)    voxelPos[2] = volumeSize[2]-1;

	int index = GetVoxelIndex( voxelPos );
	if ( voxelGrid[index].state == VOXEL_CROSS_MESH )
	{
		vector<Triangle*> voxelTris = voxelGrid[index].GetVoxelOrigTris();
		for (int m=0; m<voxelTris.size(); m++)
		{
			testTris.push_back( voxelTris[m] );
		}
	}
}




//**************************************************************************************//
//                    Compute Shape Descriptor (using Single Feature)
//**************************************************************************************//

void LVolume::ComputeDescriptor_SurfArea(vector<float> &descVec)
{
	descVec.clear();

	for (int i=0; i<voxelGrid.size(); i++)
	{
		float surfArea = voxelGrid[i].ComputeSurfaceArea();

		descVec.push_back( surfArea );
	}
}

void LVolume::ComputeDescriptor_SurfPos(vector<float> &descVec)
{
	descVec.clear();

	for (int i=0; i<voxelGrid.size(); i++)
	{
		Vector3f surfCenter = voxelGrid[i].ComputeSurfaceCenter();

		descVec.push_back( surfCenter[0] );
		descVec.push_back( surfCenter[1] );
		descVec.push_back( surfCenter[2] );
	}
}

void LVolume::ComputeDescriptor_Normal(vector<float> &descVec, Vector3f keyPtNormal)
{
	descVec.clear();

	for (int i=0; i<voxelGrid.size(); i++)
	{
		Vector3f surfNormal = voxelGrid[i].ComputeSurfaceNormal();
		float value = surfNormal DOT keyPtNormal; // Note: the value could be negative
	
		descVec.push_back( value );
	}
}

void LVolume::ComputeDescriptor_Curvature(vector<float> &descVec)
{
	descVec.clear();

	for (int i=0; i<voxelGrid.size(); i++)
	{
		float surfCurvature = voxelGrid[i].ComputeSurfaceCurvature(); //Note: the value could be negative

		descVec.push_back( surfCurvature );
	}
}

void LVolume::ComputeDescriptor_ShapeVolume(vector<float> &descVec)
{
	ComputePointGrid(); // For estimating local shape inside each voxel

	descVec.clear();
	for (int i=0; i<voxelGrid.size(); i++)
	{
		voxelGrid[i].SetSamplePoints( pointGrid->GetPointsInsideVoxel(voxelGrid[i].pos) );
		float localVolume = voxelGrid[i].ComputeShapeVolume();

		descVec.push_back( localVolume );
	}
}

void LVolume::ComputeDescriptor_ShapePos(vector<float> &descVec)
{
	ComputePointGrid(); // For estimating local shape inside each voxel

	descVec.clear();
	for (int i=0; i<voxelGrid.size(); i++)
	{
		voxelGrid[i].SetSamplePoints( pointGrid->GetPointsInsideVoxel(voxelGrid[i].pos) );
		Vector3f shapeCenter = voxelGrid[i].ComputeShapeCenter();

		descVec.push_back( shapeCenter[0] );
		descVec.push_back( shapeCenter[1] );
		descVec.push_back( shapeCenter[2] );
	}
}

void LVolume::ComputeSGC(SGCentroid& feature)
{
//	ComputePointGrid(); // For estimating local shape inside each voxel

	feature.ClearSGCentroid();
	
	for (int i = 0; i<voxelGrid.size(); i++)
	{

		SGCBin thisbin(0, 0, 0, 0);

		Vector3f centroid(0,0,0);
		int cur_non_empty = 0;

		for (int j = 0; j < volumePoints.size(); j++)
		{
			if (voxelGrid[i].isPointIn(volumePoints.at(j).pos))
			{

				centroid += volumePoints.at(j).pos - voxelGrid[i].center + 0.5f * voxelGrid[i].size;
				cur_non_empty++;
			}
		}
		if (cur_non_empty > 0)
		{
			centroid /= (1.0*cur_non_empty);
			thisbin.setValue(centroid[0], centroid[1], centroid[2], cur_non_empty);
		}

		feature.FeatureDetails.push_back(thisbin);
	}
}


void LVolume::ComputeDescriptor_Color(vector<float> &descVec)
{
	descVec.clear();

	for (int i=0; i<voxelGrid.size(); i++)
	{
		Vector3f surfColor = voxelGrid[i].ComputeSurfaceColor();

		descVec.push_back( surfColor[0] );
		descVec.push_back( surfColor[1] );
		descVec.push_back( surfColor[2] );
	}
}




//**************************************************************************************//
//                    Compute Descriptor (using Feature Combinations)
//**************************************************************************************//

void LVolume::ComputeDescriptor_Area_SurfCen(vector<float> &descVec)
{
	descVec.clear();

	for (int i=0; i<voxelGrid.size(); i++)
	{
		float surfArea = voxelGrid[i].ComputeSurfaceArea();
		Vector3f surfCenter = voxelGrid[i].ComputeSurfaceCenter();

		descVec.push_back( surfArea );
		descVec.push_back( surfCenter[0] );
		descVec.push_back( surfCenter[1] );
		descVec.push_back( surfCenter[2] );
	}
}

void LVolume::ComputeDescriptor_Area_Normal(vector<float> &descVec, Vector3f keyPtNormal)
{
	descVec.clear();

	for (int i=0; i<voxelGrid.size(); i++)
	{
		float surfArea = voxelGrid[i].ComputeSurfaceArea();
		Vector3f surfNormal = voxelGrid[i].ComputeSurfaceNormal();
		float dotValue = surfNormal DOT keyPtNormal; // Note: the value could be negative

		descVec.push_back( surfArea );
		descVec.push_back( dotValue );
	}
}

void LVolume::ComputeDescriptor_Volume_VolCen(vector<float> &descVec)
{
	ComputePointGrid(); // For estimating local shape inside each voxel

	descVec.clear();
	for (int i=0; i<voxelGrid.size(); i++)
	{
		voxelGrid[i].SetSamplePoints( pointGrid->GetPointsInsideVoxel(voxelGrid[i].pos) );
		float localVolume    = voxelGrid[i].ComputeShapeVolume();
		Vector3f shapeCenter = voxelGrid[i].ComputeShapeCenter();

		descVec.push_back( localVolume );
		descVec.push_back( shapeCenter[0] );
		descVec.push_back( shapeCenter[1] );
		descVec.push_back( shapeCenter[2] );
	}
}

void LVolume::ComputeDescriptor_SurfCen_Normal(vector<float> &descVec, Vector3f keyPtNormal)
{
	descVec.clear();

	for (int i=0; i<voxelGrid.size(); i++)
	{
		Vector3f surfCenter = voxelGrid[i].ComputeSurfaceCenter();
		Vector3f surfNormal = voxelGrid[i].ComputeSurfaceNormal();
		float dotValue = surfNormal DOT keyPtNormal; // Note: the value could be negative

		descVec.push_back( surfCenter[0] );
		descVec.push_back( surfCenter[1] );
		descVec.push_back( surfCenter[2] );
		descVec.push_back( dotValue );
	}
}

void LVolume::ComputeDescriptor_Area_SurfCen_Normal(vector<float> &descVec, Vector3f keyPtNormal)
{
	descVec.clear();

	for (int i=0; i<voxelGrid.size(); i++)
	{
		float surfArea = voxelGrid[i].ComputeSurfaceArea();
		Vector3f surfCenter = voxelGrid[i].ComputeSurfaceCenter();
		Vector3f surfNormal = voxelGrid[i].ComputeSurfaceNormal();
		float dotValue = surfNormal DOT keyPtNormal; // Note: the value could be negative

		descVec.push_back( surfArea );
		descVec.push_back( surfCenter[0] );
		descVec.push_back( surfCenter[1] );
		descVec.push_back( surfCenter[2] );
		descVec.push_back( dotValue );
	}
}




//**************************************************************************************//
//                                  Rendering Stuff
//**************************************************************************************//

void LVolume::DrawLocalBBox(float lineWidth, Vector3f color)
{
	//if ( volumeTris.size() == 0 )
	//	return;

	DrawWireCuboid( volumeMinPt, volumeMaxPt, lineWidth, color );
}

void LVolume::DrawVoxelGrid(float lineWidth, vec color)
{
	if ( volumeTris.size() == 0 )
		return;

	for (int i=0; i<voxelGrid.size(); i++)
	{
		if( voxelGrid[i].state == VOXEL_CROSS_MESH  )
		{
			voxelGrid[i].DrawVoxel(lineWidth, color);
		}
	}
}


void LVolume::DrawLocalShape()
{
	if ( volumeTris.size() == 0 )
		return;

	for (int i=0; i<voxelGrid.size(); i++)
	{
		if( voxelGrid[i].state == VOXEL_CROSS_MESH )
		{
			voxelGrid[i].DrawVoxelTris();
			//break;
		}
	}
}

void LVolume::DrawLocalShapeWire(vec color, float width)
{
	if ( volumeTris.size() == 0 )
		return;

	for (int i=0; i<voxelGrid.size(); i++)
	{
		if( voxelGrid[i].state == VOXEL_CROSS_MESH )
		{
			voxelGrid[i].DrawVoxelTrisWire(color, width);
			//break;
		}
	}
}

/*void LVolume::DrawLocalShape()
{
	if ( volumeTris.size() == 0 )
		return;

	glEnable(GL_LIGHTING);

	//for (int i=0; i<sphereTris.size(); i++)
	for (int i=0; i<volumeTris.size(); i++)
	{
		//Triangle *tri = sphereTris[i];
		Triangle *tri = volumeTris[i];

		//if ( tri->isCross == false )
		//	continue;

		glBegin(GL_TRIANGLES);
		glNormal3fv( tri->normal );
		glVertex3fv( tri->v[0] );
		glVertex3fv( tri->v[1] );
		glVertex3fv( tri->v[2] );
		glEnd();
	}
}

void LVolume::DrawLocalShapeWire(vec color, float width)
{
	if ( volumeTris.size() == 0 )
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
	//for (int i=0; i<sphereTris.size(); i++)
	for (int i=0; i<volumeTris.size(); i++)
	{
		//Triangle *tri = sphereTris[i];
		Triangle *tri = volumeTris[i];

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
}*/




//**************************************************************************************//
//                                  Draw Point Grid
//**************************************************************************************//

void LVolume::DrawPointGrid(float pointSize, Vector3f pointColor)
{
	if ( pointGrid == NULL )
		return;

	pointGrid->DrawPointGrid(pointSize, pointColor);
}

void LVolume::DrawCrossPoints()
{
	if ( pointGrid == NULL )
		return;

	pointGrid->DrawCrossPoints();
}

void LVolume::DrawVolumeSamplePts(int pointSize, Vector3f color)
{
	for (int i=0; i<voxelGrid.size(); i++)
	{
		if( voxelGrid[i].state == VOXEL_CROSS_MESH  )
		{
			voxelGrid[i].DrawSamplePoints(pointSize, color);
		}
	}
}