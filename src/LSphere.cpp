///////////////////////////////////////////////////////////////
//
// LSphere.cpp
//
//   Build Local Sphere for Describing Local Shape
//
// by Song Peng ( songpeng@ustc.edu.cn )
// 
// 14/May/2015
//
///////////////////////////////////////////////////////////////

#include "Controls.h"
#include "Helper.h"
#include "math3D.h"
#include "LSBin.h"
#include "LSphere.h"


//**************************************************************************************//
//                                   Initialization
//**************************************************************************************//

LSphere::LSphere()
{

}

LSphere::~LSphere()
{
	ClearLSphere();
}

void LSphere::ClearLSphere()
{

}

void LSphere::InitLSphere(Vector3i _gridSize, float _radius, vector<Triangle*> _shapeTris, vector<Point> _shapePts)
{
	gridSize  = _gridSize;
	radius    = _radius;
	shapeTris = _shapeTris;
	shapePts  = _shapePts;
}




//**************************************************************************************//
//                            Preprocess Local Shape
//**************************************************************************************//

void LSphere::ComputPointDensity(float localRadius)
{
	for (int i=0; i<shapePts.size(); i++)
	{
		// TODO: speed up this step by building kd-tree
		vector<Point> neiborPts = GetNeighborPoints(shapePts[i].pos, radius, shapePts);

		shapePts[i].curv = neiborPts.size(); 
	}
}

vector<Point> LSphere::GetNeighborPoints(Vector3f center, float radius, vector<Point> pointList)
{
	vector<Point> neiborPts;

	for (int i=0; i<pointList.size(); i++)
	{
		Point currPt = pointList[i];

		float dist = len(currPt.pos-center);

		if ( dist < radius )
		{
			neiborPts.push_back( currPt );
		}
	}

	return neiborPts;
}




//**************************************************************************************//
//                                Compute Voxel Grid
//**************************************************************************************//

void LSphere::ComputeBinGrid()
{
	binGrid.clear();

	float rmin = 0.1*radius;
	float rmax = radius;

	for( int k=0; k< gridSize[2]; k++ )
	for( int j=0; j< gridSize[1]; j++ )
	for( int i=0; i< gridSize[0]; i++ )
	{
		LSBin bin;
		bin.state = BIN_OUT_MESH;

		bin.pos[0] = i; 
		bin.pos[1] = j; 
		bin.pos[2] = k; 

		// Partition the radius linearly
		//bin.minRadius = (   i /(float)gridSize[0]) * radius;
		//bin.maxRadius = ((i+1)/(float)gridSize[0]) * radius;

		// Partition the radius logarithmically
		bin.minRadius = exp( log(rmin)+ (   i /(float)gridSize[0]) * log(rmax/rmin) ); 
		bin.maxRadius = exp( log(rmin)+ ((i+1)/(float)gridSize[0]) * log(rmax/rmin) ); 

		bin.minTheta  = (   j /(float)gridSize[1] ) * 180;
		bin.maxTheta  = ((j+1)/(float)gridSize[1] ) * 180;
		bin.minAlpha  = (   k /(float)gridSize[2] ) * 360;
		bin.maxAlpha  = ((k+1)/(float)gridSize[2] ) * 360;

		float maxSphereVol = 4.0/3.0*M_PI * pow(bin.maxRadius, 3);
		float minSphereVol = 4.0/3.0*M_PI * pow(bin.minRadius, 3);
		bin.volume = (maxSphereVol - minSphereVol) / (float)(gridSize[1]*gridSize[2]);
		bin.volume = bin.volume / pow(radius,3); // Normalize the bin volume

		bin.ComputeBinCenter();

		binGrid.push_back( bin );
	}

	//printf("bin grid size: %d \n", binGrid.size());

	//for (int i=0; i<binGrid.size(); i++)
	//{
	//	binGrid[i].PrintBin();
	//}
}

void LSphere::ComputBinState()
{
	if ( shapeTris.size() ==0 )
		return;

	for (int i=0; i<binGrid.size(); i++)
	{
		binGrid[i].GetBinTriangles( shapeTris );
		//binGrid[i].GetBinPoints( shapeTris );
		binGrid[i].GetBinPoints( shapePts );

		if ( binGrid[i].binTris.size() != 0 ||
			 binGrid[i].binPoints.size() != 0 )
		{
			binGrid[i].state = BIN_CROSS_MESH;
		}
	}
}




//**************************************************************************************//
//                     Compute Shape Descriptor (3DSC and SHOT)
//**************************************************************************************//

void LSphere::ComputeDescriptor_3DSC(vector<float> &descVec)
{
	descVec.clear();

	for (int i=0; i<binGrid.size(); i++)
	{
		float binValue = binGrid[i].ComputeBin3DSCValue();

		descVec.push_back( binValue );
	}

	//printf("desc size: %d \n", descVec.size());
	//for (int i=0; i<descVec.size(); i++)
	//{
	//	printf(" %.2f ", descVec[i]);
	//}
	//printf("\n\n");
}

void LSphere::ComputeDescriptor_SHOT(vector<float> &descVec, Vector3f keyPtNormal)
{
	descVec.clear();

	for (int i=0; i<binGrid.size(); i++)
	{
		vector<int> binValues = binGrid[i].ComputeBinSHOTValues( keyPtNormal );

		for (int j=0; j<binValues.size(); j++)
		{
			float descValue = binValues[j] / 100.0;  // Note: to make the value around 1.0
			descVec.push_back( descValue );
		}
		//descVec.push_back( binValue );
	}

	//printf(" desc size: %d \n", descVec.size());

	//for (int i=0; i<descVec.size(); i++)
	//{
	//	printf(" %.2f ", descVec[i]);
	//}
	//printf("\n\n");
}




//**************************************************************************************//
//                                  Rendering Stuff
//**************************************************************************************//

void LSphere::DrawLocalSphere(float lineWidth, Vector3f color)
{
	glDisable( GL_LIGHTING );
	glPointSize( 2.0 );
	glLineWidth( 3.0 );

	glColor3f(0.9, 0.3, 0.3);

	//for (int i=0; i<180; i=i+5)
	//{
	//	for (int j=0; j<360; j=j+5)
	//	{
	//		Vector3f point = sphBin.GetCartesianCoord(radius, (float)i, (float)j);

	//		glBegin(GL_POINTS);
	//		glVertex3f(point[0], point[1], point[2]);
	//		glEnd();
	//	}
	//}

	//printf("shapePt Num: %d \n", shapePts.size());

	/*if ( shapePts.size() >400 )
	{
		int index = 300;

		Point currPt = shapePts[index];
		float radius = 0.02;

		glPointSize(16.0);
		glColor3f(0.3,0.5,0.8);

		glBegin(GL_POINTS);
		glVertex3f(	currPt.pos[0], currPt.pos[1], currPt.pos[2] );
		glEnd();

		DrawWireSphere(currPt.pos, radius, Vector3f(0.3,0.6,0.3));

		vector<Point> tempPts = GetNeighborPoints(currPt.pos, radius, shapePts);

		glPointSize(18.0);
		glColor3f(0.1,0.9,0.1);
		for (int i=0; i<tempPts.size(); i++)
		{
			glBegin(GL_POINTS);
			glVertex3f(	tempPts[i].pos[0], tempPts[i].pos[1], tempPts[i].pos[2] );
			glEnd();
		}

		glPointSize(1.0);
	}*/

	for (int i=0; i<binGrid.size(); i++)
	{
		if ( binGrid[i].state == BIN_CROSS_MESH )
		{
			binGrid[i].DrawSphericalBin( 3, Vector3f(0.2,0.9,0.2) );
		}
		else
		{
			binGrid[i].DrawSphericalBin( 3, Vector3f(0.5,0.5,0.5) );
		}
	}

	//for (int i=0; i<binGrid.size(); i++)
	//{
	//	//binGrid[i].DrawBinTriangles( Vector3f(0.8,0.3,0.9));
	//	binGrid[i].DrawBinPoints( Vector3f(0.8,0.3,0.9));
	//}

	glEnable(GL_LIGHTING);
	glPointSize(1.0);
	glLineWidth(1.0);
}

void LSphere::DrawSpheBinGrid(float lineWidth, vec color)
{
	//if ( volumeTris.size() == 0 )
	//	return;

	//for (int i=0; i<voxelGrid.size(); i++)
	//{
	//	if( voxelGrid[i].state == VOXEL_CROSS_MESH  )
	//	{
	//		voxelGrid[i].DrawVoxel(lineWidth, color);
	//	}
	//}
}

void LSphere::DrawLocalShapeWire(vec color, float width)
{
	//if ( volumeTris.size() == 0 )
	//	return;

	//for (int i=0; i<voxelGrid.size(); i++)
	//{
	//	if( voxelGrid[i].state == VOXEL_CROSS_MESH )
	//	{
	//		voxelGrid[i].DrawVoxelTrisWire(color, width);
	//		//break;
	//	}
	//}
}