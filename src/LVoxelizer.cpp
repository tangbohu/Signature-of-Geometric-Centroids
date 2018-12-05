///////////////////////////////////////////////////////////////
//
// LVoxelizer.cpp
//
//   Local Voxelizer Class
//
// by Song Peng ( songpeng@ustc.edu.cn )
// 
// 13/May/2015
//
///////////////////////////////////////////////////////////////

#include "Controls.h"
#include <vector>
#include "math3D.h"
#include <GL/glut.h>
#include "LVoxelizer.h"


//**************************************************************************************//
//                                   Initialization
//**************************************************************************************//

LVoxelizer::LVoxelizer()
{

}

LVoxelizer::~LVoxelizer()
{

}

void LVoxelizer::InitLVoxelizer(Point _keyPoint, float _radius, vector<Triangle*> _LRFTris, vector<Triangle*> _shapeTris, vector<Point> _LRFPts, vector<Point> _shapePts)
{
	keyPoint  = _keyPoint;
	radius    = _radius;

	LRFTris   = _LRFTris;
	shapeTris = _shapeTris;
	LRFPts = _LRFPts;
	shapePts  = _shapePts;
}




//**************************************************************************************//
//                          Compute Descriptor using LSphere
//**************************************************************************************//

Descriptor* LVoxelizer::ComputeDescriptor_Spherical(int featureTypeID, bool isDebug)
{
	Descriptor *descriptor;

	if ( !isDebug )
	{
		descriptor = new Descriptor(); 	// Avoid memory release issue
	}

	if ( shapeTris.size() == 0 )
		return descriptor;


	///////////////////////////////////////////////////
	// 1. Build LRF using the local shape (i.e., LRFTris)

	keyPtLRF.InitLRF(keyPoint, radius, LRFTris, shapePts);
	keyPtLRF.BuildLocalFrame();


	///////////////////////////////////////////////////
	// 2. Transform local shape (i.e., shapeTris) into LRF

	vector<Triangle*> shapeTris = TransformTriangles( LRFTris );
	vector<Point> newShapePts = TransformPoints( shapePts );


	///////////////////////////////////////////////////
	// 3. Create local spherical grid based on an LRF

	Vector3i sphereDimen;
#ifdef USE_DESCRIPTOR_3DSC
	sphereDimen = Vector3i(5,4,8);
#endif
#ifdef USE_DESCRIPTOR_SHOT
	sphereDimen = Vector3i(2,2,8);
#endif

	localSphere.InitLSphere(sphereDimen, radius, shapeTris, newShapePts);
	localSphere.ComputeBinGrid();
#ifdef USE_DESCRIPTOR_3DSC
	float localPtRadius = 0.01;
	localSphere.ComputPointDensity( localPtRadius );
#endif
	localSphere.ComputBinState();


	///////////////////////////////////////////////////
	// 4. Create local spherical grid based on an LRF

#ifdef USE_DESCRIPTOR_3DSC
	localSphere.ComputeDescriptor_3DSC( descVector );
#endif 

#ifdef USE_DESCRIPTOR_SHOT
	localSphere.ComputeDescriptor_SHOT( descVector, keyPoint.nor );
#endif


	///////////////////////////////////////////////////
	// 5. Return computer descriptor

	if ( !isDebug )
	{
		descriptor->descVector = descVector;
		descriptor->keyPtLRF   = keyPtLRF;
	}

	return descriptor;
}

vector<Point> LVoxelizer::TransformPoints(vector<Point> _shapePoints)
{
	vector<Point> newShapePts;

	// Calculate the inverse of LRF matrix
	double inveLRFMat[16], tranInveLRFMat[16];
	memcpy(inveLRFMat, keyPtLRF.LRFMatrix, sizeof(double)*16);
	if(invertMat(inveLRFMat, 4)==1)  printf("Inverse Matrix Error \n");
	memcpy(tranInveLRFMat, inveLRFMat, sizeof(double)*16);
	transposeMat(tranInveLRFMat, 4);

	// Note: only transform triangles that are a bit larger than the local volume
	for (int i=0; i<_shapePoints.size(); i++)
	{
		double origVertex[3];
		origVertex[0] = _shapePoints[i].pos[0];
		origVertex[1] = _shapePoints[i].pos[1];			
		origVertex[2] = _shapePoints[i].pos[2];

		double origNormal[3];
		origNormal[0] = _shapePoints[i].nor[0];
		origNormal[1] = _shapePoints[i].nor[1];			
		origNormal[2] = _shapePoints[i].nor[2];

		Vector3f currVertex;
		currVertex[0] = dot3D(origVertex, tranInveLRFMat  ) + inveLRFMat[12];
		currVertex[1] = dot3D(origVertex, tranInveLRFMat+4) + inveLRFMat[13];
		currVertex[2] = dot3D(origVertex, tranInveLRFMat+8) + inveLRFMat[14];

		Vector3f currNormal;
		currNormal[0] = dot3D(origNormal, tranInveLRFMat  );
		currNormal[1] = dot3D(origNormal, tranInveLRFMat+4);
		currNormal[2] = dot3D(origNormal, tranInveLRFMat+8);

		Point point;
		point.pos  = currVertex;
		point.nor  = currNormal;
		point.curv = _shapePoints[i].curv; // To represent point density around the point

		newShapePts.push_back( point );
	}

	return newShapePts;
}




//**************************************************************************************//
//                           Compute Descriptor using LVolume
//**************************************************************************************//


SGCentroid * LVoxelizer::ComputeSGC()
{
	SGCentroid *descriptor= new SGCentroid(); 	// Avoid memory release issue
	

	///////////////////////////////////////////////////
	// 1. Build LRF using the local shape (i.e., LRFTris)

	keyPtLRF.InitLRF(keyPoint, radius, LRFTris, LRFPts);
	keyPtLRF.BuildLocalFrame(true);


	///////////////////////////////////////////////////
	// 2. Transform local shape (i.e., shapeTris) into LRF

	vector<Point> spherePoint = TransformPoints(shapePts);


	///////////////////////////////////////////////////
	// 3. Create local voxelization based on an LRF

	localVolume.InitLVolumeWithTransformPoints(LOCAL_VOLUME_DIMEN, vec(2 * radius, 2 * radius, 2 * radius), spherePoint);
	localVolume.ComputeVoxelGrid();
//	localVolume.ComputeVoxelState();




	///////////////////////////////////////////////////
	// 4. Construct descriptor based on the computed local volume

	localVolume.ComputeSGC(*descriptor);

	descriptor->keyPtLRF = keyPtLRF;


	///////////////////////////////////////////////////
	// 5. Return computer descriptor

	return descriptor;

}



Descriptor* LVoxelizer::ComputeDescriptor(int featureTypeID, bool isDebug)
{
	Descriptor *descriptor;

	if ( !isDebug )
	{
		descriptor = new Descriptor(); 	// Avoid memory release issue
	}

	if ( shapeTris.size() == 0 )
		return descriptor;


	///////////////////////////////////////////////////
    // 1. Build LRF using the local shape (i.e., LRFTris)

	keyPtLRF.InitLRF(keyPoint, radius, LRFTris, shapePts);
	keyPtLRF.BuildLocalFrame();


	///////////////////////////////////////////////////
	// 2. Transform local shape (i.e., shapeTris) into LRF

	vector<Triangle*> sphereTris = TransformTriangles( shapeTris );


	///////////////////////////////////////////////////
	// 3. Create local voxelization based on an LRF

	localVolume.InitLVolume(LOCAL_VOLUME_DIMEN, vec(2*radius, 2*radius, 2*radius), sphereTris);
	localVolume.ComputeVoxelGrid();
	localVolume.ComputeVoxelState();

	// TODO: clip local mesh only when it is very dense
	if ( featureTypeID == FEATURE_SURFACE_CURVATURE ||
		 featureTypeID == FEATURE_SHAPE_VOLUME ||
		 featureTypeID == FEATURE_SHAPE_CENTER )
	{
		;
	}
	else
	{
		localVolume.ClipMeshInVoxels();
	}


	///////////////////////////////////////////////////
	// 4. Construct descriptor based on the computed local volume

	if      ( featureTypeID == FEATURE_SURFACE_AREA )             localVolume.ComputeDescriptor_SurfArea( descVector );
	else if ( featureTypeID == FEATURE_SURFACE_CRNTER )           localVolume.ComputeDescriptor_SurfPos( descVector );
	else if ( featureTypeID == FEATURE_SURFACE_NORMAL )           localVolume.ComputeDescriptor_Normal( descVector, keyPoint.nor );
	else if ( featureTypeID == FEATURE_SURFACE_CURVATURE )        localVolume.ComputeDescriptor_Curvature( descVector );
	else if ( featureTypeID == FEATURE_SHAPE_VOLUME )             localVolume.ComputeDescriptor_ShapeVolume( descVector ); // Note: quite slow
	else if ( featureTypeID == FEATURE_SHAPE_CENTER )             localVolume.ComputeDescriptor_ShapePos( descVector ); // Note: quite slow
	
	else if ( featureTypeID == FEATURE_COMB_AREA_SURFCEN )        localVolume.ComputeDescriptor_Area_SurfCen( descVector );
	else if ( featureTypeID == FEATURE_COMB_AREA_NORMAL )         localVolume.ComputeDescriptor_Area_Normal( descVector, keyPoint.nor );
	else if ( featureTypeID == FEATURE_COMB_VOLUME_VOLCEN )       localVolume.ComputeDescriptor_Volume_VolCen( descVector );
	//else if ( featureTypeID == FEATURE_COMB_SURFCEN_NORMAL )      localVolume.ComputeDescriptor_SurfCen_Normal( descVector, keyPoint.nor );
	//else if ( featureTypeID == FEATURE_COMB_AREA_SURFCEN_NORM )   localVolume.ComputeDescriptor_Area_SurfCen_Normal( descVector, keyPoint.nor );


	///////////////////////////////////////////////////
	// 5. Return computer descriptor

	if ( !isDebug )
	{
		descriptor->descVector = descVector;
		descriptor->keyPtLRF   = keyPtLRF;
	}

	return descriptor;
}

vector<Triangle*> LVoxelizer::TransformTriangles(vector<Triangle*> _shapeTris)
{	
	vector<Triangle*> newShapeTris;

	// Calculate the inverse of LRF matrix
	double inveLRFMat[16], tranInveLRFMat[16];
	memcpy(inveLRFMat, keyPtLRF.LRFMatrix, sizeof(double)*16);
	if(invertMat(inveLRFMat, 4)==1)  printf("Inverse Matrix Error \n");
	memcpy(tranInveLRFMat, inveLRFMat, sizeof(double)*16);
	transposeMat(tranInveLRFMat, 4);

	// Note: only transform triangles that are a bit larger than the local volume
	for (int i=0; i<_shapeTris.size(); i++)
	{
		Triangle* tri = new Triangle();

		// Transform the vertices of the triangle
		for (int j=0; j<3; j++)
		{
			double origVertex[3];
			origVertex[0] = _shapeTris[i]->v[j][0];
			origVertex[1] = _shapeTris[i]->v[j][1];			
			origVertex[2] = _shapeTris[i]->v[j][2];

			vec currVertex;
			currVertex[0] = dot3D(origVertex, tranInveLRFMat  ) + inveLRFMat[12];
			currVertex[1] = dot3D(origVertex, tranInveLRFMat+4) + inveLRFMat[13];
			currVertex[2] = dot3D(origVertex, tranInveLRFMat+8) + inveLRFMat[14];

			tri->v[j]     = currVertex;

			tri->color[j] = _shapeTris[i]->color[j];
			tri->curv[j]  = _shapeTris[i]->curv[j];
		}

		// Triangle size remains the same
		tri->area = _shapeTris[i]->area;

		vec normal = (tri->v[1] - tri->v[0]) CROSS (tri->v[2] - tri->v[0]);
		tri->normal = normal / len(normal);

		newShapeTris.push_back( tri );
	}

	return newShapeTris;
}


