#include "SGCentroid.h"


SGCentroid::SGCentroid()
{
	ClearSGCentroid();
}

SGCentroid::~SGCentroid()
{
	ClearSGCentroid();

}

void SGCentroid::ClearSGCentroid()
{
	FeatureDetails.clear();
}

void SGCentroid::InsertBin(SGCBin bin)
{
	FeatureDetails.push_back(bin);
}

double SGCentroid::getDistanceWithOthers(SGCentroid& that)
{



	double dist = 0;
	double tempdist = 0;

	for (int i = 0; i < that.FeatureDetails.size(); i++)
	{
		double mindist = 0;
		mindist = FeatureDetails.at(i).distWith(that.FeatureDetails.at(i));

		dist += mindist;
	}
	return dist;
}


//Descriptor* SGCentroid::ComputeDescriptor(int featureTypeID, bool isDebug)
//{
//	Descriptor *descriptor;
//
//	if (!isDebug)
//	{
//		descriptor = new Descriptor(); 	// Avoid memory release issue
//	}
//
//	if (shapeTris.size() == 0)
//		return descriptor;
//
//
//	///////////////////////////////////////////////////
//	// 1. Build LRF using the local shape (i.e., LRFTris)
//
//	keyPtLRF.InitLRF(keyPoint, radius, LRFTris);
//	keyPtLRF.BuildLocalFrame();
//
//
//	///////////////////////////////////////////////////
//	// 2. Transform local shape (i.e., shapeTris) into LRF
//
//	vector<Triangle*> sphereTris = TransformTriangles(shapeTris);
//
//
//	///////////////////////////////////////////////////
//	// 3. Create local voxelization based on an LRF
//
//	localVolume.InitLVolume(LOCAL_VOLUME_DIMEN, vec(2 * radius, 2 * radius, 2 * radius), sphereTris);
//	localVolume.ComputeVoxelGrid();
//	localVolume.ComputeVoxelState();
//
//	// TODO: clip local mesh only when it is very dense
//	if (featureTypeID == FEATURE_SURFACE_CURVATURE ||
//		featureTypeID == FEATURE_SHAPE_VOLUME ||
//		featureTypeID == FEATURE_SHAPE_CENTER)
//	{
//		;
//	}
//	else
//	{
//		localVolume.ClipMeshInVoxels();
//	}
//
//
//	///////////////////////////////////////////////////
//	// 4. Construct descriptor based on the computed local volume
//
//	if (featureTypeID == FEATURE_SURFACE_AREA)             localVolume.ComputeDescriptor_SurfArea(descVector);
//	else if (featureTypeID == FEATURE_SURFACE_CRNTER)           localVolume.ComputeDescriptor_SurfPos(descVector);
//	else if (featureTypeID == FEATURE_SURFACE_NORMAL)           localVolume.ComputeDescriptor_Normal(descVector, keyPoint.nor);
//	else if (featureTypeID == FEATURE_SURFACE_CURVATURE)        localVolume.ComputeDescriptor_Curvature(descVector);
//	else if (featureTypeID == FEATURE_SHAPE_VOLUME)             localVolume.ComputeDescriptor_ShapeVolume(descVector); // Note: quite slow
//	else if (featureTypeID == FEATURE_SHAPE_CENTER)             localVolume.ComputeDescriptor_ShapePos(descVector); // Note: quite slow
//
//	else if (featureTypeID == FEATURE_COMB_AREA_SURFCEN)        localVolume.ComputeDescriptor_Area_SurfCen(descVector);
//	else if (featureTypeID == FEATURE_COMB_AREA_NORMAL)         localVolume.ComputeDescriptor_Area_Normal(descVector, keyPoint.nor);
//	else if (featureTypeID == FEATURE_COMB_VOLUME_VOLCEN)       localVolume.ComputeDescriptor_Volume_VolCen(descVector);
//	//else if ( featureTypeID == FEATURE_COMB_SURFCEN_NORMAL )      localVolume.ComputeDescriptor_SurfCen_Normal( descVector, keyPoint.nor );
//	//else if ( featureTypeID == FEATURE_COMB_AREA_SURFCEN_NORM )   localVolume.ComputeDescriptor_Area_SurfCen_Normal( descVector, keyPoint.nor );
//
//
//	///////////////////////////////////////////////////
//	// 5. Return computer descriptor
//
//	if (!isDebug)
//	{
//		descriptor->descVector = descVector;
//		descriptor->keyPtLRF = keyPtLRF;
//	}
//
//	return descriptor;
//}
