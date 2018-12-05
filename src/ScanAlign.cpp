///////////////////////////////////////////////////////////////
//
// ScanAlign.cpp
//
//   Aligning Two scanList Class
//
// by Song Peng ( songpeng@ustc.edu.cn )
// 
// 11/May/2015
//
//
///////////////////////////////////////////////////////////////

#include <vector>
#include "math3D.h"
#include <GL/glut.h>
#include "LRF.h"
#include "Helper.h"
#include "ScanAlign.h"


//**************************************************************************************//
//                                   Initialization
//**************************************************************************************//

ScanAlign::ScanAlign()
{
	::identity(alignMatrix);
}

ScanAlign::~ScanAlign()
{

}

void ScanAlign::InitScanAlign(LRF _modelLRF, LRF _dataLRF)
{
	modelLRF  =  _modelLRF;
	dataLRF   =  _dataLRF;
}




//**************************************************************************************//
//                              Compute Alignment Matrix
//**************************************************************************************//

void ScanAlign::ComputeAlignMatrix()
{
	// Calculate the inverse world matrix
	double tranDataLRFMat[16], inveDataLRFMat[16], tranInveDataLRFMat[16];
	memcpy(inveDataLRFMat, dataLRF.LRFMatrix, sizeof(double)*16);
	if(invertMat(inveDataLRFMat, 4)==1)  printf("Inverse Matrix Error \n");

	// Rotate the object
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadMatrixd( modelLRF.LRFMatrix  );
	glMultMatrixd( inveDataLRFMat );
	glGetDoublev(GL_MODELVIEW_MATRIX, alignMatrix);
	glPopMatrix();
}




//**************************************************************************************//
//                              Evaluate Alignment Matrix
//**************************************************************************************//

void ScanAlign::CompareAlignMatrix(double alignMat[16], double realMat[16], float &rotError, float &tranError)
{
	///////////////////////////////////////////////////////////////////
	// 1. Compute rotational error

	double inveRealMat[16];
	memcpy(inveRealMat, realMat, sizeof(double)*16);
	if(invert4by4(inveRealMat) != 1)  printf("Inverse Matrix Error \n");

	double outMat[16];
	MultiplyMatrix(alignMat, inveRealMat, outMat);

	double outTrace = outMat[0] + outMat[5] + outMat[10];
	rotError = acos(0.5*(outTrace-1)) * 180.0 / M_PI;


	///////////////////////////////////////////////////////////////////
	// 2. Compute translational error

	Vector3f alignVec = Vector3f(alignMat[12], alignMat[13], alignMat[14]);
	Vector3f realVec  = Vector3f(realMat[12],  realMat[13],  realMat[14] );

	tranError = len(alignVec-realVec);

	printf("=> rotError: %.2f   tranError: %.2f \n", rotError, tranError);
}
