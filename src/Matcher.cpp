///////////////////////////////////////////////////////////////
//
// Matcher.cpp
//
//   Match Model and Data Meshes
//
// by Song Peng ( songpeng@ustc.edu.cn )
// 
// 11/May/2015
//
//=======================================================
// TODO (optional)
// 1. Mesh point sampling
//    - Feature point selection 
//    - Refine sampling algorithm to support a larger number of sample points
// 2. Local shape descriptor
//    - Speed up the computation of LRF and descriptor
//=======================================================
//
///////////////////////////////////////////////////////////////

#include <windows.h>
#include "time.h"
#include <GL/glut.h>
#include "KDtree.h"
#include "math3D.h"
#include "Helper.h"
#include "HelpDefines.h"
#include "Controls.h"
#include "Scan.h"
#include "LVoxelizer.h"
#include "Model.h"
#include "Matcher.h"
#include "PoissonRecon.h"
#include "Experiment.h"


#define SCAN_FILE_POSTFIX_LENGTH         6

// Color scheme table
Vector3f colorTable[10] = { 
	Vector3f( 0.2, 0.8, 0.8 ),   //  1: Cyan
	Vector3f( 0.6, 0.6, 0.7 ),   //  2: Gray
	Vector3f( 0.8, 0.8, 0.5 ),   //  3: Yellow
	//Vector3f( 0.57, 0.74, 0.88 ),   //  3: Gray (for the overview figure)
	Vector3f( 0.9, 0.3, 0.3 ),   //  4: Red
	Vector3f( 0.3, 0.3, 0.9 ),   //  5: Blue
	Vector3f( 0.9, 0.5, 0.2 ),   //  6: Orange
	Vector3f( 0.7, 0.2, 0.9 ),   //  7: Purple
	Vector3f( 0.5, 0.3, 0.2 ),   //  8: Brown
	Vector3f( 0.9, 0.2, 0.6 ),   //  9: Pink
	Vector3f( 0.9, 0.6, 0.5 ) }; // 10:LightSalmon


//**************************************************************************************//
//                                   Initialization
//**************************************************************************************//

Matcher::Matcher()
{
	scaleFactor = 0.01;

	::identity(worldMatrix);
	::identity(worldAxesMatrix);

	tagtScanID = -1;
}

Matcher::~Matcher()
{
	for (int i=0; i<scanList.size(); i++)
	{
		delete scanList[i];
	}
}

void Matcher::ClearScans()
{
	for (int i=0; i<scanList.size(); i++)
	{
		delete scanList[i];
	}

	scanList.clear();
}

void Matcher::ResetScene()
{
	InitWorldMatrix(INIT_WORLD_POSITION, INIT_WORLD_ROT_AXIS, INIT_WORLD_ROT_ANGLE);
	InitWorldAxesMatrix(INIT_WORLD_AXES_POSITION, INIT_WORLD_AXES_ROT_AXIS, INIT_WORLD_AXES_ROT_ANGLE);

	for (int i=0; i<scanList.size(); i++)
	{
		scanList[i]->ResetPose();
	}
}




//**************************************************************************************//
//                               World and Scan Matrix
//**************************************************************************************//

void Matcher::InitWorldMatrix(Vector3f position, Vector3f rotAxis, float rotAngle)
{
	glPushMatrix();   
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(position[0], position[1], position[2]);  
	glRotatef(rotAngle, rotAxis[0], rotAxis[1], rotAxis[2]);
	glGetDoublev(GL_MODELVIEW_MATRIX, worldMatrix);
	glPopMatrix();
}

void Matcher::InitWorldAxesMatrix(Vector3f position, Vector3f rotAxis, float rotAngle)
{
	glPushMatrix();   
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(position[0], position[1], position[2]);  
	glRotatef(rotAngle, rotAxis[0], rotAxis[1], rotAxis[2]);
	glGetDoublev(GL_MODELVIEW_MATRIX, worldAxesMatrix);
	glPopMatrix();
}

void Matcher::UpdateGlobalMatrix()
{
	for (int i=0; i<scanList.size(); i++)
	{
		GetGlobalMatrix(scanList[i]->localMatrix,   scanList[i]->globalMatrix);
		GetGlobalMatrix(scanList[i]->alignLocalMat, scanList[i]->alignGlobalMat);
	}
}

void Matcher::GetGlobalMatrix(double *localMat, double *globalMat)
{
	// Transpose world matrix
	double tranWorldMat[16];
	memcpy(tranWorldMat, worldMatrix, sizeof(double)*16);
	transposeMat(tranWorldMat,4);

	// Calculate object's global matrix: GlobalMatrix = WorldMatrix * LocalMatrix
	double objTranMat[16], resultMat[16];   
	memcpy(objTranMat, localMat, sizeof(double)*16);
	transposeMat(objTranMat,4);
	multMat(tranWorldMat, objTranMat,resultMat,4);
	transposeMat(resultMat, 4);
	memcpy(globalMat, resultMat, sizeof(double)*16);
}




//**************************************************************************************//
//                                  Load Scan Meshes
//**************************************************************************************//

void Matcher::LoadScans(char *scanFilePath)
{
	ClearScans();

	vector<string> scanFileNames = GetScanFileNames( scanFilePath );
	for (int i=0; i<scanFileNames.size(); i++)
	{
		Scan *scan = new Scan();

		// Load 3D model of the scan
		scan->LoadMesh( scanFileNames[i].c_str() );
		scan->InitPose();

		scan->scanID = scanList.size();
	
		if ( i == 0 )
		{
			scaleFactor = scan->GetScaleFactor();
			printf("dimension scale factor: %.2f \n", scaleFactor);
			//printf("curvature scale factor: %.2f \n", curvScale);
		}

		scan->ProcessMesh( Vector3f(0.1,0.7,0.1), scaleFactor );
		scanList.push_back( scan );
	}

	SetScanModelsMTL();
}

void Matcher::SetScanModelsMTL()
{
	if ( scanList.size() == 0 )
		return;

	// Set color for pieceIDs from [0, 9]
	for (int i=0; i<_MIN((int)scanList.size(), 10); i++)
	{
		for (int j=0; j<3; j++)
		{
			scanList[i]->mtlDiffuse[j]  = colorTable[i][j];
			scanList[i]->mtlSpecular[j] = colorTable[i][j];
		}
	}

	// Set color for pieceIDs from [10, K-1]
	int colorNum = 10; // These pieces cannot have similar color to the pieces with index less than colorNum
	if( scanList.size() > 11 )
	{
		for (int i=10; i<scanList.size(); i++)
		{
			Vector3f randColor = GenerateRandomMTL( colorTable, colorNum );
			scanList[i]->SetDiffuseMTL( randColor );
		}
	}
}

Vector3f Matcher::GenerateRandomMTL(Vector3f existColors[10], int colorNum)
{
	const float colorDistThres = 0.25;
	const float colorSatuThres = 0.80;

	while( true ) 
	{
		Vector3f randColor;
		randColor[0] = rand()/(RAND_MAX+1.0);
		randColor[1] = rand()/(RAND_MAX+1.0);
		randColor[2] = rand()/(RAND_MAX+1.0);

		float maxValue = _MAX(randColor[0], _MAX(randColor[1], randColor[2]));
		float minValue = _MIN(randColor[0], _MIN(randColor[1], randColor[2]));
		float saturation = (maxValue-minValue) / maxValue;

		float minColorDist = GetMinColorDistance(existColors, colorNum, randColor);	

		if ( minColorDist > colorDistThres && 
			saturation   < colorSatuThres )
		{
			return randColor;
		}
	}
}

float Matcher::GetMinColorDistance(Vector3f existColors[10], int colorNum, Vector3f randColor)
{
	float minColorDist = MAX_FLOAT;

	for (int i=0; i<colorNum; i++)
	{
		float tempDist = dist(existColors[i], randColor);
		if ( tempDist < minColorDist)
		{
			minColorDist = tempDist;
		}
	}

	return minColorDist;
}




//**************************************************************************************//
//                               Get Scan File Names
//**************************************************************************************//

vector<string> Matcher::GetScanFileNames(char *scanFilePath)
{
	// Find the first scan file 
	WIN32_FIND_DATA FindFileData;
	string firstFileName = GetFirstFileName( scanFilePath );
	HANDLE hFind = FindFirstFile(firstFileName.c_str(), &FindFileData);

	// Get all scan file names in the folder
	vector<string> scanFileNames;
	do
	{
		string sScanFilePath  = string( FindFileData.cFileName );
		scanFileNames.push_back( sScanFilePath );
	}
	while (FindNextFile(hFind, &FindFileData));

	// Close handle
	FindClose(hFind);

	// Sort scan file names according to their postfix number
	SortScanFileNames( scanFileNames, true );

	return scanFileNames;
}

string Matcher::GetFirstFileName(const char scanFilePath[])
{
	string sScanFilePath  = string( scanFilePath );

	// Subtract scan file folder path
	size_t foundFolder = sScanFilePath.find_last_of("\\");
	string scanFolderPath = sScanFilePath.substr( 0, foundFolder+1 );

	// Subtract scan file name
	string scanFileName  = sScanFilePath.substr( foundFolder+1 );
	size_t findExtension = scanFileName.find(".puz");
	scanFileName = scanFileName.substr( 0, findExtension );

	// Construct first scan file name
	string firstFileName = scanFolderPath.append("*.ply"); // TODO: support OBJ
	//printf( firstFileName.c_str() );
	//printf("\n");

	return firstFileName;
}

void Matcher::SortScanFileNames(vector<string> &scanFileNames, bool isPrint)
{
	// Save postfix number for each scan file
	vector<int> fileNumbers;
	for (int i=0; i<scanFileNames.size(); i++)
	{
		int fileNumber = GetScanFileNumber( scanFileNames[i] );
		fileNumbers.push_back( fileNumber );
	}

	// Print initial scan file names
	//for (int i=0; i<scanFileNames.size(); i++)
	//{
	//	printf("i=%2d ", i);
	//	printf(scanFileNames[i].c_str());
	//	printf("  %d \n", fileNumbers[i]);
	//}
	//printf("\n");

	// Sort the postfix numbers
	vector<int> sortedIndices = BubbleSort(fileNumbers, true);
	vector<string> origScanFileNames = scanFileNames;

	// Reorder scan file names according to the postfix numbers
	scanFileNames.clear();
	for (int i=0; i<sortedIndices.size(); i++)
	{
		int index = sortedIndices[i];
		scanFileNames.push_back( origScanFileNames[index] );
	}

	// Print sorted scan file names
	if ( isPrint )
	{
		printf("number of scan files: %d \n", scanFileNames.size());
		for (int i=0; i<scanFileNames.size(); i++)
		{
			printf("i=%2d ", i);
			printf(scanFileNames[i].c_str());
			printf("  %d \n", fileNumbers[i]);
		}
	}
}

int Matcher::GetScanFileNumber(string scanFileName)
{
	// Subtract scan file postfix string (with 6 characters)
	int foundName = scanFileName.find_last_of(".");
	int start  = _MAX(foundName-SCAN_FILE_POSTFIX_LENGTH, 0                       );
	int length = _MIN(foundName-start,                    SCAN_FILE_POSTFIX_LENGTH);
	string postfix = scanFileName.substr( start , length );

	// Subtract numbers only in the postfix
	string postfixNum;
	for (int i=0; i<postfix.size(); i++)
	{
		if( isdigit(postfix[i]) !=0 )
		{
			postfixNum.append( 1, postfix[i] );
		}
	}

	// Convert the postfix string into an integer
	int scanFileNum = atoi( postfixNum.c_str() );
	if ( scanFileNum < 0)
	{
		printf("Warning: The scan file postfix number (%d) should not be smaller than 0. \n", scanFileNum);
	}

	return scanFileNum;
}




//**************************************************************************************//
//                            Reconstruct 3D Model from Scans
//**************************************************************************************//

void Matcher::AlignAllScans(int startScanID, int endScanID)
{
	if ( (startScanID < 0 || startScanID >= scanList.size()) ||
		 (endScanID < 0   || endScanID >= scanList.size()) )
		return;

	myRegister.InitRegister( scanList );
	myRegister.AlignAllScans( startScanID, endScanID );

	UpdateGlobalMatrix();
}

Point Matcher::PickScanSurfPoint(int scanID, int winX, int winY)
{
	Point surfPt;

	if ( scanID >= 0 && scanID < scanList.size() )
	{
		surfPt = scanList[scanID]->PickMeshSurfPoint(winX, winY);
	}

	return surfPt;
}

void Matcher::ComputeScanDescriptor(int scanID, Point keyPoint, float radius)
{
	if ( scanID < 0 || scanID > (int)(scanList.size()-1) )
		return;

//#ifdef USE_DESCRIPTOR_OURS
//	scanList[scanID]->BuildLocalVoxelizer(keyPoint, radius, true);
//#else
//	scanList[scanID]->BuildLocalSphere(keyPoint, radius, true);
//#endif

	scanList[scanID]->BuildSGC(keyPoint, radius, false);
}




//**************************************************************************************//
//                             Poisson Reconstruction
//**************************************************************************************//

void Matcher::PoissonReconstruction(char *modelFileName, int PSRLevel)
{
	clock_t beginTime = clock();

	printf("=> Start reconstructing 3D model using PSR...\n");

	SummarizeScanPoints();

	// Write oriented points into file
	string sNPTSFileName = GetNPTSFileName( modelFileName );
	char *nptsFileName = (char *)sNPTSFileName.c_str();
	WriteModelPoints( nptsFileName ); 

	// Reconstruct the 3D model from the points
	PoissonRecon PSR;
	PSR.MPoisson(nptsFileName, modelFileName, PSRLevel);	

	// Load the 3D model from the file
	ReadModel( modelFileName );

	clock_t endTime = clock();
	float elapsed = ((float) (endTime - beginTime)) / CLOCKS_PER_SEC; // Time unit: second
	printf("\n=> PSR taken %.2f secs \n", elapsed);
}

void Matcher::SummarizeScanPoints()
{
	modelPoints.clear();

	for (int i=0; i<scanList.size(); i++)
	{
		double scanAlignRotMat[16];
		memcpy(scanAlignRotMat, scanList[i]->alignLocalMat, sizeof(double)*16);
		scanAlignRotMat[12] = 0.0;
		scanAlignRotMat[13] = 0.0;
		scanAlignRotMat[14] = 0.0;
		
		for (int j=0; j<scanList[i]->verList.size(); j++)
		{
			//printf("i=%d j=%d \n", i, j);
			Point scanVertex;

			MultiplyVector(scanList[i]->verList[j].pos, scanList[i]->alignLocalMat, scanVertex.pos);
			MultiplyVector(scanList[i]->verList[j].nor, scanAlignRotMat,            scanVertex.nor);

			modelPoints.push_back( scanVertex );
		}
	}

	printf("modelPtNum: %d \n", modelPoints.size());
}

string Matcher::GetNPTSFileName(const char modelFileName[])
{
	string sModelFilePath  = string( modelFileName );

	// Construct npts file name
	size_t findExtension = sModelFilePath.find_last_of(".");
	string nptsFileName  = sModelFilePath.substr(0, findExtension);
	nptsFileName = nptsFileName.append(".npts");

	printf( nptsFileName.c_str() );

	return nptsFileName;
}




//**************************************************************************************//
//                                Functions for Debug
//**************************************************************************************//

void Matcher::Function_Sample(int startScanID, int endScanID)
{
	for (int i=startScanID; i<=_MIN((int)scanList.size()-1, endScanID); i++)
	{
		// TODO: 1) speed up; 2) call this once only when loading scans

		clock_t beginTime = clock();
		printf("computing sample points for scan %2d ... ", i);

		scanList[i]->SampleScan();
		//scanList[i]->DescribeScan();

		clock_t endTime = clock();
		float elapsed = ((float) (endTime - beginTime)) / CLOCKS_PER_SEC; // Time unit: second
		printf("taken %.2f secs \n", elapsed);
	}
}

void Matcher::Function_Evaluate(int startScanID, int endScanID)
{
	for (int i=startScanID; i<=_MIN((int)scanList.size()-1, endScanID); i++)
	{
		clock_t beginTime = clock();
		printf("computing sample points for scan %2d ... ", i);

		scanList[i]->SampleScan();

		clock_t endTime = clock();
		float elapsed = ((float) (endTime - beginTime)) / CLOCKS_PER_SEC; // Time unit: second
		printf("taken %.2f secs \n", elapsed);
	}

	for (int i=startScanID; i<_MIN((int)scanList.size()-1, endScanID); i++)
	{
		clock_t beginTime = clock();
		printf("evaluating alignment for scan pair [%2d %2d] ... ", i, i+1);

		myRegister.EvaluateAlign(scanList[i], scanList[i+1]);

		clock_t endTime = clock();
		float elapsed = ((float) (endTime - beginTime)) / CLOCKS_PER_SEC; // Time unit: second
		printf("taken %.2f secs \n", elapsed);
	}
}

void Matcher::Function_Test(int startScanID, int endScanID)
{
	myRegister.InitRegister( scanList );

	myRegister.DescribeScans( startScanID, endScanID );

	myRegister.SetFixedScan( scanList[startScanID] );
	myRegister.AlignScan( scanList[startScanID+1] );
	//myRegister.RefineAlign_ICP( scanList[startScanID], scanList[startScanID+1] );
	myRegister.EvaluateAlign(scanList[startScanID], scanList[startScanID+1]);
}


void Matcher::Function_Match(int dataScanID, int endScanID)
{
	if ( dataScanID < 0 || dataScanID > endScanID || dataScanID > (int)(scanList.size()-1) )
		return;

	// Select a neighbor of the picked scan as the target scan
	if ( dataScanID == endScanID )
		tagtScanID = dataScanID - 1;
	else
		tagtScanID = dataScanID + 1;

	myRegister.Function_HighMatch(scanList[dataScanID], scanList[tagtScanID]);

	UpdateGlobalMatrix();
}

void Matcher::Function_Experiment(int startScanID, int endScanID)
{
	///////////////////////////////////////////////////////
	// Experiment 1: Descriptor parameters

	if ( scanList.size() != 2 )
		return;

	Experiment myExperiment;
	myExperiment.InitScanPair(scanList[0], scanList[1]);
	myExperiment.MatchScanPair();


	///////////////////////////////////////////////////////
	// Experiment 2: Descriptor accuracy
	
	//if ( scanList.size() < 2 )
	//	return;

	//myRegister.InitRegister( scanList );
	//myRegister.Experiment_Accuracy(startScanID, endScanID);
	//UpdateGlobalMatrix();


	///////////////////////////////////////////////////////
	// Experiment 3: Scan pair overlap

	//if ( scanList.size() == 0 )
	//	return;

	//myRegister.InitRegister( scanList );
	//myRegister.Experiment_Overlap(startScanID, endScanID);
	//UpdateGlobalMatrix();
}

void Matcher::Function_MeshClipping()
{
	if ( scanList.size() == 0 )
		return;

	for (int i=0; i<scanList[0]->voxelizer.localVolume.voxelGrid.size(); i++)
	{
		if( scanList[0]->voxelizer.localVolume.voxelGrid[i].state == VOXEL_CROSS_MESH )
		{
			scanList[0]->voxelizer.localVolume.voxelGrid[i].ClipLocalMesh();
			//break;
		}
	}
}

void Matcher::Function_MeshTriNum()
{
	if ( scanList.size() == 0 )
		return;

	int totalTriNum = 0;
	for (int i=0; i<scanList.size(); i++)
	{
		int scanTriNum = scanList[i]->triList.size();
		totalTriNum += scanTriNum;

		printf(" i=%d  scanTriNum: %d \n", i, scanTriNum);
	}

	int avgTriNum = totalTriNum / (float)(scanList.size());
	printf("=> TotalTriNum: %d   ScanNum: %d   AvgTriNum: %d \n", totalTriNum, scanList.size(), avgTriNum);
}

void Matcher::Function_PairwiseAlign()
{
	if ( scanList.size() != 2 )
		return;

	myRegister.InitRegister( scanList );
	myRegister.AlignScanPair();
	UpdateGlobalMatrix();
}

void Matcher::Function_RefineAlign()
{
	if ( scanList.size() != 2 )
		return;

	myRegister.InitRegister( scanList );
	myRegister.RefineAlign_ICP( scanList[0], scanList[1] );
}




//**************************************************************************************//
//                                     Draw scanList 
//**************************************************************************************//

void Matcher::WriteAlignMatrices(const char *fileName)
{
	FILE *fp;
	if ( (fp=fopen(fileName, "w+")) == NULL )
	{
		printf("Error: file not exists! \n");	
	}
	else
	{
		///////////////////////////////////////////////
		// 1. Write part info

		// Scan number
		fprintf( fp, "Scan Num: %d \n\n", scanList.size());

		// Align local matrix of each Scan 
		for (int i=0; i<scanList.size(); i++)
		{
			fprintf(fp, "Scan ID: %d \n", scanList[i]->scanID);
			for (int j=0; j<16; j++)
			{
				fprintf(fp, "%.9f ", scanList[i]->alignLocalMat[j]);
				if( (j+1)%4 == 0 )
					fprintf(fp, "\n");
			}
			fprintf(fp, "\n");	
		}

		fprintf(fp, "\n");
	}

	fclose(fp); 
}

void Matcher::ReadAlignMatrices(const char *fileName)
{
	FILE *fp;

	///////////////////////////////////////////////////////
	// 1. File does not exist

	if ( (fp=fopen(fileName, "r")) == NULL )
	{
		printf("Error: file doe not exist! \n");	
		fclose(fp);
		return;
	}


	///////////////////////////////////////////////////////
	// 2. Read align matrices file (in MATS format)

	int dummyScanNum;
	fscanf( fp, "Scan Num: %d \n\n", &dummyScanNum);

	for (int i=0; i<scanList.size(); i++)
	{
		int dummyScanID;
		fscanf(fp, "Scan ID: %d \n", &dummyScanID);

		for (int j=0; j<16; j++)
		{
			float matValue;
			fscanf(fp, "%f ", &matValue);
			scanList[i]->alignLocalMat[j] = matValue;

			if( (j+1)%4 == 0 )
				fscanf(fp, "\n");
		}
		fscanf(fp, "\n");		

	}

	//for (int i=0; i<scanList.size(); i++)
	//{
	//	PrintMatrix( scanList[i]->alignLocalMat );
	//}
}

void Matcher::WriteModelPoints(const char *fileName)
{
	FILE *fp;
	if ( (fp=fopen(fileName, "w+")) == NULL )
	{
		printf("Error: file not exists! \n");	
	}
	else
	{
		// Align local matrix of each Scan 
		for (int i=0; i<modelPoints.size(); i++)
		{
			Point pt = modelPoints[i];
			fprintf(fp, "%.9f %.9f %.9f %.9f %.9f %.9f\n", pt.pos[0], pt.pos[1], pt.pos[2], pt.nor[0], pt.nor[1], pt.nor[2]);	
		}

		fprintf(fp, "\n");
	}

	fclose(fp); 
}

void Matcher::ReadModel(const char *fileName)
{
	FILE *fp;

	///////////////////////////////////////////////////////
	// 1. File does not exist

	if ( (fp=fopen(fileName, "r")) == NULL )
	{
		printf("Error: file doe not exist! \n");	
		fclose(fp);
		return;
	}

	model.LoadMesh( fileName );
	//model.InitPose();
	model.ProcessMesh( Vector3f(0.1,0.7,0.1), 1.0 );
}




//**************************************************************************************//
//                                     Draw scanList 
//**************************************************************************************//

void Matcher::DrawScans(int startScanID, int endScanID, int mode)
{
	for (int i=startScanID; i<=_MIN((int)scanList.size()-1, endScanID); i++)
	{
		if ( mode == GL_SELECT ) 
		{
			glLoadName(i);
			glPushName(i);
		}

		scanList[i]->DrawScan( false );

		if ( mode == GL_SELECT ) 
		{
			glPopName();
		}
	}
}

void Matcher::DrawScansWire(int startScanID, int endScanID)
{
	for (int i=startScanID; i<=_MIN((int)scanList.size()-1, endScanID); i++)
	{
		scanList[i]->DrawScanWire(Vector3f(0.1,0.1,0.1), 1.0);
	}
}

void Matcher::DrawScansBBox(int startScanID, int endScanID)
{
	for (int i=startScanID; i<=_MIN((int)scanList.size()-1, endScanID); i++)
	{
		scanList[i]->DrawScanBBox();
	}
}

void Matcher::DrawScansFeaturePoints(int startScanID, int endScanID)
{
	for (int i=startScanID; i<=_MIN((int)scanList.size()-1, endScanID); i++)
	{
		scanList[i]->DrawFeaturePoints( false );
	}
}

void Matcher::DrawScansSamplePoints(int startScanID, int endScanID)
{
	for (int i=startScanID; i<=_MIN((int)scanList.size()-1, endScanID); i++)
	{
		scanList[i]->DrawSamplePoints( false );
	}
}


void Matcher::DrawAlignedScans(int startScanID, int endScanID)
{
	for (int i=startScanID; i<=_MIN((int)scanList.size()-1, endScanID); i++)
	{
		scanList[i]->DrawScan( true );
	}
}




//**************************************************************************************//
//                                    Draw Descriptor 
//**************************************************************************************//

void Matcher::DrawScansRay(int startScanID, int endScanID)
{
	for (int i=startScanID; i<=_MIN((int)scanList.size()-1, endScanID); i++)
	{
		scanList[i]->DrawRay();
	}
}

void Matcher::DrawScansKeyPoint(int startScanID, int endScanID)
{
	for (int i=startScanID; i<=_MIN((int)scanList.size()-1, endScanID); i++)
	{
		scanList[i]->DrawKeyPoint();
	}
}

void Matcher::DrawScansLRFShape(int startScanID, int endScanID)
{
	for (int i=startScanID; i<=_MIN((int)scanList.size()-1, endScanID); i++)
	{
		scanList[i]->DrawLRFShape();
	}
}

void Matcher::DrawScansLRFSphere(int startScanID, int endScanID)
{
	for (int i=startScanID; i<=_MIN((int)scanList.size()-1, endScanID); i++)
	{
		scanList[i]->DrawLRFSphere();
	}
}

void Matcher::DrawScansLRF(int startScanID, int endScanID)
{
	for (int i=startScanID; i<=_MIN((int)scanList.size()-1, endScanID); i++)
	{
		scanList[i]->DrawLRF();
	}
}

void Matcher::DrawScansLocalShape(int startScanID, int endScanID)
{
	for (int i=startScanID; i<=_MIN((int)scanList.size()-1, endScanID); i++)
	{
		scanList[i]->DrawLocalShape();
	}
}

void Matcher::DrawScansLocalShapeWire(int startScanID, int endScanID)
{
	for (int i=startScanID; i<=_MIN((int)scanList.size()-1, endScanID); i++)
	{
		scanList[i]->DrawLocalShapeWire();
	}
}

void Matcher::DrawScansLocalBBox(int startScanID, int endScanID)
{
	for (int i=startScanID; i<=_MIN((int)scanList.size()-1, endScanID); i++)
	{
		scanList[i]->DrawLocalBBox();
	}
}

void Matcher::DrawScansLocalGrid(int startScanID, int endScanID)
{
	for (int i=startScanID; i<=_MIN((int)scanList.size()-1, endScanID); i++)
	{
		scanList[i]->DrawLocalGrid();
	}
}

void Matcher::DrawTargetScanHighMatches()
{
	if ( tagtScanID < 0 || tagtScanID > (int)(scanList.size()-1) )
		return;

	myRegister.DrawHighMatches( scanList[tagtScanID] );
}

void Matcher::DrawScansLocalSphere(int startScanID, int endScanID)
{
	for (int i=startScanID; i<=_MIN((int)scanList.size()-1, endScanID); i++)
	{
		scanList[i]->DrawLocalSphere();
	}
}




//**************************************************************************************//
//                                 Draw Scan Alignment 
//**************************************************************************************//

void Matcher::DrawClosestPoint()
{
	if ( scanList.size() == NULL )
		return;

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadMatrixd( worldMatrix );
	//glLoadIdentity();

	myRegister.DrawClosestPoint();

	glPopMatrix();
}

void Matcher::ShowPickedScan(int pickScanID)
{
	if ( pickScanID < 0 || pickScanID > int(scanList.size()-1) )
		return;
		
	scanList[pickScanID]->DrawScanBBox();
}

void Matcher::DrawModelPoints()
{
	if ( scanList.size() == NULL )
		return;

	glDisable( GL_LIGHTING );
	glPointSize( 3.0 );

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadMatrixd( worldMatrix );
	//glLoadIdentity();

	glColor3f(0.1, 0.1, 0.6);
	glBegin(GL_POINTS);
	for (int i=0; i<modelPoints.size(); i=i+10)
	//for (int i=0; i<modelPoints.size(); i++)
	{
		glVertex3f(modelPoints[i].pos[0], modelPoints[i].pos[1], modelPoints[i].pos[2]);
	}
	glEnd();

	glColor3f(0.6, 0.1, 0.1);
	glBegin(GL_LINES);
	for (int i=0; i<modelPoints.size(); i=i+10)
	//for (int i=0; i<modelPoints.size(); i++)
	{
		Vector3f startPt = modelPoints[i].pos;
		Vector3f endPt   = modelPoints[i].pos + 0.05f*modelPoints[i].nor;
		glVertex3f(startPt[0], startPt[1], startPt[2]);
		glVertex3f(endPt[0],   endPt[1],   endPt[2]  );
	}
	glEnd();

	glPopMatrix();


	glPointSize( 1.0 );
	glEnable(GL_LIGHTING);
}




//**************************************************************************************//
//                                  Draw 3D Model
//**************************************************************************************//

void Matcher::DrawModel()
{
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadMatrixd( worldMatrix );

	model.DrawModel( false );

	glPopMatrix();
}

void Matcher::DrawModelWire()
{
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadMatrixd( worldMatrix );

	model.DrawModelWire(Vector3f(0.2,0.2,0.2), 1.0);

	glPopMatrix();
}

void Matcher::DrawModelBBox()
{
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadMatrixd( worldMatrix );

	model.DrawModelBBox( Vector3f(0.2,0.2,0.2), 2.0 );

	glPopMatrix();
}




//**************************************************************************************//
//                               Draw World Axes and Planes
//**************************************************************************************//

void Matcher::DrawWorldAxes(int winW, int winH, float currFovy)
{
	glDisable(GL_LIGHTING);
	glLineWidth( 2.0 );

	int viewW = 0.3 * winW;
	int viewH = 0.3 * winH;
	glViewport(0, 0, viewW, viewH);

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluPerspective(currFovy, viewW/((double)viewH), 1.0, 5000.0);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadMatrixd(worldAxesMatrix);
	glBegin(GL_LINES);
	// light red x axis arrow
	glColor3f(1.f,0.5f,.5f);
	glVertex3f(0.0f,0.0f,0.0f);
	glVertex3f(1.0f,0.0f,0.0f);
	// light green y axis arrow
	glColor3f(.5f,1.f,0.5f);
	glVertex3f(0.0f,0.0f,0.0f);
	glVertex3f(0.0f,1.0f,0.0f);
	// light blue z axis arrow
	glColor3f(.5f,.5f,1.f);
	glVertex3f(0.0f,0.0f,0.0f);
	glVertex3f(0.0f,0.0f,1.0f);
	glEnd();
	glBegin(GL_LINES);
	// x letter & arrowhead
	glColor3f(1.f,0.5f,.5f);
	glVertex3f(1.1f,0.1f,0.0f);
	glVertex3f(1.3f,-0.1f,0.0f);
	glVertex3f(1.3f,0.1f,0.0f);
	glVertex3f(1.1f,-0.1f,0.0f);
	glVertex3f(1.0f,0.0f,0.0f);
	glVertex3f(0.9f,0.1f,0.0f);
	glVertex3f(1.0f,0.0f,0.0f);
	glVertex3f(0.9f,-0.1f,0.0f);
	// y letter & arrowhead
	glColor3f(.5f,1.f,0.5f);
	glVertex3f(-0.1f,1.3f,0.0f);
	glVertex3f(0.f,1.2f,0.0f);
	glVertex3f(0.1f,1.3f,0.0f);
	glVertex3f(0.f,1.2f,0.0f);
	glVertex3f(0.f,1.2f,0.0f);
	glVertex3f(0.f,1.1f,0.0f);
	glVertex3f(0.0f,1.0f,0.0f);
	glVertex3f(0.1f,0.9f,0.0f);
	glVertex3f(0.0f,1.0f,0.0f);
	glVertex3f(-0.1f,0.9f,0.0f);
	// z letter & arrowhead
	glColor3f(.5f,.5f,1.f);
	glVertex3f(0.0f,-0.1f,1.3f);
	glVertex3f(0.0f,0.1f,1.3f);
	glVertex3f(0.0f,0.1f,1.3f);
	glVertex3f(0.0f,-0.1f,1.1f);
	glVertex3f(0.0f,-0.1f,1.1f);
	glVertex3f(0.0f,0.1f,1.1f);
	glVertex3f(0.0f,0.0f,1.0f);
	glVertex3f(0.0f,0.1f,0.9f);
	glVertex3f(0.0f,0.0f,1.0f);
	glVertex3f(0.0f,-0.1f,0.9f);
	glEnd();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix(); 

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glEnable(GL_LIGHTING);
	glLineWidth( 1.0 );
	glViewport(0, 25, winW, winH);
}

void Matcher::RotateWorldAxes(vec rotAxis, float rotAngle)
{
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glTranslatef(worldAxesMatrix[12], worldAxesMatrix[13], worldAxesMatrix[14]);
	glRotatef(rotAngle, rotAxis[0], rotAxis[1], rotAxis[2]);
	glTranslatef(-worldAxesMatrix[12], -worldAxesMatrix[13], -worldAxesMatrix[14]); 
	glMultMatrixd(worldAxesMatrix);
	glGetDoublev(GL_MODELVIEW_MATRIX, worldAxesMatrix);
	glPopMatrix();
}




//**************************************************************************************//
//                                   Transformation
//**************************************************************************************//

//**********************Manipulate Model**********************//

void Matcher::TranslateScan(int scanID, Vector3f transVec)
{
	// Calculate the inverse world matrix
	double tranWorldMat[16], inveWorldMat[16], tranInveWorldMat[16];
	memcpy(tranWorldMat, worldMatrix, sizeof(double)*16);
	transposeMat(tranWorldMat,4);
	memcpy(inveWorldMat, worldMatrix, sizeof(double)*16);
	if(invertMat(inveWorldMat, 4)==1)  printf("Inverse Matrix Error \n");
	memcpy(tranInveWorldMat, inveWorldMat, sizeof(double)*16);
	transposeMat(tranInveWorldMat, 4);

	// ObjTransVector = Inv(WorldMatrix) * WorldTransVector (only rotate part)
	double tempTransVec[3];
	vec objTransVector;
	tempTransVec[0] = transVec[0];
	tempTransVec[1] = transVec[1];
	tempTransVec[2] = transVec[2];      
	objTransVector[0] = dot3D(tempTransVec, tranInveWorldMat  );
	objTransVector[1] = dot3D(tempTransVec, tranInveWorldMat+4);
	objTransVector[2] = dot3D(tempTransVec, tranInveWorldMat+8);

	// Translate the object
	scanList[scanID]->localMatrix[12] += objTransVector[0];
	scanList[scanID]->localMatrix[13] += objTransVector[1];
	scanList[scanID]->localMatrix[14] += objTransVector[2];
}

void Matcher::RotateScan(int scanID, Vector3f rotAxis, float rotAngle)
{
	// Calculate the inverse world matrix
	double tranWorldMat[16], inveWorldMat[16], tranInveWorldMat[16];
	memcpy(tranWorldMat, worldMatrix, sizeof(double)*16);
	transposeMat(tranWorldMat,4);
	memcpy(inveWorldMat, worldMatrix, sizeof(double)*16);
	if(invertMat(inveWorldMat, 4)==1)  printf("Inverse Matrix Error \n");
	memcpy(tranInveWorldMat, inveWorldMat, sizeof(double)*16);
	transposeMat(tranInveWorldMat, 4);

	// ObjRotAxis = Inv(WorldMatrix) * WorldRotAxis	(only rotate part)
	double tempRotAxis[3];
	vec objRotAxis;
	tempRotAxis[0] = rotAxis[0];
	tempRotAxis[1] = rotAxis[1];
	tempRotAxis[2] = rotAxis[2];
	objRotAxis[0] = dot3D(tempRotAxis, tranInveWorldMat  );
	objRotAxis[1] = dot3D(tempRotAxis, tranInveWorldMat+4);
	objRotAxis[2] = dot3D(tempRotAxis, tranInveWorldMat+8);

	// Object rotation center (scan bbox center)
	//vec objRotCenter;
	//objRotCenter[0] = scanList[scanID]->localMatrix[12]; 
	//objRotCenter[1] = scanList[scanID]->localMatrix[13];
	//objRotCenter[2] = scanList[scanID]->localMatrix[14];

	// Object rotation center (picked LRF origin)
	Vector3f shiftVec;
	Vector3f keyPtPos = scanList[scanID]->voxelizer.keyPoint.pos;
	MultiplyVector(keyPtPos, scanList[scanID]->localMatrix, shiftVec);

	vec objRotCenter;
	objRotCenter[0] = shiftVec[0];
	objRotCenter[1] = shiftVec[1];
	objRotCenter[2] = shiftVec[2];


	// Rotate the object
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glTranslatef(objRotCenter[0], objRotCenter[1], objRotCenter[2]);
	glRotatef(rotAngle, objRotAxis[0], objRotAxis[1], objRotAxis[2]);
	glTranslatef(-objRotCenter[0], -objRotCenter[1], -objRotCenter[2]); 
	glMultMatrixd(scanList[scanID]->localMatrix);
	glGetDoublev(GL_MODELVIEW_MATRIX, scanList[scanID]->localMatrix);
	glPopMatrix();
}

void Matcher::ScaleWorld_AroundScan(int scanID, Vector3f scaleVec)
{
	Vector3f rotCenter;
	Vector3f keyPtPos = scanList[scanID]->voxelizer.keyPoint.pos;
	MultiplyVector(keyPtPos, scanList[scanID]->globalMatrix, rotCenter);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glTranslatef(rotCenter[0], rotCenter[1], rotCenter[2]);
	glScalef(scaleVec[0], scaleVec[1], scaleVec[2]);
	glTranslatef(-rotCenter[0], -rotCenter[1], -rotCenter[2]);
	glMultMatrixd(worldMatrix);
	glGetDoublev(GL_MODELVIEW_MATRIX, worldMatrix);
	glPopMatrix();
}


//**********************Manipulate World**********************//

void Matcher::TranslateWorld(Vector3f transVec)
{
	worldMatrix[12] += transVec[0];
	worldMatrix[13] += transVec[1];
	worldMatrix[14] += transVec[2];
}

void Matcher::RotateWorld(Vector3f rotAxis, float rotAngle)
{
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glTranslatef(worldMatrix[12], worldMatrix[13], worldMatrix[14]);
	glRotatef(rotAngle, rotAxis[0], rotAxis[1], rotAxis[2]);
	glTranslatef(-worldMatrix[12], -worldMatrix[13], -worldMatrix[14]); 
	glMultMatrixd(worldMatrix);
	glGetDoublev(GL_MODELVIEW_MATRIX, worldMatrix);
	glPopMatrix();
}

void Matcher::ScaleWorld(Vector3f scaleVec)
{
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glTranslatef(worldMatrix[12], worldMatrix[13], worldMatrix[14]);
	glScalef(scaleVec[0], scaleVec[1], scaleVec[2]);
	glTranslatef(-worldMatrix[12], -worldMatrix[13], -worldMatrix[14]); 
	glMultMatrixd(worldMatrix);
	glGetDoublev(GL_MODELVIEW_MATRIX, worldMatrix);
	glPopMatrix();
}