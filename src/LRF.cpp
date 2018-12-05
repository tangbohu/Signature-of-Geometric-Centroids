///////////////////////////////////////////////////////////////
//
// LRF.cpp
//
//   Local Reference Frame Class
//
// by Song Peng ( songpeng@ustc.edu.cn )
// 
// 27/Mar/2015
//
///////////////////////////////////////////////////////////////

#include "Controls.h"
#include "LRF.h"
#include <GL/glut.h>
#include "lineqn.h"
#include "math3D.h"
#include "Helper.h"
#include "HelpDefines.h"


//**************************************************************************************//
//                                   Initialization
//**************************************************************************************//

LRF::LRF()
{
	::identity(LRFMatrix);
}

LRF::~LRF()
{

}

LRF& LRF::operator=(const LRF &lrf)
{
	if( this == &lrf )
		return *this;

	this->keyPoint = lrf.keyPoint;
	this->radius   = lrf.radius;
	this->LRFTris  = lrf.LRFTris;

	for (int i=0; i<3; i++)
	{
		this->localAxes[i] = lrf.localAxes[i];
	}

	for (int i=0; i<16; i++)
	{
		this->LRFMatrix[i] = lrf.LRFMatrix[i];
	}

	return *this;
};

void LRF::InitLRF(Point _keyPoint, float _radius, vector<Triangle*> _LRFTris, vector<Point> _LRFpoints)
{
	keyPoint = _keyPoint;
	radius   = _radius;
	LRFTris  = _LRFTris;
	LRFPoints = _LRFpoints;
}




//**************************************************************************************//
//                            Build Local Reference Frame (LRF)
//**************************************************************************************//

void LRF::BuildLocalFrame(bool usePoint)
{
	double matrix[9];

	if (usePoint)
		ComputeScatterMatrix(LRFPoints, keyPoint.pos, radius, matrix);
	else
		ComputeScatterMatrix(LRFTris, keyPoint.pos, radius, matrix);


	if (usePoint)
	{
		LRFTris.clear();
	}

	// Compute local frame axes
	ComputeLocalFrameAxes( matrix, radius );

	// Compute LRF matrix
	ComputeLRFMatrix();
}


void LRF::ComputeLRFMatrix()
{
	::identity(LRFMatrix);

	LRFMatrix[0]  = localAxes[0][0];  LRFMatrix[1]  = localAxes[0][1];  LRFMatrix[2]  = localAxes[0][2];
	LRFMatrix[4]  = localAxes[1][0];  LRFMatrix[5]  = localAxes[1][1];  LRFMatrix[6]  = localAxes[1][2];
	LRFMatrix[8]  = localAxes[2][0];  LRFMatrix[9]  = localAxes[2][1];  LRFMatrix[10] = localAxes[2][2];
	LRFMatrix[12] = keyPoint.pos[0];  LRFMatrix[13] = keyPoint.pos[1];  LRFMatrix[14] = keyPoint.pos[2];

	//PrintMatrix(LRFMatrix);
}




//**************************************************************************************//
//                       Compute Scatter Matrix for the Triangles
//**************************************************************************************//

void LRF::ComputeScatterMatrix(vector<Point> LRFpoints, Vector3f keyPointPos, float radius, double matrix[9])
{
	// Initialize the total scatter matrix
	for (int i = 0; i<9; i++)
	{
		matrix[i] = 0;
	}

	float dist = 0;
	float all_dist = 0;
	for (int i = 0; i<LRFpoints.size(); i++)
	{

		double tempMat[9];
		vec P = LRFpoints.at(i).pos - keyPointPos;
		ComputeMatrix(P, P, tempMat);

		dist = radius - sqrt(P[0] * P[0] + P[1] * P[1] + P[2] * P[2]);
		all_dist += dist;

		for (int k = 0; k<9; k++)
		{
			matrix[k] += dist*tempMat[k];
		}		
	}
	for (int k = 0; k<9; k++)
	{
		matrix[k] /= all_dist;
	}
}



void LRF::ComputeScatterMatrix(vector<Triangle*> LRFTris, Vector3f keyPointPos, float radius, double matrix[9])
{
	// Compute local surface area for the area weight
	float totalArea = ComputeLocalShapeArea(LRFTris);

	// Initialize the total scatter matrix
	for (int i=0; i<9; i++)
	{
		matrix[i] = 0;
	}
	if (LRFTris.size() == 0)
	{
		cout << "big error" << endl;
	}

	// Compute the total scatter matrix
	for (int i=0; i<LRFTris.size(); i++)
	{
		float w1 = GetAreaWeight(*LRFTris[i], totalArea);
		float w2 = GetDistWeight(*LRFTris[i], keyPointPos, radius);

		double tempMat[9];
		GetScatterMatrix(*LRFTris[i], keyPointPos, tempMat);

		//printf("i=%d  w1: %.5f  w2: %.5f \n", i, w1, w2);

		for (int j=0; j<9; j++)
		{
			matrix[j] += w1*w2*tempMat[j];
		}
	}
}

float LRF::ComputeLocalShapeArea(vector<Triangle*> LRFTris)
{
	float totalArea = 0;

	for (int i=0; i<LRFTris.size(); i++)
	{
		totalArea += LRFTris[i]->area;
	}

	return totalArea;
}

float LRF::GetAreaWeight(Triangle tri, float totalArea)
{
	float areaWeight = tri.area / totalArea;
	return areaWeight;
}

float LRF::GetDistWeight(Triangle tri, vec keyPointPos, float radius)
{
	vec tempVec = keyPointPos - (tri.v[0]+tri.v[1]+tri.v[2])/3.0f;
	float distWeight = pow((radius - len(tempVec)), 2);
	return distWeight;
}


void LRF::GetScatterMatrix(Triangle tri, vec keyPointPos, double matrix[9])
{
	///////////////////////////////////////////////////////////
	// Compute the first part of the scatter matrix

	double matrixA[9];
	for (int i=0; i<9; i++)
	{
		matrixA[i] = 0;
	}

	for (int i=0; i<3; i++)
	{
		for (int j=0; j<3; j++)
		{
			double tempMat[9];
			vec P = tri.v[i] - keyPointPos;
			vec Q = tri.v[j] - keyPointPos;
			ComputeMatrix(P, Q, tempMat);

			//printf("i=%d j=%d : \n", i, j);
			//PrintMatrix( tempMat );

			for (int k=0; k<9; k++)
			{
				matrixA[k] += tempMat[k];
			}
		}
	}


	///////////////////////////////////////////////////////////
	// Compute the second part of the scatter matrix

	double matrixB[9];
	for (int i=0; i<9; i++)
	{
		matrixB[i] = 0;
	}

	for (int i=0; i<3; i++)
	{
		double tempMat[9];
		vec P = tri.v[i] - keyPointPos;
		vec Q = tri.v[i] - keyPointPos;
		ComputeMatrix(P, Q, tempMat);

		for (int k=0; k<9; k++)
		{
			matrixB[k] += tempMat[k];
		}
	}


	///////////////////////////////////////////////////////////
	// Compute the output matrix

	for (int i=0; i<9; i++)
	{
		matrix[i] = 1/12.0 * (matrixA[i] + matrixB[i]);
	}

	// TODO: need to check the scatter matrix later
	//for (int i=0; i<9; i++)
	//{
	//	matrix[i] = 1/12.0 * matrixB[i];
	//}

	//printf("=> Output matrix is: \n");
	//PrintMatrix( matrix );
}

void LRF::ComputeMatrix(vec P, vec Q, double mat[9])
{
	mat[0] = P[0]*Q[0];    mat[1] = P[0]*Q[1];    mat[2] = P[0]*Q[2]; 
	mat[3] = P[1]*Q[0];    mat[4] = P[1]*Q[1];    mat[5] = P[1]*Q[2]; 
	mat[6] = P[2]*Q[0];    mat[7] = P[2]*Q[1];    mat[8] = P[2]*Q[2];

	//PrintMatrix( mat );
}




//**************************************************************************************//
//                             Compute Local Frame Axes
//**************************************************************************************//

void LRF::ComputeLocalFrameAxes(double matrix[9], float radius)
{
	// Convert scatter matrix into a 3x3 matrix
	double A[3][3];
	for (int k=0; k<9; k++)
	{
		int i = k/3;
		int j = k - i*3;

		A[i][j] = matrix[k];
	}

	double d[3];
	trimesh::eigdc<double,3>(A, d);

	//printf("\n===============Local Frame Axes=================\n");

	//printf("d   [%.6f %.6f %.6f] \n", d[0], d[1], d[2]);
	//printf("A0  [%.4f %.4f %.4f] \n", A[0][0], A[1][0], A[2][0]);
	//printf("A1  [%.4f %.4f %.4f] \n", A[0][1], A[1][1], A[2][1]);
	//printf("A2  [%.4f %.4f %.4f] \n", A[0][2], A[1][2], A[2][2]);

	localAxes[0] = vec(A[0][0], A[1][0], A[2][0]);
	localAxes[1] = vec(A[0][1], A[1][1], A[2][1]);
	localAxes[2] = vec(A[0][2], A[1][2], A[2][2]);

#ifdef INPUT_HAS_NORMAL
	// Disambiguate the sign of X-axis when surface normals are available
	GetLocalAxisSign_Normal(LRFTris, LRFPoints, keyPoint, radius, localAxes[0]);
#else
	// Disambiguate the sign of X-axis using method in [Guo et al., 2013]
	GetLocalAxisSign(LRFTris, keyPoint.pos, radius, localAxes[0]);
#endif

	// TODO: 1) the sign disambiguation for Z-axis is not stable
	// Disambiguate the sign of Z-axis using method in [Guo et al., 2013]
	GetLocalAxisSign(LRFTris,LRFPoints, keyPoint.pos, radius, localAxes[2]);

	//float dotP = localAxes[0] DOT localAxes[1];
	//printf(" dotP:  %.6f \n ", dotP);

	//localAxes[2] = localAxes[0] CROSS localAxes[1];
	localAxes[1] = localAxes[2] CROSS localAxes[0];

	//printf("X-axis  [%.4f %.4f %.4f] \n", localAxes[0][0], localAxes[1][0], localAxes[2][0]);
	//printf("Y-axis  [%.4f %.4f %.4f] \n", localAxes[0][1], localAxes[1][1], localAxes[2][1]);
	//printf("Z-axis  [%.4f %.4f %.4f] \n", localAxes[0][2], localAxes[1][2], localAxes[2][2]);
}	

void LRF::GetLocalAxisSign(vector<Triangle*> LRFTris,  vector<Point> LRFPoints, vec keyPointPos, float radius, vec &localAxis)
{
	if (LRFTris.size() > 0 )
	{
		// Compute local surface area for the area weight
		float totalArea = ComputeLocalShapeArea(LRFTris);

		// Compute h to disambigulate the sign
		float h = 0;
		for (int i = 0; i < LRFTris.size(); i++)
		{
			Triangle tri = *LRFTris[i];

			float w1 = GetAreaWeight(tri, totalArea);
			float w2 = GetDistWeight(tri, keyPointPos, radius);

			float value = 0;
			for (int j = 0; j < 3; j++)
			{
				value += (tri.v[j] - keyPointPos) DOT localAxis;
			}

			h += w1*w2*value / 6.0f;
		}

		//printf(" h = %.6f \n", h);

		//if ( h < 0 )
		//	localAxis = -1.0f*localAxis;

		// For easily visualzing the x-axis
		if (h > 0)
			localAxis = -1.0f*localAxis;
	}
	else
	{
		float value = 0;

		for (int i = 0; i < LRFPoints.size(); i++)
		{
			

			if ( ((LRFPoints.at(i).pos - keyPointPos) DOT localAxis) >0)
			{
				value ++;
			}
			else
			{
				value--;
			}
		}

		if (value<0)
			localAxis = -1.0f*localAxis;
	}
}

void LRF::GetLocalAxisSign_Normal(vector<Triangle*> LRFTris, vector<Point> LRFPoints, Point keyPoint, float radius, vec &localAxis)
{
	// Get averaging local surface normal
	Vector3f avgNormal = Vector3f(0,0,0);

	if (LRFTris.size() > 0)
	{
		for (int i = 0; i < LRFTris.size(); i++)
		{
			avgNormal += LRFTris[i]->normal;
		}
		avgNormal = avgNormal / (float)LRFTris.size();

	}
	else
	{

		for (int i = 0; i < LRFPoints.size(); i++)
		{
			avgNormal += LRFPoints[i].nor;
		}

		avgNormal = keyPoint.nor;
	}


	// The x-axis should be consistent with the averaging normal
	float dotResult = avgNormal DOT localAxis;
	if ( dotResult < 0 )
	{
		localAxis = -1.0f*localAxis;
	}
}




//**************************************************************************************//
//                                     Draw LRF
//**************************************************************************************//

void LRF::DrawKeyPoint(Vector3f color)
{
	if ( keyPoint.pos == Vector3f(0,0,0) )
		return;

	glDisable( GL_LIGHTING );
	glColor3f(color[0], color[1], color[2]);
	glPointSize( 16.0 );

	// Set the offset value
	const float epsilon = 5e-4;
	float projMat[16];
	glGetFloatv(GL_PROJECTION_MATRIX, projMat);

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glTranslatef(0.0f, 0.0f, -epsilon);
	glMultMatrixf(projMat);

	// Draw the key point on the model
	glBegin(GL_POINTS);
	glVertex3f(keyPoint.pos[0], keyPoint.pos[1], keyPoint.pos[2]);
	glEnd();

	glPopMatrix();

	glMatrixMode(GL_MODELVIEW);

	glPointSize( 1.0 );
	glEnable(GL_LIGHTING);
}

void LRF::DrawLRFSphere(Vector3f color)
{
	if ( keyPoint.pos == Vector3f(0,0,0) )
		return;

	glLineWidth(1.5f);
	glDisable(GL_LIGHTING);

	// Draw the bounding sphere 
	DrawWireSphere(keyPoint.pos, radius, color);
	//DrawWireSphere(keyPoint.pos, SQUARE_ROOT_TWO*radius, color);

	glEnable(GL_LIGHTING);
	glLineWidth(1.0f);
}

void LRF::DrawLRF()
{
	if ( keyPoint.pos == Vector3f(0,0,0) )
		return;

	const float bboxMaxLen = 2.0;
	const float axisLen = 0.2 * bboxMaxLen;

	glPointSize(12.0);
	glColor3f(0.8,0.2,0.2);
	glDisable(GL_LIGHTING);

	// Draw the center of the sphere
	//glBegin(GL_POINTS);
	//glVertex3f(keyPoint.pos[0], keyPoint.pos[1], keyPoint.pos[2]);
	//glEnd();

	glLineWidth(5.0);
	glBegin(GL_LINES);
	for (int i=0; i<3; i++)
	{
		if ( i == 0 )         glColor3f(0.9f, 0.1f , 0.1f);
		else if ( i == 1 )    glColor3f(0.1f, 0.9f , 0.1f);
		else if ( i == 2 )    glColor3f(0.1f, 0.1f , 0.9f);

		glVertex3f(keyPoint.pos[0], keyPoint.pos[1], keyPoint.pos[2]);
		glVertex3f(keyPoint.pos[0]+localAxes[i][0]*axisLen, keyPoint.pos[1]+localAxes[i][1]*axisLen, keyPoint.pos[2]+localAxes[i][2]*axisLen);
	}
	glEnd();

	glEnable(GL_LIGHTING);
	glLineWidth(1.0f);
	glPointSize(1.0f);
}

void LRF::DrawLRFShape()
{
	if ( keyPoint.pos == Vector3f(0,0,0) )
		return;

	glEnable(GL_LIGHTING);

	for (int i=0; i<LRFTris.size(); i++)
	{
		Triangle *tri = LRFTris[i];

		glBegin(GL_TRIANGLES);
		glNormal3fv( tri->normal );
		glVertex3fv( tri->v[0] );
		glVertex3fv( tri->v[1] );
		glVertex3fv( tri->v[2] );
		glEnd();
	}
}