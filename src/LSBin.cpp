///////////////////////////////////////////////////////////////
//
// LSBin.cpp
//
//   Local Spherical Bin Class
//
// by Song Peng ( songpeng@ustc.edu.cn )
// 
// 14/May/2015
//
///////////////////////////////////////////////////////////////

#include "Controls.h"
#include "Helper.h"
#include "LSBin.h"
#include "HelpDefines.h"

#define POINT_OUT_LIST   -1


//**************************************************************************************//
//                                   Initialization
//**************************************************************************************//

LSBin::LSBin()
{
	volume = 0.0;
}

LSBin::~LSBin()
{

}

void LSBin::InitBin(float _minRadius, float _maxRadius, float _minTheta, float _maxTheta, float _minAlpha, float _maxAlpha)
{
	minRadius = _minRadius;
	maxRadius = _maxRadius; 
	minTheta  = _minTheta; 
	maxTheta  = _maxTheta;
	minAlpha  = _minAlpha; 
	maxAlpha  = _maxAlpha;
}

void LSBin::PrintBin()
{
	printf("Bin: [%.3f %.3f] [%.2f %.2f] [%.2f %.2f] \n", 
		  minRadius, maxRadius, minTheta, maxTheta, minAlpha, maxAlpha);
}




//**************************************************************************************//
//                               Compute Bin Variables
//**************************************************************************************//

void LSBin::ComputeBinCenter()
{
	Vector3f corners[8];
	corners[0] = GetCartesianCoord(minRadius, minTheta, minAlpha);
	corners[1] = GetCartesianCoord(minRadius, minTheta, maxAlpha);
	corners[2] = GetCartesianCoord(minRadius, maxTheta, minAlpha);
	corners[3] = GetCartesianCoord(minRadius, maxTheta, maxAlpha);

	corners[4] = GetCartesianCoord(maxRadius, minTheta, minAlpha);
	corners[5] = GetCartesianCoord(maxRadius, minTheta, maxAlpha);
	corners[6] = GetCartesianCoord(maxRadius, maxTheta, minAlpha);
	corners[7] = GetCartesianCoord(maxRadius, maxTheta, maxAlpha);

	cen = Vector3f(0,0,0);
	for (int i=0; i<8; i++)
	{
		cen += corners[i];
	}

	cen = cen / 8.0f;
}

float LSBin::ComputeBin3DSCValue()
{
	float binValue = 0;

	for (int i=0; i<binPoints.size(); i++)
	{
		binValue += 1.0 / binPoints[i].curv;  // Note "curv" saves the number of neighboring points
	}

	//binValue = binValue / pow(volume, 1/3.0f);  // TODO: this need to be tuned 
	binValue = binValue / sqrt(volume);
	//binValue = binValue * 50;

	//printf("[%d %d %d] %.2f \n", pos[0], pos[1], pos[2], binValue);

	return binValue;
}

vector<int> LSBin::ComputeBinSHOTValues(Vector3f keyPtNormal)
{
	const int histBinNum = 10;
	vector<int> histogram;
	for (int i=0; i<histBinNum; i++)
	{
		histogram.push_back( 0 );
	}

	for (int i=0; i<binPoints.size(); i++)
	{
		float dotp = binPoints[i].nor DOT keyPtNormal;	

		//printf("dot: %.2f ", dotp);

		int index = ((dotp+1)/2.0f) * 10;
		if ( index < 0 )  index = 0;
		if ( index > 9 )  index = 9;

		histogram[index] += 1;
	}
	//printf("\n");

	//for (int i=0; i<histogram.size(); i++)
	//{
	//	printf("i=%d %d  ", i, histogram[i]);
	//}
	//printf("\n");

	return histogram;
}




//**************************************************************************************//
//                               Construct Spherical Bin
//**************************************************************************************//

void LSBin::GetBinTriangles(vector<Triangle*> localTris)
{
	for (int i=0; i<localTris.size(); i++)
	{
		Triangle *tri = localTris[i];

		if ( IsTriangleInBin( tri->v ) == true )
		{
			binTris.push_back( tri );
		}
	}
}

void LSBin::GetBinPoints(vector<Point> localPts)
{
	for (int i=0; i<localPts.size(); i++)
	{
		if ( IsPointInBin(localPts[i].pos ) == true )
		{
			binPoints.push_back( localPts[i] );
		}
	}
}

/*void LSBin::GetBinPoints(vector<Triangle*> localTris)
{
	for (int i=0; i<localTris.size(); i++)
	{
		Triangle *tri = localTris[i];

		for (int j=0; j<3; j++)
		{
			if ( IsPointInBin(tri->v[j]) == true )
			{
				Point point;
				point.pos = tri->v[j];
				point.nor = tri->normal;
				binPoints.push_back( point );
			}
		}
	}

	UpdateBinPoints( binPoints );
}

void LSBin::UpdateBinPoints(vector<Point> &binPtList)
{
	vector<Vector3f> uniquePts;
	vector<vector<Vector3f>> uniquePtsNors;
	vector<Vector3f> uniqueNors; 

	/////////////////////////////////////////////////////////////
	// 1. Save all unique points and related normals for each point
	for (int i=0; i<binPtList.size(); i++)
	{
		Point pt = binPtList[i];
		int ptIndex = GetPointIndexInList(pt.pos, uniquePts);
		if ( ptIndex == POINT_OUT_LIST )
		{
			vector<Vector3f> ptNors;
			ptNors.push_back( pt.nor );

			uniquePts.push_back( pt.pos );
			uniquePtsNors.push_back( ptNors );
		}
		else
		{
			uniquePtsNors[ptIndex].push_back( pt.nor );
		}
	}

	/////////////////////////////////////////////////////////////
	// 2. Computer average normal for each unique point
	for (int i=0; i<uniquePtsNors.size(); i++)
	{
		Vector3f avgNor = Vector3f(0,0,0);

		if ( uniquePtsNors[i].size() != 0 )
		{
			for (int j=0; j<uniquePtsNors[i].size(); j++)
			{
				avgNor += uniquePtsNors[i][j];
			}

			avgNor = avgNor / (float)uniquePtsNors[i].size();
		}

		uniqueNors.push_back( avgNor );
		//printf("i=%d  norNum: %d \n", i, uniquePtsNors[i].size());
	}

	/////////////////////////////////////////////////////////////
	// 3. Only save each unique point with an average normal
	binPtList.clear();
	for (int i=0; i<uniquePts.size(); i++)
	{
		Point pt;
		pt.pos = uniquePts[i];
		pt.nor = uniqueNors[i];

		binPtList.push_back( pt );
	}
}*/




//**************************************************************************************//
//                                  Utility Functions
//**************************************************************************************//

bool LSBin::IsTriangleInBin(Vector3f vertices[3])
{
	bool triInBin = true;

	for (int i=0; i<3; i++)
	{
		if ( IsPointInBin( vertices[i] ) == false )
		{
			triInBin = false;
		}
	}

	return triInBin;
}

bool LSBin::IsPointInBin(Vector3f point)
{
	Vector3f spheCoord = GetSphericalCoord( point );

	if ( spheCoord[0] < maxRadius && spheCoord[0] > minRadius &&
		 spheCoord[1] < maxTheta  && spheCoord[1] > minTheta  &&
		 spheCoord[2] < maxAlpha  && spheCoord[2] > minAlpha )
	{
		return true;
	}
	else
	{
		return false;
	}
}

Vector3f LSBin::GetSphericalCoord(Vector3f point)
{
	float x = point[0];
	float y = point[1];
	float z = point[2];

	Vector3f coord;

	float radius = sqrt(x*x+y*y+z*z);
	float theta  = acos(z/radius) * (180.0 / M_PI); 
	float alpha  = atan2(y, x) * (180.0 / M_PI);

	// The range of alpha is [0, 360)
	if ( alpha < 0 )
		alpha += 360.0;

	coord[0] = radius;
	coord[1] = theta;
	coord[2] = alpha;

	return coord;
}

Vector3f LSBin::GetCartesianCoord(float radius, float theta, float alpha)
{
	theta = theta * (M_PI/180.0);
	alpha = alpha * (M_PI/180.0);

	Vector3f point;
	point[0] = radius*sin(theta)*cos(alpha);
	point[1] = radius*sin(theta)*sin(alpha);
	point[2] = radius*cos(theta);

	return point;
}

int LSBin::GetPointIndexInList(Vector3f tagtPt, vector<Vector3f> ptList)
{
	for (int i=0; i<ptList.size(); i++)
	{
		if ( ptList[i] == tagtPt )
		{
			return i;
		}
	}

	return POINT_OUT_LIST;
}




//**************************************************************************************//
//                               Draw Shape In Bin 
//**************************************************************************************//

void LSBin::DrawBinTriangles(Vector3f color)
{
	glDisable(GL_LIGHTING);
	glLineWidth(5.0);
	glColor3f(color[0], color[1], color[2]);

	for (int i=0; i<binTris.size(); i++)
	{
		glBegin(GL_LINE_LOOP);
		glVertex3f( binTris[i]->v[0][0], binTris[i]->v[0][1], binTris[i]->v[0][2] );
		glVertex3f( binTris[i]->v[1][0], binTris[i]->v[1][1], binTris[i]->v[1][2] );
		glVertex3f( binTris[i]->v[2][0], binTris[i]->v[2][1], binTris[i]->v[2][2] );
		glEnd();
	}

	glEnable(GL_LIGHTING);
	glLineWidth(1.0);
}

void LSBin::DrawBinPoints(Vector3f color)
{
	glDisable(GL_LIGHTING);
	glPointSize(12.0);
	glLineWidth(3.0);

	for (int i=0; i<binPoints.size(); i++)
	{
		Point pt = binPoints[i];
		Vector3f end = pt.pos + 0.01f*pt.nor;

		//printf("i=%d curv:%.1f \n", i, binPoints[i].curv);

		color[0] = binPoints[i].curv / 30.0;
		color[1] = 0.1;
		color[2] = 0.1;

		glColor3f(color[0], color[1], color[2]);
		glBegin(GL_POINTS);
		glVertex3f(pt.pos[0], pt.pos[1], pt.pos[2]);
		glEnd();

		glColor3f(0.3,0.3,0.8);
		glBegin(GL_LINES);
		glVertex3f(pt.pos[0], pt.pos[1], pt.pos[2]);
		glVertex3f(end[0],    end[1],    end[2]);
		glEnd();
	}

	glEnable(GL_LIGHTING);
	glPointSize(1.0);
	glLineWidth(1.0);
}




//**************************************************************************************//
//                                   Draw Bin 
//**************************************************************************************//

void LSBin::DrawSphericalBin(int drawStep, Vector3f color)
{
	glColor3f(color[0], color[1], color[2]);

	//glColor3f(binValue/100.0, binValue/100.0, binValue/100.0);

	DrawBinCenter(Vector3f(0.8,0.8,0.2), 5.0);

	DrawSphericalCurve_Alpha(minRadius, minTheta, minAlpha, maxAlpha, drawStep);
	DrawSphericalCurve_Alpha(minRadius, maxTheta, minAlpha, maxAlpha, drawStep);
	DrawSphericalCurve_Theta(minRadius, minAlpha, minTheta, maxTheta, drawStep);
	DrawSphericalCurve_Theta(minRadius, maxAlpha, minTheta, maxTheta, drawStep);

	DrawSphericalCurve_Alpha(maxRadius, minTheta, minAlpha, maxAlpha, drawStep);
	DrawSphericalCurve_Alpha(maxRadius, maxTheta, minAlpha, maxAlpha, drawStep);
	DrawSphericalCurve_Theta(maxRadius, minAlpha, minTheta, maxTheta, drawStep);
	DrawSphericalCurve_Theta(maxRadius, maxAlpha, minTheta, maxTheta, drawStep);

	DrawSphericalLine_Radius(minTheta, minAlpha, minRadius, maxRadius);
	DrawSphericalLine_Radius(minTheta, maxAlpha, minRadius, maxRadius);
	DrawSphericalLine_Radius(maxTheta, minAlpha, minRadius, maxRadius);
	DrawSphericalLine_Radius(maxTheta, maxAlpha, minRadius, maxRadius);
}

void LSBin::DrawBinCenter(Vector3f color, int pointSize)
{
	glDisable(GL_LIGHTING);
	glColor3f(color[0], color[1], color[2]);
	glPointSize(pointSize);

	glBegin(GL_POINTS);
	glVertex3f(cen[0], cen[1], cen[2]);
	glEnd();

	glEnable(GL_LIGHTING);
	glPointSize(1.0);
}

void LSBin::DrawSphericalCurve_Alpha(float radius, float theta, float alphaStart, float alphaEnd, float alphaStep)
{
	for (float i=alphaStart; i<alphaEnd; i=i+alphaStep)
	{
		float alpha1 = i;
		float alpha2 = i + alphaStep;

		Vector3f point1 = GetCartesianCoord(radius, theta, alpha1);
		Vector3f point2 = GetCartesianCoord(radius, theta, alpha2);

		glBegin(GL_LINES);
		glVertex3f(point1[0], point1[1], point1[2]);
		glVertex3f(point2[0], point2[1], point2[2]);
		glEnd();
	}
}

void LSBin::DrawSphericalCurve_Theta(float radius, float alpha, float thetaStart, float thetaEnd, float thetaStep)
{
	for (float i=thetaStart; i<thetaEnd; i=i+thetaStep)
	{
		float theta1 = i;
		float theta2 = i + thetaStep;

		Vector3f point1 = GetCartesianCoord(radius, theta1, alpha);
		Vector3f point2 = GetCartesianCoord(radius, theta2, alpha);

		glBegin(GL_LINES);
		glVertex3f(point1[0], point1[1], point1[2]);
		glVertex3f(point2[0], point2[1], point2[2]);
		glEnd();
	}
}

void LSBin::DrawSphericalLine_Radius(float theta, float alpha, float radiusStart, float radiusEnd)
{
	Vector3f point1 = GetCartesianCoord(radiusStart, theta, alpha);
	Vector3f point2 = GetCartesianCoord(radiusEnd,   theta, alpha);

	glBegin(GL_LINES);
	glVertex3f(point1[0], point1[1], point1[2]);
	glVertex3f(point2[0], point2[1], point2[2]);
	glEnd();
}