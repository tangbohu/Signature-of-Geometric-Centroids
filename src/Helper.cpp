///////////////////////////////////////////////////////////////
//
// Helper.cpp
//
//   Utility Tool Functions
//
// by Song Peng ( songpeng@ustc.edu.cn )
// 
// 07/May/2015
//
///////////////////////////////////////////////////////////////

#include "HelpDefines.h"
#include "HelpStructs.h"
#include "Helper.h"
#include "math3D.h"


///////////////////////////////////////////////////////////////
// Global Variables
///////////////////////////////////////////////////////////////

// World axes matrix
double worldAxesMatrix[16];    

// Render 3D Primitive
GLUquadricObj *quadObj = gluNewQuadric();


//**************************************************************************************//
//                                  Random Stuff
//**************************************************************************************//

//vector<int> GetRandomObjIndexList(int objNum, int sampleNum)
//{
//	vector<int> objIndexList;
//	if ( sampleNum > objNum || sampleNum <= 0  || objNum<=0 )
//	{
//		printf("Warning: The input objNum is not correct  objNum: %d  sampleNum: %d\n\n", objNum, sampleNum);
//		return objIndexList;
//	}
//
//	while( objIndexList.size() < sampleNum )
//	{
//		int tempIndex = GetRandomObjIndex( objNum );
//		if ( std::find(objIndexList.begin(), objIndexList.end(), tempIndex) == objIndexList.end() )
//		{
//			objIndexList.push_back(tempIndex);
//		}
//	}
//
//	return objIndexList;
//}
//
//int GetRandomObjIndex(int objNum)
//{
//	float randValue = objNum * ( rand()/(float)RAND_MAX );
//
//	return randValue;
//}


//vector<int> GetRandomObjIndexList(vector<float> possibList, float alpha, int sampleNum)
//{
//	vector<int> objIndexList;
//	if ( sampleNum > possibList.size() || sampleNum <= 0 )
//	{
//		printf("Warning: The input objNum is not correct  real: %d!  objNum: %d\n\n", possibList.size(), sampleNum);
//		return objIndexList;
//	}
//
//	while( objIndexList.size() < sampleNum )
//	{
//		int tempIndex = GetRandomObjIndex(possibList, alpha);
//		if ( std::find(objIndexList.begin(), objIndexList.end(), tempIndex) == objIndexList.end() )
//		{
//			objIndexList.push_back(tempIndex);
//		}
//
//		printf(" %d ", objIndexList.size());
//	}
//
//	return objIndexList;
//}
//
//// Input:  the selection possibility value of each object
//// Output: randomly selected object index
//int GetRandomObjIndex(vector<float> possibList, float alpha)
//{
//	if ( possibList.size() == 0)
//		return -1;
//
//	if ( possibList.size() == 1)
//		return 0;
//
//	vector<float> possibMapList = PossibExpMapping(possibList, alpha);
//
//	// Compute possibility regions for objects with possibility value [P0, P1, P2, ..., P(k-1), Pk]
//	// which is [P0, P0+P1, P0+P1+P2, ..., P0+P1+..+P(k-1), P0+P1+..+P(k-1)+Pk]
//	vector<float> possibRegions;
//	for (int i=0; i<possibMapList.size(); i++)
//	{
//		float possibSum = 0;
//		for (int j=0; j<=i; j++)
//		{
//			possibSum += possibMapList[j];
//		}
//		possibRegions.push_back(possibSum);
//	}
//
//	// Generate a random value in range [0, P0+P1+..+P(k-1)+Pk)
//	int lastObjIndex = possibMapList.size() - 1;
//	float randValue = ( rand()/(RAND_MAX+1.0) ) * possibRegions[lastObjIndex];
//
//	// Return object index by finding which region the random value falls into
//	for (int i=0; i<possibRegions.size(); i++)
//	{
//		// Find each object's possibility range
//		float regionMinValue, regionMaxValue;
//		if ( i == 0 ) 
//		{
//			regionMinValue = 0;
//			regionMaxValue = possibRegions[0];
//		}
//		else           
//		{
//			regionMinValue = possibRegions[i-1]; 
//			regionMaxValue = possibRegions[i];
//		}
//
//		// Return the randomly selected object index
//		if ( randValue > regionMinValue &&
//			randValue < regionMaxValue )
//		{
//			return i;
//		}
//	}
//}
//
//vector<float> PossibExpMapping(vector<float> possibList, float alpha)
//{
//	// Check if the input possibility values are correct
//	for (int i=0; i<possibList.size(); i++)
//	{
//		if ( possibList[i] < 0 || possibList[i] > MAX_INT )
//		{
//			printf("Warning: The input possibility values [i=%2d  P=%.2f] are not correct! \n\n", i, possibList[i]);
//		}
//	}
//
//	// Possibility value exponential mapping 
//	vector<float> possibMapList;
//	for (int i=0; i<possibList.size(); i++)
//	{
//		float tempPossib = pow(possibList[i], alpha);
//		possibMapList.push_back(tempPossib);
//		//printf("i=%d  P1 %.2f  P2 %.2f \n", i, possibList[i], possibMapList[i]);
//	}
//
//	// Calculate the sum of all the possibility values in each list
//	//float totalPossib1 = 0;
//	//float totalPossib2 = 0;
//	//for (int i=0; i<possibList.size(); i++)
//	//{
//	//	totalPossib1 += possibList[i];
//	//	totalPossib2 += possibMapList[i];
//	//}
//	////printf("Total Possib 1: %.2f \n", totalPossib1);
//	////printf("Total Possib 2: %.2f \n", totalPossib2);
//
//	//// Normalize the possibility values in each list
//	//for (int i=0; i<possibList.size(); i++)
//	//{
//	//	possibList[i]    /= totalPossib1;
//	//	possibMapList[i] /= totalPossib2;
//	//	//printf("i=%d  P1 %.2f  P2 %.2f \n", i, possibList[i], possibMapList[i]);
//	//}
//
//	return possibMapList;
//}
//
//// Return a float value in [0 1)
//float GetRandomNumber(int seedIndex)
//{
//	if ( seedIndex < 1 )
//		printf("Warning: seedIndex cannot be smaller than 1! \n\n");
//
//	float randValue = seedIndex * ( rand()/(RAND_MAX+1.0) );
//
//	if( randValue > 1.0 )
//		randValue = randValue - floor(randValue);
//	else if( randValue == 1.0 )
//		randValue = 0.999999;
//
//	return randValue;
//}




//**************************************************************************************//
//                             Vector, Matrix and Array
//**************************************************************************************//

void PrintMatrix(double matrix[])
{
	for (int i=0; i<16; i++)
	{
		printf("%.2f ", matrix[i]);
		if( (i+1)%4 == 0 )
			printf("\n");
	}
	printf("\n");	
}

void PrintMatrix_3x3(double matrix[])
{
	for (int i=0; i<9; i++)
	{
		printf("%.6f ", matrix[i]);
		if( (i+1)%3 == 0 )
			printf("\n");
	}
	printf("\n");	
}

vector<int> BubbleSort(vector<int> &Array, bool isAscend)
{
	vector<int> Indices;
	for (int i=0; i<Array.size(); i++)
		Indices.push_back(i);

	//printf("Before Sorting: ");
	//for (int i=0; i<Array.size(); i++)
	//	printf(" %d: %.2f ", Indices[i], Array[i]);
	//printf("\n");

	int i, j, flag = 1; // Set flag to 1 to start first pass
	//float tempValue;    // Holding variable
	int tempValue;    // Holding variable
	int tempIndex;      // Holding variable index 
	int num = Array.size(); 
	for(i = 1; (i <= num) && flag; i++)
	{
		flag = 0;
		for (j=0; j < (num -1); j++)
		{
			if  ( ( isAscend && Array[j+1] < Array[j]) ||
				(!isAscend && Array[j+1] > Array[j]) )
			{ 
				// Swap the values in the array
				tempValue = Array[j];   
				Array[j] = Array[j+1];
				Array[j+1] = tempValue;

				// Swap the index of these two values

				tempIndex = Indices[j];
				Indices[j] = Indices[j+1];
				Indices[j+1] = tempIndex;

				flag = 1;              
			}
		}
	}

	//printf("After Sorting:  ");
	//for (int i=0; i<Array.size(); i++)
	//	printf(" %d: %.2f ", Indices[i], Array[i]);
	//printf("\n");

	return Indices;
}

vector<int> BubbleSort(vector<float> &Array, bool isAscend)
{
	vector<int> Indices;
	for (int i=0; i<Array.size(); i++)
		Indices.push_back(i);

	//printf("Before Sorting: ");
	//for (int i=0; i<Array.size(); i++)
	//	printf(" %d: %.2f ", Indices[i], Array[i]);
	//printf("\n");

	int i, j, flag = 1; // Set flag to 1 to start first pass
	float tempValue;    // Holding variable
	//double tempValue;    // Holding variable
	int tempIndex;      // Holding variable index 
	int num = Array.size(); 
	for(i = 1; (i <= num) && flag; i++)
	{
		flag = 0;
		for (j=0; j < (num -1); j++)
		{
			if  ( ( isAscend && Array[j+1] < Array[j]) ||
				(!isAscend && Array[j+1] > Array[j]) )
			{ 
				// Swap the values in the array
				tempValue = Array[j];   
				Array[j] = Array[j+1];
				Array[j+1] = tempValue;

				// Swap the index of these two values

				tempIndex = Indices[j];
				Indices[j] = Indices[j+1];
				Indices[j+1] = tempIndex;

				flag = 1;              
			}
		}
	}

	//printf("After Sorting:  ");
	//for (int i=0; i<Array.size(); i++)
	//	printf(" %d: %.2f ", Indices[i], Array[i]);
	//printf("\n");

	return Indices;
}


vector<int> BubbleSort(vector<double> &Array, bool isAscend)
{
	vector<int> Indices;
	for (int i=0; i<Array.size(); i++)
		Indices.push_back(i);

	//printf("Before Sorting: ");
	//for (int i=0; i<Array.size(); i++)
	//	printf(" %d: %.2f ", Indices[i], Array[i]);
	//printf("\n");

	int i, j, flag = 1; // Set flag to 1 to start first pass
	//float tempValue;    // Holding variable
	double tempValue;    // Holding variable
	int tempIndex;      // Holding variable index 
	int num = Array.size(); 
	for(i = 1; (i <= num) && flag; i++)
	{
		flag = 0;
		for (j=0; j < (num -1); j++)
		{
			if  ( ( isAscend && Array[j+1] < Array[j]) ||
				(!isAscend && Array[j+1] > Array[j]) )
			{ 
				// Swap the values in the array
				tempValue = Array[j];   
				Array[j] = Array[j+1];
				Array[j+1] = tempValue;

				// Swap the index of these two values

				tempIndex = Indices[j];
				Indices[j] = Indices[j+1];
				Indices[j+1] = tempIndex;

				flag = 1;              
			}
		}
	}

	//printf("After Sorting:  ");
	//for (int i=0; i<Array.size(); i++)
	//	printf(" %d: %.2f ", Indices[i], Array[i]);
	//printf("\n");

	return Indices;
}




//**************************************************************************************//
//                              Geometrical Calculation
//**************************************************************************************//

// Note: Input matrix is an OpenGL style matrix (column major)
void MultiplyVector(vec inVec, double inMat[16], vec &outVec)
{
	double tempVec[3] = {inVec[0], inVec[1], inVec[2]};

	double inTranMat[16];
	memcpy(inTranMat, inMat, sizeof(double)*16);
	transposeMat(inTranMat, 4);

	outVec[0] = inMat[12] + dot3D(tempVec, inTranMat  );
	outVec[1] = inMat[13] + dot3D(tempVec, inTranMat+4);
	outVec[2] = inMat[14] + dot3D(tempVec, inTranMat+8);	
}

// Note: Input and output matrices are OpenGL style matrix (column major)
void MultiplyMatrix(double inLefMat[16], double inRigMat[16], double outMat[16])
{
	// Transpose left matrix
	double inLefTranMat[16];
	memcpy(inLefTranMat, inLefMat, sizeof(double)*16);
	transposeMat(inLefTranMat, 4);

	// Transpose right matrix
	double inRigTranMat[16];
	memcpy(inRigTranMat, inRigMat, sizeof(double)*16);
	transposeMat(inRigTranMat, 4);

	// Multiply two transposed matrices
	multMat(inLefTranMat, inRigTranMat, outMat, 4);
	transposeMat(outMat, 4);
}

bool EqualMatrix(double matA[16], double matB[16])
{
	bool isEqual = true;
	for (int i=0; i<16; i++)
	{
		if ( matA[i] != matB[i] )
		{
			isEqual = false;
			break;
		}
	}

	return isEqual;
}




//**************************************************************************************//
//                                 Area Calculation
//**************************************************************************************//

/*float Area2DPolygon(vector<vec> vertices, int ignoreCoord)
{
	// a degenerate polygon
	int n = vertices.size();
	if ( n < 3 ) 
		return 0; 

	//////////////////////////////////////////////
	// Get (n+1) 2D vertices for the polygon (v2D[n] = v2D[0] )

	vector<vec> vertices2D;
	if ( ignoreCoord == COORD_X )
	{
		for (int i=0; i<vertices.size(); i++)
		{
			vec tempV = vec(vertices[i].y, vertices[i].z, 0.0);
			vertices2D.push_back( tempV );
		}
		vec vn = vec(vertices[0].y, vertices[0].z, 0.0); 
		vertices2D.push_back(vn);
	}

	else if ( ignoreCoord == COORD_Y )
	{
		for (int i=0; i<vertices.size(); i++)
		{
			vec tempV = vec(vertices[i].x, vertices[i].z, 0.0);
			vertices2D.push_back( tempV );
		}
		vec vn = vec(vertices[0].x, vertices[0].z, 0.0 ); 
		vertices2D.push_back(vn);
	}

	else if ( ignoreCoord == COORD_Z )
	{
		for (int i=0; i<vertices.size(); i++)
		{
			vec tempV = vec(vertices[i].x, vertices[i].y, 0.0);
			vertices2D.push_back( tempV );
		}
		vec vn = vec(vertices[0].x, vertices[0].y, 0.0);
		vertices2D.push_back(vn);
	}


	//////////////////////////////////////////////
	// Compute area for the polygon 

	float area = 0;
	int  i, j, k;   // indices

	for (i=1, j=2, k=0; i<n; i++, j++, k++) 
	{
		area += vertices2D[i].x * (vertices2D[j].y - vertices2D[k].y);

		//printf("area = %.3f \n", area);
	}
	area += vertices2D[n].x * (vertices2D[1].y - vertices2D[n-1].y);  // wrap-around term

	return area / 2.0;
}*/




//**************************************************************************************//
//                              Intersection Calculation
//**************************************************************************************//

bool IsPointInBox(Vector3f point, Vector3f boxCen, Vector3f boxSize)
{
	Vector3f boxMatPt = boxCen + 0.5f * boxSize;
	Vector3f boxMinPt = boxCen - 0.5f * boxSize;

	if ( point[0] > boxMinPt[0] & point[0] < boxMatPt[0] &
		 point[1] > boxMinPt[1] & point[1] < boxMatPt[1] &
	     point[2] > boxMinPt[2] & point[2] < boxMatPt[2])
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool IsPlaneTriangleIntersect(vec planeNor, vec planePt, Triangle triangle)
{
	float dotR0 = planeNor DOT (triangle.v[0]-planePt);
	float dotR1 = planeNor DOT (triangle.v[1]-planePt);
	float dotR2 = planeNor DOT (triangle.v[2]-planePt);

	if ( (dotR0>0 && dotR1>0 && dotR2>0) || 
		 (dotR0<0 && dotR1<0 && dotR2<0) )
	{
		return false;
	}
	else
	{
		return true;
	}
}

//vector<vec> CircleTriangleInterset(Circle circle, Triangle triangle)
//{
//	triangle.normal = (triangle.v[2] - triangle.v[0]) CROSS (triangle.v[1] - triangle.v[0]);
//
//	vec hitPt0, hitPt1;
//	bool isIntersect = CirclePlaneInterset(circle, triangle.normal, triangle.v[0], hitPt0, hitPt1);
//
//	vector<vec> hitPointList; 
//	if ( isIntersect )
//	{
//		if ( IsPointInsideTriangle(hitPt0, triangle))
//			hitPointList.push_back( hitPt0 );
//
//		if ( IsPointInsideTriangle(hitPt1, triangle))
//			hitPointList.push_back( hitPt1 );
//	}
//
//	return hitPointList;
//}



bool CirclePlaneInterset(Circle testCircle, vec planeNor, vec planePt, vec &hitCoord0, vec &hitCoord1)
{
	vec lineDir, linePt;
	bool isIntersect = PlanePlaneInterset(testCircle.normal, testCircle.center, planeNor, planePt, lineDir, linePt);
	
	if ( !isIntersect )
		return false;

	bool isHit = LineHitSphere(testCircle.center, testCircle.radius, linePt, lineDir, hitCoord0, hitCoord1);
	return isHit;	
}

float LinePointDistance(vec linePt0, vec linePt1, vec tagtPt)
{
	vec line1 = linePt1 - linePt0;
	vec line2 = tagtPt  - linePt0;

	vec crossResult = line1 CROSS line2;
	float dist = len(crossResult) / len(line1);
	return dist;
}

double LineHitTriangle(vec rayOrg, vec rayDir, Triangle triangle, vec &hitCoord, bool isPrint)
{
	//double temp[3],s1[3],s2[3],s3[3],m,d1,d2,d3;
	//double v1[3],v2[3],v3[3],y1[3],y2[3],y3[3];

	double m,d1,d2,d3;
	vec temp;

	temp = triangle.v[0] - rayOrg;
	d1 = rayDir DOT temp;
	temp = triangle.v[1] - rayOrg;
	d2 = rayDir DOT temp;
	temp = triangle.v[2] - rayOrg;
	d3 = rayDir DOT temp;
	if (d1 <= 0.0 && d2 <= 0.0 && d3 <= 0.0)
	{
		if ( isPrint )
			printf("case 0: \n");
		return -1.0;
	}

	vec vec0 = triangle.v[2] - triangle.v[0] ;
	vec vec1 = triangle.v[1] - triangle.v[0] ;
	vec hitNormal = vec0 CROSS vec1;

	temp = triangle.v[0] - rayOrg;
	m = (hitNormal DOT temp) / (hitNormal DOT rayDir);

	// check if hitNormal need to reverse in direction 
	if ( (hitNormal DOT rayDir) >= 0) 
	{
		hitNormal[0] = -hitNormal[0];
		hitNormal[1] = -hitNormal[1];
		hitNormal[2] = -hitNormal[2];
	}
	hitCoord[0] = rayOrg[0] + m * rayDir[0];
	hitCoord[1] = rayOrg[1] + m * rayDir[1];
	hitCoord[2] = rayOrg[2] + m * rayDir[2];

	if ( IsPointInsideTriangle(hitCoord, triangle, isPrint) )
	{
		if ( isPrint )
			printf("Case 2: hit triangle. \n");
		return m;
	}
	else
	{
		if ( isPrint )
			printf("Case 1: Out triangle \n");
		return -1.0;
	}
}

bool IsPointInsideTriangle(vec point, Triangle triangle, bool isPrint)
{
	// Consider the numerical error (FLOAT_ERROR value depends on triangle scale)
	//const float FLOAT_ERROR = 0.0000001;
	//const float FLOAT_ERROR = 0.000001; 
	const float FLOAT_ERROR = 0.000002; // For Ring with voxelSize = 0.15

	vec vec0 = triangle.v[2] - triangle.v[0] ;
	vec vec1 = triangle.v[1] - triangle.v[0] ;
	vec vec2 =      point  - triangle.v[0] ;

	// Note: DotResult should include the length of vec2 (when point is very close to triangle.v[0]).
	vec normal = vec0 CROSS vec1;
	normal = normal / len(normal);
	float dotResult = normal DOT vec2; 

	if ( fabs(dotResult) > FLOAT_LARGE_ERROR )
	{
		if( isPrint )
		{
			printf("Warning: The point is not on the triangle's plane. \n");
			printf("error:  %.8f \n\n", fabs(dotResult));
		}
		return false;
	}

	float dot00 = vec0 DOT vec0 ;
	float dot01 = vec0 DOT vec1 ;
	float dot02 = vec0 DOT vec2 ;
	float dot11 = vec1 DOT vec1 ;
	float dot12 = vec1 DOT vec2 ;

	float inverDeno = 1 / (dot00 * dot11 - dot01 * dot01) ;

	float u = (dot11 * dot02 - dot01 * dot12) * inverDeno ;
	//if (u < 0 || u > 1) // if u out of range, return directly
	if ( u < 0-FLOAT_ERROR || u > 1+FLOAT_ERROR  )
	{
		if( isPrint )
			printf("Warning: u=%.12f is out of range \n", u);
		return false ;
	}

	float v = (dot00 * dot12 - dot01 * dot02) * inverDeno ;
	//if (v < 0 || v > 1) // if v out of range, return directly
	if ( v < 0-FLOAT_ERROR || v > 1+FLOAT_ERROR  )
	{
		if( isPrint )
			printf("Warning: v=%.12f is out of range \n", v);
		return false ;
	}

	if( isPrint )
		printf( "u+v = %.12f \n", u+v);

	return u + v <= 1+FLOAT_ERROR ;
}


bool LineHitSphere(vec center, float radius, vec rayOrg, vec rayDir, vec &hitCoord0, vec &hitCoord1)
{
	double A,B,C,discriminant,t,m1,m2;
	vec temp, origin=vec(0,0,0);

	origin[0] = origin[1] = origin[2] = 0;
	A = pow(dist(rayDir, origin), 2);
	temp = rayOrg - center;
	B = 2 * (rayDir DOT temp);
	C = pow( dist(temp, origin), 2) - radius*radius;

	discriminant = B*B-4*A*C;

	if (discriminant < 0)
		return false;

	t  = sqrt(discriminant);
	m1 = (-B + t) / (2*A);
	m2 = (-B - t) / (2*A);
	//printf("m1 %.3f   m2 %.3f \n", m1, m2);

	hitCoord0[0] = rayOrg[0] + m1*rayDir[0];
	hitCoord0[1] = rayOrg[1] + m1*rayDir[1];
	hitCoord0[2] = rayOrg[2] + m1*rayDir[2];

	hitCoord1[0] = rayOrg[0] + m2*rayDir[0];
	hitCoord1[1] = rayOrg[1] + m2*rayDir[1];
	hitCoord1[2] = rayOrg[2] + m2*rayDir[2];

	return true;
}


bool PlanePlaneInterset(vec normal0, vec point0, vec normal1, vec point1, vec &lineDir, vec &linePt)
{
	// If the two plane's normals are parallel, then no intersection line.
	normal0 = normal0 / len(normal0);
	normal1 = normal1 / len(normal1);
	//printf("nor0 [%.3f %.3f %.3f] \n", normal0.x, normal0.y, normal0.z);
	//printf("nor1 [%.3f %.3f %.3f] \n", normal1.x, normal1.y, normal1.z);

	if ( normal0 == normal1 )
	{
		printf("Warning: The two planes are parallel! \n\n");
		return false;
	}

	// Calculate the intersection line direction and mainAxisID
	lineDir = normal0 CROSS normal1;
	int mainAxisID = GetMaxAbsoluteIndex( lineDir );
	//printf("lineDir: [%.3f %.3f %.3f]    MainAxisID: %d\n", lineDir.x, lineDir.y, lineDir.z, mainAxisID);

	// Calculate one point on the intersection line; choose the coordinate with the largest absolute value
	// of lineDir, as this will give the most robust computations. 
	float c0 = normal0 DOT point0;
	float c1 = normal1 DOT point1;
	switch( mainAxisID )
	{
	case MAIN_AXIS_X:
		linePt[0] = 0;
		linePt[1] = (c0*normal1[2] - c1*normal0[2]) / (normal0[1]*normal1[2] - normal0[2]*normal1[1]);
		linePt[2] = (c1*normal0[1] - c0*normal1[1]) / (normal0[1]*normal1[2] - normal0[2]*normal1[1]);
		break;

	case MAIN_AXIS_Y:
		linePt[0] = (c1*normal0[2] - c0*normal1[2]) / (normal0[2]*normal1[0] - normal0[0]*normal1[2]);
		linePt[1] = 0;
		linePt[2] = (c0*normal1[0] - c1*normal0[0]) / (normal0[2]*normal1[0] - normal0[0]*normal1[2]);
		break;

	case MAIN_AXIS_Z:
		linePt[0] = (c0*normal1[1] - c1*normal0[1]) / (normal0[0]*normal1[1] - normal0[1]*normal1[0]);
		linePt[1] = (c1*normal0[0] - c0*normal1[0]) / (normal0[0]*normal1[1] - normal0[1]*normal1[0]);
		linePt[2] = 0;
		break;

	default:
		break;
	}

	//printf("linePt: [%.3f %.3f %.3f] \n", linePt.x, linePt.y, linePt.z);
	return true;
}

int GetMaxAbsoluteIndex(vec lineDir)
{
	float absX = fabs( lineDir[0] );
	float absY = fabs( lineDir[1] );
	float absZ = fabs( lineDir[2] );

	if ( absX >= absY && absX >= absZ)	  return MAIN_AXIS_X;
	if ( absY >= absX && absY >= absZ)	  return MAIN_AXIS_Y;
	if ( absZ >= absX && absZ >= absY)	  return MAIN_AXIS_Z;
}




//**************************************************************************************//
//                              Draw Plane Intersection
//**************************************************************************************//

void DrawTriangle(vec triV0, vec triV1, vec triV2, vec color)
{
	glDisable( GL_LIGHTING );
	glLineWidth( 2.0 );
	glColor3f( color[0], color[1], color[2]);

	glBegin(GL_TRIANGLES);
	glVertex3f( triV0[0], triV0[1], triV0[2] );
	glVertex3f( triV1[0], triV1[1], triV1[2] );
	glVertex3f( triV2[0], triV2[1], triV2[2] );
	glEnd();

	glLineWidth( 1.0 );
	glEnable( GL_LIGHTING );
}

void DrawPlane(vec normal, vec point, float scale, vec color)
{
	vec axisU = normal CROSS vec(1.0, 0.1, 0.1);
	vec axisV = normal CROSS axisU;

	axisU = axisU / len(axisU);
	axisV = axisV / len(axisV);

	glPointSize(15.0);
	glDisable(GL_LIGHTING);
	glColor3f(color[0], color[1], color[2]);

	glBegin(GL_QUADS);
	for (int i=0; i<360; i=i+90)
	{
		float angle = i * (2*M_PI/360.0);
		vec planePt = point + axisU*scale*cos(angle) + axisV*scale*sin(angle);
		glVertex3f(planePt[0], planePt[1], planePt[2]);
	}
	glEnd();

	glEnable(GL_LIGHTING);
	glPointSize(2.0);
}

void DrawCircle(vec center, vec normal, float radius, vec color)
{
	vec axisU = normal CROSS vec(1.0, 0.1, 0.1);
	vec axisV = normal CROSS axisU;

	axisU = axisU / len(axisU);
	axisV = axisV / len(axisV);

	glLineWidth( 3.0 );
	glColor3f(color[0], color[1], color[2]);
	glDisable(GL_LIGHTING);

	glBegin(GL_LINE_LOOP);
	for (int i=0; i<360; i++)
	{
		float angle = i * (2*M_PI/360.0);
		vec circlePt = center + axisU*radius*cos(angle) + axisV*radius*sin(angle);
		glVertex3f(circlePt[0], circlePt[1], circlePt[2]);
	}
	glEnd();

	glEnable(GL_LIGHTING);
	glLineWidth( 2.0 );
}




//**************************************************************************************//
//                                Draw Basic Element
//**************************************************************************************//

void DrawSolidCuboid(vec minPt, vec maxPt, vec ambient, vec diffuse, vec specular)
{
	// Get material properties
	const GLfloat mtlAmbient[]   = { ambient[0],  ambient[1],  ambient[2],  1.0f };
	const GLfloat mtlDiffuse[]   = { diffuse[0],  diffuse[1],  diffuse[2],  1.0f };
	const GLfloat mtlSpecular[]  = { specular[0], specular[1], specular[2], 1.0f };
	const GLfloat mtlEmission[]  = { 0.2f, 0.2f, 0.2f, 1.0f };
	const GLfloat mtlShininess[] = { 76.8f };

	glEnable(GL_LIGHTING);

	// Push back  material properties
	glPushAttrib(GL_LIGHTING_BIT);

	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,   mtlAmbient);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,   mtlDiffuse);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR,  mtlSpecular);
	glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION,  mtlEmission);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mtlShininess);

	glBegin(GL_QUADS);
	glNormal3f( 0, 0, -1);
	glVertex3f(minPt[0], minPt[1], minPt[2]); // 0
	glNormal3f( 0, 0, -1);
	glVertex3f(minPt[0], maxPt[1], minPt[2]); // 3
	glNormal3f( 0, 0, -1);
	glVertex3f(maxPt[0], maxPt[1], minPt[2]); // 2
	glNormal3f( 0, 0, -1);
	glVertex3f(maxPt[0], minPt[1], minPt[2]); // 1

	glNormal3f( 0, 0, 1);
	glVertex3f(minPt[0], minPt[1], maxPt[2]); // 4
	glNormal3f( 0, 0, 1);
	glVertex3f(maxPt[0], minPt[1], maxPt[2]); // 5
	glNormal3f( 0, 0, 1);
	glVertex3f(maxPt[0], maxPt[1], maxPt[2]); // 6
	glNormal3f( 0, 0, 1);
	glVertex3f(minPt[0], maxPt[1], maxPt[2]); // 7

	glNormal3f( 0, -1, 0);
	glVertex3f(minPt[0], minPt[1], minPt[2]); // 0
	glNormal3f( 0, -1, 0);
	glVertex3f(maxPt[0], minPt[1], minPt[2]); // 1
	glNormal3f( 0, -1, 0);
	glVertex3f(maxPt[0], minPt[1], maxPt[2]); // 5
	glNormal3f( 0, -1, 0);
	glVertex3f(minPt[0], minPt[1], maxPt[2]); // 4

	glNormal3f( 0, 1, 0);
	glVertex3f(minPt[0], maxPt[1], minPt[2]); // 3
	glNormal3f( 0, 1, 0);
	glVertex3f(minPt[0], maxPt[1], maxPt[2]); // 7
	glNormal3f( 0, 1, 0);
	glVertex3f(maxPt[0], maxPt[1], maxPt[2]); // 6
	glNormal3f( 0, 1, 0);
	glVertex3f(maxPt[0], maxPt[1], minPt[2]); // 2

	glNormal3f( -1, 0, 0);
	glVertex3f(minPt[0], minPt[1], minPt[2]); // 0
	glNormal3f( -1, 0, 0);
	glVertex3f(minPt[0], minPt[1], maxPt[2]); // 4
	glNormal3f( -1, 0, 0);
	glVertex3f(minPt[0], maxPt[1], maxPt[2]); // 7
	glNormal3f( -1, 0, 0);
	glVertex3f(minPt[0], maxPt[1], minPt[2]); // 3

	glNormal3f( 1, 0, 0);
	glVertex3f(maxPt[0], minPt[1], minPt[2]); // 1
	glNormal3f( 1, 0, 0);
	glVertex3f(maxPt[0], maxPt[1], minPt[2]); // 2
	glNormal3f( 1, 0, 0);
	glVertex3f(maxPt[0], maxPt[1], maxPt[2]); // 6
	glNormal3f( 1, 0, 0);
	glVertex3f(maxPt[0], minPt[1], maxPt[2]); // 5
	glEnd();

	glPopAttrib();
}

void DrawWireCuboid(vec minPt, vec maxPt, float lineWidth, vec color)
{
	glDisable(GL_LIGHTING);
	glColor3f(color[0], color[1], color[2]);
	glLineWidth(lineWidth);

	glBegin(GL_LINE_LOOP);
	glVertex3f(minPt[0], minPt[1], minPt[2]); // 0
	glVertex3f(maxPt[0], minPt[1], minPt[2]); // 1
	glVertex3f(maxPt[0], maxPt[1], minPt[2]); // 2
	glVertex3f(minPt[0], maxPt[1], minPt[2]); // 3
	glEnd();

	glBegin(GL_LINE_LOOP);
	glVertex3f(minPt[0], minPt[1], maxPt[2]); // 4
	glVertex3f(minPt[0], maxPt[1], maxPt[2]); // 7
	glVertex3f(maxPt[0], maxPt[1], maxPt[2]); // 6
	glVertex3f(maxPt[0], minPt[1], maxPt[2]); // 5
	glEnd();

	glBegin(GL_LINE_LOOP);
	glVertex3f(minPt[0], minPt[1], minPt[2]); // 0
	glVertex3f(minPt[0], minPt[1], maxPt[2]); // 4
	glVertex3f(maxPt[0], minPt[1], maxPt[2]); // 5
	glVertex3f(maxPt[0], minPt[1], minPt[2]); // 1
	glEnd();

	glBegin(GL_LINE_LOOP);
	glVertex3f(minPt[0], maxPt[1], minPt[2]); // 3
	glVertex3f(maxPt[0], maxPt[1], minPt[2]); // 2
	glVertex3f(maxPt[0], maxPt[1], maxPt[2]); // 6
	glVertex3f(minPt[0], maxPt[1], maxPt[2]); // 7
	glEnd();

	glEnable(GL_LIGHTING);
	glLineWidth(1.0);
}

void DrawWireSphere(vec position, float radius, vec color)
{
	glDisable(GL_LIGHTING);

	glColor3f(color[0], color[1], color[2]);

	//gluQuadricDrawStyle(quadObj, GLU_FILL);
	//gluQuadricNormals(quadObj, GLU_SMOOTH);

	glMatrixMode(GL_MODELVIEW);

	glPushMatrix();
	glTranslatef(position[0], position[1], position[2]);
	// TODO: memory leak issue of glutSolidSphere
	//gluSphere(quadObj, radius, 10, 5);
	glutWireSphere(radius, 20, 10);
	glPopMatrix();

	//glPopAttrib();

	glEnable(GL_LIGHTING);
}


void DrawSphere(vec position, float radius, vec ambient, vec diffuse, vec specular, vec emission)
{
	glEnable(GL_LIGHTING);

	glPushAttrib( GL_LIGHTING_BIT );

	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
	glMaterialfv( GL_FRONT_AND_BACK, GL_EMISSION,  emission);
	GLfloat shininess[] = { 76.8f };
	glMaterialfv( GL_FRONT_AND_BACK, GL_SHININESS, shininess);

	gluQuadricDrawStyle(quadObj, GLU_FILL);
	gluQuadricNormals(quadObj, GLU_SMOOTH);

	glPushMatrix();
	glTranslatef(position[0], position[1], position[2]);
	// TODO: memory leak issue of glutSolidSphere
	//glutSolidSphere(radius, 10, 5);
	gluSphere(quadObj, radius, 10, 5);
	//gluSphere(quadObj, radius, 100, 50);
	glPopMatrix();

	glPopAttrib();

	//glDisable(GL_LIGHTING);
}

void DrawCylinder(vec p1, vec p2, float radius, vec ambient, vec diffuse, vec specular, vec emission) 
{
	glEnable(GL_LIGHTING);

	glPushAttrib( GL_LIGHTING_BIT );

	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
	glMaterialfv( GL_FRONT_AND_BACK, GL_EMISSION,  emission);
	GLfloat shininess[] = { 76.8f };
	glMaterialfv( GL_FRONT_AND_BACK, GL_SHININESS, shininess);

	float vx = p2[0] - p1[0];
	float vy = p2[1] - p1[1];
	float vz = p2[2] - p1[2];

	//handle the degenerate case with an approximation
	if (vz == 0) vz = .00000001;

	float v = sqrt(vx * vx + vy * vy + vz * vz);
	float ax = 57.2957795 * acos(vz / v);

	if ( vz < 0.0 ) ax = -ax;
	float rx = -vy * vz;
	float ry =  vx * vz;

	gluQuadricDrawStyle(quadObj, GLU_FILL);
	gluQuadricNormals(quadObj, GLU_SMOOTH);

	glPushMatrix();
	glTranslatef(p1[0], p1[1], p1[2]);
	glRotatef(ax, rx, ry, 0.0);
	//gluCylinder(quadObj, radius, radius, v, 16, 1);
	gluCylinder(quadObj, radius, radius, v, 8, 1);
	//gluQuadricOrientation(q,GLU_INSIDE);
	//draw the first cap
	//gluDisk( q, 0.0, radius, 32, 1);
	//glTranslatef( 0,0,v );
	//draw the second cap
	//gluQuadricOrientation(q,GLU_OUTSIDE);
	//gluDisk( q, 0.0, radius, 32, 1);
	glPopMatrix();

	glPopAttrib();
}


float Draw2DTextAt(char *str, float offsetX, float offsetY, float size, int alignment, int spaceShift, bool isDrawPoint)
{
	float pixelShift;
	int   i;

	///////////////////////////////////////////
	// 1. Need the width of text?

	pixelShift = 0.0;

	// Additional space shift?
	pixelShift += glutStrokeWidth(GLUT_STROKE_ROMAN,'-') * spaceShift;

	///////////////////////////////////////////
	// 2. Reset Matrices
#define SCREEN_RATIO 1.6

	// reset the projection
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluOrtho2D(0.0,SCREEN_RATIO,0.0,1.0);

	// reset the modelview
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glTranslatef(offsetX*SCREEN_RATIO+pixelShift*size,1.0-offsetY,0.0f);
	glScalef(size,size,size);

	// TODO: this is hardcode
	if( isDrawPoint )
	{
		glPointSize(25.0);
		glBegin(GL_POINTS);
		glVertex3f(-100.0, 50.0, 0.0);
		glEnd();
		glPointSize(2.0);
	}

	///////////////////////////////////////////
	// 3. Render the text

	glDepthFunc(GL_ALWAYS);
	for (i=0; i<((int)strlen(str)); i++) {
		if (str[i] != ' ')
			glutStrokeCharacter(GLUT_STROKE_ROMAN,str[i]);

		else
			glTranslatef(200000*size,0.0,0.0);
	}
	glDepthFunc(GL_LESS);

	///////////////////////////////////////////
	// 4. Reset the matrices

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();

	return pixelShift ;
}




//**************************************************************************************//
//                                Draw Ground and Axes
//**************************************************************************************//

void DrawGround()
{
	glDisable(GL_LIGHTING);
	glColor3f(0.4f, 0.4f, 0.4f);
	glLineWidth( 1.0 );

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	//glLoadIdentity();
	//glTranslatef(0.0, -0.2, -3.2);  
	//glRotatef(15.0, 1.0, 0.0, 0.0);
	glScalef(1.5, 1.5, 1.5);

	for (int i=0; i<=20; i++)
	{
		glBegin(GL_LINES);
		glVertex3f( -2.0+0.2*i,0.0,  -2.0 );
		glVertex3f( -2.0+0.2*i, 0.0, 2.0 );
		glEnd();

		glBegin(GL_LINES);
		glVertex3f( -2.0, 0.0, -2.0+0.2*i );
		glVertex3f( 2.0, 0.0, -2.0+0.2*i );
		glEnd();
	}
	glPopMatrix();

	glEnable(GL_LIGHTING);
	glLineWidth( 1.0 );
}

void DrawWorldAxes(int winW, int winH, float currFovy)
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

void InitWorldAxesMatrix(vec scale, vec position, vec rotAxis, float rotAngle)
{
	glPushMatrix();   
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(position[0], position[1], position[2]);  
	glRotatef(rotAngle, rotAxis[0], rotAxis[1], rotAxis[2]);
	glScalef(scale[0], scale[1], scale[2]);
	glGetDoublev(GL_MODELVIEW_MATRIX, worldAxesMatrix);
	glPopMatrix();
}

void RotateWorldAxes(vec rotAxis, float rotAngle)
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
