///////////////////////////////////////////////////////////////
//
// Helper.h
//
//   Utility Tool Functions
//
// by Song Peng ( songpeng@ustc.edu.cn )
// 
// 07/May/2015
//
///////////////////////////////////////////////////////////////

#ifndef _HELPER_H
#define _HELPER_H

#include "Vec.h"
#include <vector>
#include <GL/glut.h>


using namespace std;

///////////////////////////////////////////////////////////////
// Define Variables
///////////////////////////////////////////////////////////////

// Draw Mode
#define DRAW_POINT              0
#define DRAW_LINE               1
#define DRAW_FLAT               2
#define DRAW_SMOOTH             3

// Main Axis 
#define MAIN_AXIS_X             0
#define MAIN_AXIS_Y             1
#define MAIN_AXIS_Z             2

// Main Coordinate 
#define COORD_X                 0
#define COORD_Y                 1
#define COORD_Z                 2

#ifndef _MAX
#define _MAX(a,b)		      ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef _MIN
#define _MIN(a,b)		      ( ((a) < (b)) ? (a) : (b) )
#endif


using namespace trimesh;

///////////////////////////////////////////////////////////////
// Define Structs
///////////////////////////////////////////////////////////////


struct Triangle;
struct Circle;


///////////////////////////////////////////////////////////////
// Function Declaration
///////////////////////////////////////////////////////////////

// Random Stuff
//vector<int> GetRandomObjIndexList(int objNum, int sampleNum);
//int GetRandomObjIndex(int objNum);

//vector<int> GetRandomObjIndexList(vector<float> possibList, float alpha, int sampleNum);
//int GetRandomObjIndex(vector<float> possibList, float alpha);
//vector<float> PossibExpMapping(vector<float> possibList, float alpha);
//float GetRandomNumber(int seedIndex);

// Vector, Matrix and Array
void PrintMatrix(double matrix[]);
void PrintMatrix_3x3(double matrix[]);
vector<int> BubbleSort(vector<double> &Array, bool isAscend);  // TODO: simplify these similar functions
vector<int> BubbleSort(vector<float> &Array, bool isAscend);
vector<int> BubbleSort(vector<int> &Array, bool isAscend);

// Geometrical Calculation
void MultiplyVector(vec inVec, double inMat[16], vec &outVec);
void MultiplyMatrix(double inLefMat[16], double inRigMat[16], double outMat[16]);
bool EqualMatrix(double matA[16], double matB[16]);

// Area Calculation
//float Area2DPolygon(vector<vec> vertices, int ignoreCoord);

// Intersection Calculation
bool IsPointInBox(Vector3f point, Vector3f boxCen, Vector3f boxSize);
bool IsPlaneTriangleIntersect(vec planeNor, vec planePt, Triangle triangle);
vector<vec> CircleTriangleInterset(Circle circle, Triangle triangle);
bool CirclePlaneInterset(Circle testCircle, vec planeNor, vec planePt, vec &hitCoord0, vec &hitCoord1);
bool LineHitSphere(vec center, float radius, vec rayOrg, vec rayDir, vec &hitCoord0, vec &hitCoord1);
bool PlanePlaneInterset(vec normal0, vec point0, vec normal1, vec point1, vec &lineDir, vec &linePt);
int GetMaxAbsoluteIndex(vec lineDir);
double LineHitTriangle(vec rayOrg, vec rayDir, Triangle triangle, vec &hitCoord, bool isPrint);
bool IsPointInsideTriangle(vec point, Triangle triangle, bool isPrint);
float LinePointDistance(vec linePt0, vec linePt1, vec tagtPt);

// Intersection Drawing
void DrawTriangle(vec triV0, vec triV1, vec triV2, vec color);
void DrawPlane(vec normal, vec point, float scale, vec color);
void DrawCircle(vec center, vec normal, float radius, vec color);

// Draw Basic Element
void DrawSolidCuboid(vec minPt, vec maxPt, vec ambient, vec diffuse, vec specular);
void DrawWireCuboid(vec minPt, vec maxPt, float lineWidth, vec color);
void DrawWireSphere(vec position, float radius, vec color);
void DrawSphere(vec position, float radius, vec ambient, vec diffuse, vec specular, vec emission);
void DrawCylinder(vec p1, vec p2, float radius, vec ambient, vec diffuse, vec specular, vec emission); 
float Draw2DTextAt(char *str, float offsetX, float offsetY, float size, int alignment, int spaceShift, bool isDrawPoint);

// Draw Ground and Axes
void DrawGround();
void DrawWorldAxes(int winW, int winH, float currFovy);
void InitWorldAxesMatrix(vec scale, vec position, vec rotAxis, float rotAngle);
void RotateWorldAxes(vec rotAxis, float rotAngle);

#endif