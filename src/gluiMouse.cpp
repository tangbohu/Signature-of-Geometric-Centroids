///////////////////////////////////////////////////////////////
//
// gluiMouse.cpp
//
// - handle mouse event and idle callback
//
// by Song Peng ( song0083@ntu.edu.sg )
//
// 24/Apr/2015
//
// 
///////////////////////////////////////////////////////////////


#include <stdio.h>
#include <limits.h>

#ifdef WIN32
#include <windows.h>
#include <float.h>
#include <GL/gl.h>
#include <GL/glut.h>

#endif

#ifdef __linux
#include <values.h>
#endif

#if defined (__APPLE__) || defined(MACOSX) /* GLUT/ on Mac OS X */
#include <sys/time.h>
#include <sys/resource.h>
#include <OPENGL/gl.h>
#include <GLUT/glut.h>
#endif

#ifdef _WIN32
#include <sys/timeb.h>
#endif
#if !defined(_WIN32) && !defined (__APPLE__) && !defined(MACOSX)
#include <sys/resource.h>
#include <sys/times.h>
#endif

#include "glui.h"
#include "gluiBuild.h"
#include "SurfaceMatcher.h"
#include "Controls.h"
#include "Helper.h"
#include "HelpDefines.h"
#include "math3D.h"
#include "Vec.h"
#include <vector>
#include "Scan.h"
#include "Matcher.h"

using namespace std;
using namespace trimesh;


///////////////////////////////////////////////////////////////
// Global Variables
///////////////////////////////////////////////////////////////

//////////////////////////////////////
// Internal

// mouse status
int    mouseButton;
int    mouseModifiers;
int    preMouseX,preMouseY;
GLint  preMouse2X,preMouse2Y;			// for idle rotation
double mouseMotionSensitivity;		// sensitivity in mouse motion

// idle rotation
int    idleRotModel;				// rotate the model when idle?
int    idleRotView;				// rotate the view - modify viewing vector?
float  idleRot_nx,idleRot_ny,idleRot_nz;	// rotational axis
float  idleRot_angle;				// angle left for the idle rotation (stop when it is less than zero)
float  idleRot_speed;				// rotational speed (degree per second)
int    idleRot_aroundOrigin;			// rotate around origin?
double idleRot_time;				// time at which we previously make an idle rotation or start to do so
double timeLastCallMotion;			// time at which we previously call motion func


//////////////////////////////////////
// External

// window
extern int mainWindowID,winW,winH;

// GLUI
extern GLUI *glui;

// viewing
extern float currFovy;
extern float currAspect;

extern Matcher myMatcher;
extern int startScanID;
extern int endScanID;

extern int pickScanID;
extern int ctrlMode;


///////////////////////////////////////////////////////////////
// External Functions
///////////////////////////////////////////////////////////////

extern void Translate(vec transVec);
extern void Rotate(vec rotAxis, float rotAngle);
extern void Scale(vec scaleVec);

int PickScan(int x, int y);


//**************************************************************************************//
//                                 Internal Functions
//**************************************************************************************//

int PickScan(int x, int y) 
{
	///////////////////////////////////////////////////////
	// 1. Offset the mouseY position for object picking
	//    Note: this offset value depends on screen resolution

	const int headBarHeight = 22; 
	y = y + headBarHeight;

	if (x > winW || x < 0 || y > winH || y < 0) 
		return false;


	///////////////////////////////////////////////////////
	// 2. Get all the hit objects and hit depths

#define BUF_SIZE 100
	GLuint buffer[BUF_SIZE];
	GLint viewport[4];

	glGetIntegerv(GL_VIEWPORT, viewport);
	glSelectBuffer(BUF_SIZE, buffer);
	glRenderMode(GL_SELECT);
	glInitNames();

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluPickMatrix((GLdouble)x, (GLdouble)y, 10.0, 10.0, viewport);
	gluPerspective(currFovy, currAspect, _Z_NEAR, _Z_FAR);
	myMatcher.DrawScans(startScanID, endScanID, GL_SELECT);
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glFlush();

	glMatrixMode(GL_MODELVIEW);

	GLint hitNum = glRenderMode(GL_RENDER);
	GLuint names, *ptr;
	float zDepth[20];
	int objIndex[20];

	//printf("hitNum=%d\n", hitNum);
	ptr = (GLuint*)buffer;
	for (int i = 0; i < (int)hitNum; i++) 
	{
		names = *ptr;
		ptr++;
		zDepth[i] = *ptr/(float)0x7fffffff;
		ptr++;
		ptr++;
		objIndex[i] = *ptr;
		for (int j = 0; j < (int)names; j++) 
		{
			ptr++;
		}
	}


	///////////////////////////////////////////////////////
	// 3. Pick nothing if there is no hit object

	if ( hitNum == 0 )
	{
		int pickScanID = -1;
		ctrlMode   = MANIPU_WORLD;
		glui->sync_live();
		return pickScanID;
	}


	///////////////////////////////////////////////////////
	// 4. Pick the hit object with nearest hit depth

	vector<double> zDepthList;
	for (int i=0; i<hitNum; i++)
	{
		//printf("i=%d  depth=%.3f   scanID: %d \n", i, zDepth[i], objIndex[i]);
		zDepthList.push_back( zDepth[i] );
	}
	vector<int> sortedIndices = BubbleSort(zDepthList, true);
	int nearestIndex = sortedIndices[0];

	int pickScanID = objIndex[nearestIndex];
	//controlMode = MANIPU_OBJECT;
	//glui->sync_live();
	printf("pickScanID: %d \n", pickScanID);

	return pickScanID;
}

static void updateFovyBy(double change)
{
    double x;

    #define _FOVY_K 1.0


    ///////////////////////////////////////////
    // Compute currFovy

    // Convert : f -> x
    if (currFovy < _FOVY_K)
	x = log10(currFovy) + _FOVY_K - log(_FOVY_K);
    else
	x = currFovy;

    // adding in the x space
    x += change;

    // Convert : x -> f
    if (x > 0.0) {
	if (x > 179.9)
	    x = 179.9;
    } else {
	x = pow(10.0,x-_FOVY_K+log(_FOVY_K));
	if (x < 1e-7)
	    x = 1e-7;
    }

    currFovy = x;
}


// It returns seconds elapsed since start of program, as a double

double myTime()
{
  #ifdef _WIN32

    struct _timeb timebuffer;
    _ftime( &timebuffer );
    return (timebuffer.time + timebuffer.millitm * 0.001);

  #else

    struct rusage t;
    double procTime;

    /// (1) Get the rusage data structure at this moment
    getrusage(0,&t);

    // (2) What is the elapsed time?
    //     - CPU time = User time + system time

      // (2a) Get the seconds
      procTime = t.ru_utime.tv_sec + t.ru_stime.tv_sec;

      // (2b) More precisely! Get the microseconds too! :)
      procTime += (t.ru_utime.tv_usec + t.ru_stime.tv_usec) * 1e-6;

    // (3) Return the time!
    return procTime;

  #endif
}




//**************************************************************************************//
//                         Mouse Button and Motion Callback
//**************************************************************************************//

void mouse ( int button, int state, int x, int y )
{
    y = winH - 1 - y;

	//GLint viewport[4];
	//glGetIntegerv( GL_VIEWPORT, viewport );
	//printf(" %d %d \n", winH, viewport[3]);

    /////////////////////////////////////////////////////////
    // since users might just use left/middle button to 
    // close the pop-up menu, we cannot simply use XOR
    // here to compute mouseButton, we have to use two 
    // cases here:

    mouseModifiers = glutGetModifiers();


    // since users might just use left/middle button to 
    // close the pop-up menu, we cannot simply use XOR
    // here to compute mouseButton, we have to use two 
    // cases here:
    if (state == GLUT_DOWN)
        mouseButton = mouseButton | (1<<(button)) ;
    else
        mouseButton = mouseButton & (~(1<<(button))) ;

    mouseModifiers = glutGetModifiers();


    // for idle rotation
    if (state == GLUT_DOWN) {

	////////////////////////////////////////////
	// no more idle rotation if any

	idleRotModel = _FALSE;
	idleRotView  = _FALSE;
	preMouse2X   = -1;

    } 
	
	else
	{
		double diffTime;

		////////////////////////////////////////////
		// set idle rotation?

		diffTime = myTime() - timeLastCallMotion;

		if ( button == GLUT_LEFT_BUTTON
		  && mouseModifiers == 0
		  && mouseButton == 0
		  && diffTime < _IDLE_ROTATION_TIME_LIMIT
		  && preMouse2X != -1 )
		{
			float scale,nx,ny;

			nx = -(y - preMouse2Y) * mouseMotionSensitivity;
			ny =  (x - preMouse2X) * mouseMotionSensitivity;

			scale = sqrt(nx*nx + ny*ny);

			if (diffTime <= 1e-3)
			diffTime  = 1e-3;

			if (scale > 0.0f) 
			{
				idleRot_speed = (scale * 0.2f) / diffTime / 100.0f;
				/*
				if (mouseModifiers & GLUT_ACTIVE_SHIFT)
					idleRot_speed *= _SHIFT_ACCELERATION;
				if (mouseModifiers & GLUT_ACTIVE_CTRL)
					idleRot_speed /= _SHIFT_ACCELERATION;
				*/

				// if speed is too small (smaller than 1 degree per second),
				// no need to do idle rotation
				if (idleRot_speed >= 1.0f) 
				{

					// stop any current idle rotation
					idleRotView = _FALSE;

					idleRotModel  = _TRUE;
					idleRot_nx    = nx / scale;
					idleRot_ny    = ny / scale;
					idleRot_nz    = 0.0f;
					idleRot_time  = myTime();
					idleRot_angle = FLT_MAX;
				}
			}
		}
    }

	////////////////////////////////////////////
	// Pick scan or pick point on scan surface

	if ( button == GLUT_LEFT_BUTTON && mouseButton == 0 )
	{
		// Pick scan
		if ( mouseModifiers == GLUT_ACTIVE_ALT )
		{
			pickScanID = PickScan( x, y );
		}

		// Pick point on the selected scan surface
		else if ( mouseModifiers == GLUT_ACTIVE_CTRL )
		{
			if ( pickScanID >= 0 )
			{
				Point scanKeyPoint = myMatcher.PickScanSurfPoint(pickScanID, x, y);
				myMatcher.ComputeScanDescriptor(pickScanID, scanKeyPoint, LOCAL_SHAPE_RADIUS);
			}
		}

		glutPostRedisplay();
	}

    // record preMouse position
    preMouseX = x;
    preMouseY = y;
}


void motion ( int x, int y )
{
	GLfloat mat[16];
	GLfloat nx,ny,tx,ty,tz,scale,angle,dx,dy,accelerator;
	int myMouseModifiers,winID;

	///////////////////////////////////////////
	// Initialize mouse positions
	y = winH - 1 - y;

	if (preMouseX == x && preMouseY == y)
		return;

	dx = (x - preMouseX) * mouseMotionSensitivity;
	dy = (y - preMouseY) * mouseMotionSensitivity;

	// for idle rotation
	timeLastCallMotion = myTime();

	preMouse2X = preMouseX;
	preMouse2Y = preMouseY;

	preMouseX = x;
	preMouseY = y;

	///////////////////////////////////////////
	// Make sure pointing to the mainWindow
	winID = glutGetWindow();
	if (winID != mainWindowID)
		glutSetWindow(mainWindowID);

	///////////////////////////////////////////
	// General Mode
	accelerator = 1.0f;
	if (mouseModifiers & GLUT_ACTIVE_SHIFT)
		accelerator *= _SHIFT_ACCELERATION;
	if (mouseModifiers & GLUT_ACTIVE_CTRL)		// SHIFT + CTRL will cancel each other
		accelerator *= (1.0f / _SHIFT_ACCELERATION);
	//myMouseModifiers = mouseModifiers & (~GLUT_ACTIVE_SHIFT) & (~GLUT_ACTIVE_CTRL);

	switch (mouseButton) 
	{
		////////////////////////////////////////////////////////////
		// LEFT BUTTON (Model Rotation/Translation/Scaling)
		////////////////////////////////////////////////////////////
	case 0x1:
		switch (mouseModifiers) 
		{
			/////////////////////////////////////////
			// 1) X-Y Plane Translation (Object/World/Camera)
		case GLUT_ACTIVE_CTRL:
			{
				float transDx = 0.0015*dx;
				float transDy = 0.0015*dy;
				Translate(vec(transDx, transDy, 0.0));
				glutPostRedisplay();
			}
			break;

			/////////////////////////////////////////
			// 2) Z Axis Translation (Object/World/Camera)
		case GLUT_ACTIVE_SHIFT:
			{
				float transDz = -0.0015*dy;
				Translate(vec(0.0, 0.0, transDz));
				glutPostRedisplay();
			}
			break;

			/////////////////////////////////////////
			// 3) Rotation (Object/World/Camera)
		default :
			nx = -dy;
			ny =  dx;
			scale = sqrt(nx*nx + ny*ny);
			if (scale > 0.0f) 
			{
				//if (accelerator > 1.0001f)
				//	resetStatusBar("Arcball rotation about the model (accelerated).");
				//else if (accelerator < 0.9999f)
				//    resetStatusBar("Arcball rotation about the model (moderate).");
				//else
				//    resetStatusBar("Arcball rotation about the model.");

				nx    = nx / scale;
				ny    = ny / scale;
				angle = scale * 0.2f * accelerator;
				Rotate(vec(nx, ny, 0.0), angle);
				glutPostRedisplay();
			}
		}

		break;

		////////////////////////////////////////////////////////////
		// MIDDLE BUTTON (Control TRANSLATION AND FOVY)
		////////////////////////////////////////////////////////////
	case 0x2:
		{
			float transDx = 0.0015*dx;
			float transDy = 0.0015*dy;
			Translate(vec(transDx, transDy, 0.0));
			glutPostRedisplay();
		}
		break;


		////////////////////////////////////////////////////////////
		// LEFT+MIDDLE BUTTON (Control TRANSLATION AND FOVY)
		////////////////////////////////////////////////////////////
	case 0x3:
		break;

		////////////////////////////////////////////////////////////
		// Right BUTTON (World Rotation/Translation/Scaling)
		////////////////////////////////////////////////////////////
	case 0x4:
		switch (mouseModifiers) 
		{
			/////////////////////////////////////////
			// 1) 
		case GLUT_ACTIVE_CTRL:
			break;

			/////////////////////////////////////////
			// 2) 
		case GLUT_ACTIVE_SHIFT:
			break;

			/////////////////////////////////////////
			// 3) Scaling (Object/World)
		default :
			if( dy > 0 )
				scale = 1 + 0.001*sqrt(dx*dx + dy*dy);
			else
				scale = 1 - 0.001*sqrt(dx*dx + dy*dy);

			Scale(vec(scale, scale, scale));
			glutPostRedisplay();
			break;
		}
	}

	if (winID != mainWindowID)
		glutSetWindow(winID);
}





//**************************************************************************************//
//                                   Idle Callback
//**************************************************************************************//

void idle ( void )
{
    ////////////////////////////////////////////////////////////////////
    // If GLUI is used, need to reset window ID

    if ( glutGetWindow() != mainWindowID )
        glutSetWindow(mainWindowID);


    ////////////////////////////////////////////////////////////////////
    // Idle Rotation

    if (idleRotModel) {

	float  mat[16];
	double currTime = myTime();

	// modify the current MV matrix
	{
	    glGetFloatv(GL_MODELVIEW_MATRIX,mat);
	    glLoadIdentity();

	    glTranslatef(mat[12],mat[13],mat[14]);
	    glRotatef( idleRot_speed*(currTime-idleRot_time),
	               idleRot_nx, idleRot_ny, 0.0f );
	    glTranslatef(-mat[12],-mat[13],-mat[14]);
	    glMultMatrixf(mat);

	    glutPostRedisplay();
	}

	// update the timer
	idleRot_time = currTime;

	//fprintf(stderr,"idle rotate by axis (%f,%f,%f) and angle %f\n",
	//    idleRot_nx, idleRot_ny, 0.0f, idleRot_speed*(currTime-idleRot_time) );

    }
    else
    if (idleRotView) {

	double mvmat[16];
	double currTime = myTime();

	// modify the current MV matrix
	if (idleRot_angle > 0.0) {

	    float angle = _MIN( idleRot_angle, (float)(idleRot_speed*(currTime-idleRot_time)) );

	    glGetDoublev(GL_MODELVIEW_MATRIX,mvmat);

	    // update update view
	    glLoadIdentity();
	    if (idleRot_aroundOrigin) glTranslated(mvmat[12],mvmat[13],mvmat[14]);
	    glRotatef( angle, idleRot_nx, idleRot_ny, idleRot_nz );
	    if (idleRot_aroundOrigin) glTranslated(-mvmat[12],-mvmat[13],-mvmat[14]);

	    glMultMatrixd(mvmat);

	    glutPostRedisplay();

	    // we have rotated this amount of degree
	    idleRot_angle -= angle;

	} else
	    idleRotView = _FALSE;

	// update the timer
	idleRot_time = currTime;

    } else
	if (preMouse2X != -1)
	    preMouse2X  = -1;		// preMouse2X serves as a flag
}


void stopIdleMotion()
{
    idleRotModel = _FALSE ;
    idleRotView  = _FALSE ;
    preMouse2X   = -1     ;
}
