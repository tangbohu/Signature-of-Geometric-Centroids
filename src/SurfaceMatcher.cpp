///////////////////////////////////////////////////////////////
//
// SurfaceMatcher.cpp
//
//   Match Two Mesh Surfaces
//
// by Song Peng ( songpeng@ustc.edu.cn )
// 
// 11/May/2015
//
///////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#ifdef WIN32
#include <windows.h>
#include <float.h>
#else
#include <unistd.h>
#endif

#ifdef __linux
#include <values.h>
#endif

#include <string.h>

#if defined (__APPLE__) || defined(MACOSX) /* GLUT/ on Mac OS X */
#include <OPENGL/gl.h>
#include <GLUT/glut.h>                     /* for GLUT, GL, GLU */
#else                                      /* GL/ on IRIX etc. */
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif

#include "Controls.h"
#include "HelpDefines.h"
#include "gluiBuild.h"
#include "gluiKey.h"
#include "gluiMouse.h"
#include "SurfaceMatcher.h"
#include "Helper.h"
#include "Scan.h"
#include "Matcher.h"

#include "math3D.h"

///////////////////////////////////////////////////////////////
// Global Variables
///////////////////////////////////////////////////////////////


//////////////////////////
// INTERNAL

// window stuff
int  mainWindowID, winW, winH;
char programName[_MAX_STR_SIZE],winTitle[_MAX_STR_SIZE],aboutStr[_MAX_STR_SIZE];
int  backgndColor[3];

// Show/Hide for Scans
int showScan;
int showScanWire;
int showScanBBox;
int showFeatPts;
int showSampPts;
int showEvalPts;
int showAlignScan;
int showPickBBox;

// Show/Hide for Voxelizer
int showKeyPt;
int showLRFShape;
int showLRFSphere;
int showLRF;
int showLShape;
int showLShapeWire;
int showLBBox;
int showLGrid;
int showHighCandis;

// Show/Hide for 3D Model
int showModel;
int showModelWire;
int showModelBBox;

// Show/Hide for 3D Scene
int showAxes;
int showRay;

// viewing
float currFovy;
float currAspect;
int showStatusBar;

// Surface matcher
Matcher myMatcher;
int startScanID;
int endScanID;
int pickScanID;
int PSRLevel;
int seedNum;

// Various modes
int drawMode; 
int ctrlMode;

// lighting 
GLfloat lightAmbient[]  = { 0.4f, 0.4f, 0.4f, 1.0f };
GLfloat lightDiffuse[]  = { 0.8f, 0.8f, 0.8f, 1.0f };
GLfloat lightSpecular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
GLfloat lightPosition[] = { 0.0f, 0.0f, 1.0f, 0.0f };
float light_rot_mat[16] ;
int showLightDir;


//////////////////////////
// EXTERNAL

extern GLUI * glui, *gluiStatusBar   ;
extern double mouseMotionSensitivity ;



///////////////////////////////////////////////////////////////
// Function Declarations
///////////////////////////////////////////////////////////////

void myInit();
void glInit();

void resetProj();
void reshape(int w, int h);

void display(void);
void DrawLightDir();


//**************************************************************************************//
//                            Initialization and Cleanup
//**************************************************************************************//

void myInit()
{
    ////////////////////////////////////////////////////////////////////
    // Note: No OpenGL initialization should be put here, glInit instead
    ////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////
    // Setup the window

    winW = WIN_W;
    winH = WIN_H;

    // window title
    sprintf(programName,"%s %s.%s",PROGRAM_NAME,VERSION_MAJOR,VERSION_MINOR);
    sprintf(winTitle,"%s",programName);

    // about string
    strcat(aboutStr, PROGRAM_ABOUT);


    ////////////////////////////////////////////////////////////////////
    // Initial Viewing

    currFovy   = _DEFAULT_FOVY;
    currAspect = (float) winW / (float) winH;


    ////////////////////////////////////////////////////////////////////
    // Status

	showScan       = _TRUE;
	showScanWire   = _FALSE;
	showScanBBox   = _FALSE;
	showFeatPts    = _FALSE;
	showSampPts    = _FALSE;
	showEvalPts    = _FALSE;
	showAlignScan  = _TRUE;

	showPickBBox   = _TRUE;

	showKeyPt      = _TRUE;
	showLRFShape   = _FALSE;
	showLRFSphere  = _FALSE;
	showLRF        = _TRUE;
	//showLShape     = _TRUE;
	showLShape     = _FALSE;
	showLShapeWire = _FALSE;
	showLBBox      = _TRUE;
	//showLGrid      = _TRUE;
	showLGrid      = _FALSE;
	showHighCandis = _TRUE;

	showModel      = _TRUE;
	showModelWire  = _FALSE;
	showModelBBox  = _FALSE;

	showAxes       = _TRUE;
	showRay        = _FALSE;

	startScanID    = 0;
	endScanID      = 1;
	pickScanID     = -1;
	PSRLevel       = 11;

	seedNum        = 100;


    ////////////////////////////////////////////////////////////////////
    // Mouse motion sensitivity

    mouseMotionSensitivity = 1.0;


    ////////////////////////////////////////////////////////////////////
    // Status bar

    showStatusBar = _TRUE;
    resetStatusBar("Ready.");


    ////////////////////////////////////////////////////////////////////
    // background color

    backgndColor[0] = 255 ;
    backgndColor[1] = 255 ;
    backgndColor[2] = 255 ;

	////////////////////////////////////////////////////////////////////
	//  headlight initially

	light_rot_mat[0] = 1.0f ; light_rot_mat[4] = 0.0f ; light_rot_mat[8]  = 0.0f ; light_rot_mat[12] = 0.0f ; 
	light_rot_mat[1] = 0.0f ; light_rot_mat[5] = 1.0f ; light_rot_mat[9]  = 0.0f ; light_rot_mat[13] = 0.0f ; 
	light_rot_mat[2] = 0.0f ; light_rot_mat[6] = 0.0f ; light_rot_mat[10] = 1.0f ; light_rot_mat[14] = 0.0f ; 
	light_rot_mat[3] = 0.0f ; light_rot_mat[7] = 0.0f ; light_rot_mat[11] = 0.0f ; light_rot_mat[15] = 1.0f ; 

	showLightDir  = _FALSE ;

    ///////////////////////////////////////
    // no idle rotation or motion

    stopIdleMotion();

	///////////////////////////////////////
	// Surface model related

	drawMode = DRAW_SMOOTH;
	ctrlMode = MANIPU_WORLD;
}


void glInit()
{
    ////////////////////////////////////////////////////////////////////
    // 1. Color and Lighting

    glClearColor( backgndColor[0]/255.0f, backgndColor[1]/255.0f, backgndColor[2]/255.0f, 0.0f );


    ////////////////////////////////////////////////////////////////////
    // 2. various status

    // may not work on all machines: improves TEXTURE specular shading
    #ifdef GL_VERSION_1_2
    glLightModeli(GL_LIGHT_MODEL_COLOR_CONTROL, GL_SEPARATE_SPECULAR_COLOR);
    #endif
    glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER,GL_TRUE);

    glEnable(GL_DEPTH_TEST);
    glClearDepth(1.0f);
    glDepthFunc(GL_LESS);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

    // Anti-aliasing
    glHint(GL_POINT_SMOOTH_HINT,GL_NICEST);
    glEnable(GL_POINT_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);
    glEnable(GL_LINE_SMOOTH);

    // Auto-Normalization
    glEnable(GL_NORMALIZE);

    // Cull the back face (speedup and transparency)
    glCullFace( GL_BACK );


    ////////////////////////////////////////////////////////////////////
    // 3. set the projection

    resetProj();


    ////////////////////////////////////////////////////////////////////
    // 4. Lighting

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
	glLightfv(GL_LIGHT0, GL_AMBIENT,  lightAmbient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE,  lightDiffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, lightSpecular);
	glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
    glEnable(GL_LIGHT0);
    glPopMatrix();


	////////////////////////////////////////////////////////////////////
	// 5. Object material

	float MatAmbientBack[]  = {1.0f, 0.0f, 0.0f, 1.0f};	// back material
	glMaterialfv( GL_BACK,GL_AMBIENT,	  MatAmbientBack);

	myMatcher.InitWorldMatrix(INIT_WORLD_POSITION, INIT_WORLD_ROT_AXIS, INIT_WORLD_ROT_ANGLE);
	myMatcher.InitWorldAxesMatrix(INIT_WORLD_AXES_POSITION, INIT_WORLD_AXES_ROT_AXIS, INIT_WORLD_AXES_ROT_ANGLE);
}




//**************************************************************************************//
//                           Reset view / projection / Reshape
//**************************************************************************************//

void resetProj()
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(currFovy, currAspect, _Z_NEAR, _Z_FAR);
    glMatrixMode(GL_MODELVIEW);
}

void reshape(int w, int h)
{
	int tx,ty,tw,th;

	GLUI_Master.get_viewport_area( &tx, &ty, &tw, &th );
	glViewport(tx,ty,tw,th);

	winW = tw;
	winH = th;

	//    if (statusBar)
	//	statusBar->set_w( winW-10 );

	if (showStatusBar)
		glViewport(0,25,tw,th);
	else
		glViewport(0,0,tw,th+25);

	currAspect = (float) tw / (float) th;

	resetProj();
}




//**************************************************************************************//
//                                   Rendering Functions
//**************************************************************************************//

void display(void)
{
	////////////////////////////////////////////////////////////////////
	// (1) Set lighting

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadMatrixf( light_rot_mat );
	glLightfv(GL_LIGHT0, GL_AMBIENT,  lightAmbient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE,  lightDiffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, lightSpecular);
	glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
	glPopMatrix();


	////////////////////////////////////////////////////////////////////
	// (2) Set rendering mode

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	switch( drawMode )
	{
	case DRAW_POINT:
		glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);
		glShadeModel(GL_SMOOTH);
		break;
	case DRAW_LINE:
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glShadeModel(GL_SMOOTH);
		break;
	case DRAW_SMOOTH:
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glShadeModel(GL_SMOOTH);
		break;
	case DRAW_FLAT:
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glShadeModel(GL_FLAT);
		break;
	default:
		break;
	}


	////////////////////////////////////////////////////////////////////
	// (3) Draw various stuff

	myMatcher.UpdateGlobalMatrix();

	if ( showPickBBox )    myMatcher.ShowPickedScan( pickScanID );

	if ( showLightDir )     DrawLightDir();
	if ( showAxes )         myMatcher.DrawWorldAxes(winW, winH, currFovy);
	if ( showRay )          myMatcher.DrawScansRay( startScanID, endScanID );

	if ( showScan )         myMatcher.DrawScans( startScanID, endScanID, GL_RENDER );
	if ( showScanWire )     myMatcher.DrawScansWire( startScanID, endScanID );
	//if ( showScanBBox )     myMatcher.DrawScansBBox( startScanID, endScanID );
	if ( showFeatPts )      myMatcher.DrawScansFeaturePoints( startScanID, endScanID );
	if ( showSampPts )      myMatcher.DrawScansSamplePoints( startScanID, endScanID );
	if ( showAlignScan )    myMatcher.DrawAlignedScans( startScanID, endScanID );
	
	if ( showKeyPt )        myMatcher.DrawScansKeyPoint( startScanID, endScanID );
	if ( showLRFShape )     myMatcher.DrawScansLRFShape( startScanID, endScanID );
	if ( showLRFSphere )    myMatcher.DrawScansLRFSphere( startScanID, endScanID );
	if ( showLRF )          myMatcher.DrawScansLRF( startScanID, endScanID );
	if ( showLShape )       myMatcher.DrawScansLocalShape( startScanID, endScanID );
	if ( showLShapeWire )   myMatcher.DrawScansLocalShapeWire( startScanID, endScanID );
	if ( showLBBox )        myMatcher.DrawScansLocalBBox( startScanID, endScanID );
	if ( showLGrid )        myMatcher.DrawScansLocalGrid( startScanID, endScanID );
#ifdef USE_DESCRIPTOR_OURS
	if ( showHighCandis)    myMatcher.DrawTargetScanHighMatches();
#else
	if ( showHighCandis)    myMatcher.DrawScansLocalSphere(startScanID, endScanID);	
#endif

	if ( showModel )        myMatcher.DrawModel();
	if ( showModelWire )    myMatcher.DrawModelWire();
	if ( showModelBBox )    myMatcher.DrawModelBBox();

	//myMatcher.DrawClosestPoint();
	//myMatcher.DrawModelPoints();


	// Refresh
	glutSwapBuffers();
	glFlush();
}

void DrawLightDir()
{
	double mvmat[16] ;

	glGetDoublev(GL_MODELVIEW_MATRIX,mvmat);
	glPushMatrix();
	glLoadIdentity();
	glTranslatef( mvmat[12] , mvmat[13] , mvmat[14] );

	glLineWidth( 2.0 );
	glDisable( GL_LIGHTING );
	glBegin(GL_LINES);
	glColor3f(0.9f,0.9f,0.0f);
	glVertex3f( -light_rot_mat[2]  * 40.0f, -light_rot_mat[6]  * 40.0f , light_rot_mat[10] * 60.0f );
	glVertex3f(0.0f,0.0f,0.0f);
	glEnd();
	glPopMatrix();

	glLineWidth( 1.0 );
	glEnable( GL_LIGHTING );
}




//**************************************************************************************//
//                                    Main Program
//**************************************************************************************//

void usage(char *cmd)
{
    fprintf(stderr,"Usage : %s\n\n",cmd);
}

int main(int argc, char **argv)
{
    ////////////////////////////////////////////////////////////////////
    // 1. Initialization

    // initialize various non-GL stuff
    myInit();


    ////////////////////////////////////////////////////////////////////
    // 2. GLUT Initialization

    // initialization
    glutInit(&argc,argv);
    glutInitWindowSize(winW, winH);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);

    // create a new window and give a title to it
    mainWindowID = glutCreateWindow(winTitle);

    // initialize OpenGL stuff
    glInit();


    ////////////////////////////////////////////////////////////////////
    // 3. GLUI Creation

    initGLUI();


    ////////////////////////////////////////////////////////////////////
    // 4. GLUT Callbacks

    glutDisplayFunc(display);

    GLUI_Master.set_glutReshapeFunc(reshape);
    GLUI_Master.set_glutMouseFunc(mouse);
    GLUI_Master.set_glutKeyboardFunc(keyboard);
    GLUI_Master.set_glutSpecialFunc(specKey);

    GLUI_Master.set_glutIdleFunc(idle);

    glutMotionFunc(motion);

    glutMainLoop();

    return(EXIT_SUCCESS);
}


// Note:
// -----
// If we are running in Windows non-console mode, the starting point for 
// for program is WinMain instead of main().

#if (!defined(_CONSOLE)) && defined(WIN32)

int WINAPI
WinMain(HINSTANCE hInst, HINSTANCE hPrevInst, PSTR cmd, int showmode)
{
    char  *args[1024];			/* arbitrary limit, but should suffice */
    char  *p, *q, **argv = args;
    int    argc = 0;

    argv[argc++] = programName;
    p = cmd;
    for (;;) {
	if (*p == ' ')
	    while (*++p == ' ')
		;
	/* now p points at the first non-space after some spaces */
	if (*p == '\0')
	    break;    /* nothing after the spaces:  done */
	argv[argc++] = q = p;
	while (*q && *q != ' ')
	    ++q;
	/* now q points at a space or the end of the string */
	if (*q == '\0')
	    break;    /* last argv already terminated; quit */
	*q = '\0';    /* change space to terminator */
	p = q + 1;
    }
    argv[argc] = NULL;   /* terminate the argv array itself */

    return main(argc,argv);
}

#endif
