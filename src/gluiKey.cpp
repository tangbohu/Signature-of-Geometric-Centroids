///////////////////////////////////////////////////////////////
//
// gluiKey.cpp
//
// - handle keyboard event
//
// by Song Peng ( song0083@ntu.edu.sg )
//
// 29/Oct/2012
//
// 
///////////////////////////////////////////////////////////////


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef WIN32
#include <windows.h>
#endif

#ifdef WIN32
#include <windows.h>
#include <GL/gl.h>
#include <GL/glut.h>
#endif

#ifdef __APPLE__
#include <OPENGL/gl.h>
#include <GLUT/glut.h>
#endif // __APPLE__

#include "HelpDefines.h"
#include "gluiBuild.h"
#include <vector>
#include "Vec.h"
#include "Controls.h"
#include "Scan.h"
#include "Matcher.h"

using namespace std;

///////////////////////////////////////////////////////////////
// Global Variables
///////////////////////////////////////////////////////////////

extern int mainWindowID ;
extern int showAxes ;

extern Matcher myMatcher;
extern int ctrlMode;
extern int pickScanID;


///////////////////////////////////////////////////////////////
// Internal Functions
///////////////////////////////////////////////////////////////




//**************************************************************************************//
//                                  Keyboard Callbacks
//**************************************************************************************//

void keyboard ( unsigned char key, int x, int y )
{
    static int waitForSecondKey = _FALSE;
    static char str[_MAX_STR_SIZE];
    int winID,r,i,l;

    // rollout list
    extern int  nRolloutList;
    extern char rolloutChar[_MAX_ROLLOUT];

    // external function
    void myExit ( int ) ;
    void myHelp (     ) ;


    ///////////////////////////////////////////
    // Make sure pointing to the mainWindow

    winID = glutGetWindow();
    if (winID != mainWindowID)
	glutSetWindow(mainWindowID);


    ///////////////////////////////////////////
    // For second key

    if (waitForSecondKey) 
	{
		// rotate the panels
		if ('1' <= key && key <= '1'+nRolloutList-1)
			r = rotatePanelsTo(key-'1');
		else
			r = rotatePanelsTo(key);

		waitForSecondKey = _FALSE;

		// if success, just return
		if (r) 
		{
			resetStatusBar("control panel rearranged.");
			return;
		}
    }


    switch (key) 
	{


	////////////////////////////////////////////////////////////////////
	// Quit (27 -> ESC)

	case 'q' : case 'Q' : case 27 :
		myExit(EXIT_SUCCESS);
		break;

    case 'h' : case 'H' :
		myHelp();
		break;

    case 'm' :
		if (strncmp(str,"Please",6))
		{
			// init. str once only
			sprintf(str,"Rearrange subpanels by keys: %c",rolloutChar[0]);
			l = strlen(str);
			for (i=1; i<nRolloutList; i++) {
			str[l++] = ',';
			str[l++] = ' ';
			str[l++] = rolloutChar[i];
			}
			str[l++] = '.';
			str[l]   = '\0';
			strcat(str, " ( M - close all subpanels )");
		}
		resetStatusBar(str);
		waitForSecondKey = _TRUE;	// wait for the next key
		break;

    case 'M' :
		closeAllRollouts();
		resetStatusBar("all panels closed.");
		break;

    ////////////////////////////////////////////////////////////////////
    // basic control

    //case 'r' : case 'R' :
	//	callbackGLUI(ID_RESETVIEW);
	//	break;

    ////////////////////////////////////////////////////////////////////
    // status

    case 'i' : case 'I' :
		//float mat[16];
		/*
		fprintf(stderr,"Mouse button = %d\n",mouseButton);
		glGetFloatv(GL_MODELVIEW_MATRIX,mat);
		fprintf(stderr,"mvmat = [\n");
		fprintf(stderr,"  %6.3f %6.3f %6.3f %6.3f\n",  mat[0],mat[4],mat[8], mat[12]);
		fprintf(stderr,"  %6.3f %6.3f %6.3f %6.3f\n",  mat[1],mat[5],mat[9], mat[13]);
		fprintf(stderr,"  %6.3f %6.3f %6.3f %6.3f\n",  mat[2],mat[6],mat[10],mat[14]);
		fprintf(stderr,"  %6.3f %6.3f %6.3f %6.3f ]\n",mat[3],mat[7],mat[11],mat[15]);
		*/
		break;

    case 'a' : case 'A' :
		showAxes = 1 - showAxes;
		callbackGLUI( ID_SHOW_AXES );
		break;

    default:
		break;
    }


    if (winID != mainWindowID)
	glutSetWindow(winID);
}


void specKey ( int key, int x, int y )
{
}




//**************************************************************************************//
//                                Transformation Functions
//**************************************************************************************//

void Translate(vec transVec)
{
	switch(ctrlMode)
	{
	case MANIPU_SCAN:
		if ( pickScanID >= 0 )
		{
			myMatcher.TranslateScan(pickScanID, transVec);
		}
		break;

	case MANIPU_WORLD:
		myMatcher.TranslateWorld( transVec );
		break;

	default:
		break;
	}
}

void Rotate(vec rotAxis, float rotAngle)
{
	switch(ctrlMode)
	{
	case MANIPU_SCAN:
		if ( pickScanID >= 0 )
		{
			myMatcher.RotateScan(pickScanID, rotAxis, rotAngle);
		}
		break;

	case MANIPU_WORLD:
		myMatcher.RotateWorldAxes(rotAxis, rotAngle);
		myMatcher.RotateWorld(rotAxis, rotAngle);
		break;

	default:
		break;
	}
}

void Scale(vec scaleVec)
{
	switch(ctrlMode)
	{
	case MANIPU_SCAN:
		if ( pickScanID >= 0 )
		{
			myMatcher.ScaleWorld_AroundScan(pickScanID, scaleVec);
		}
		break;

	case MANIPU_WORLD:
		myMatcher.ScaleWorld(scaleVec);
		break;
	}
}
