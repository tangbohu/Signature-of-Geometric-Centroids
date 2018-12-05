///////////////////////////////////////////////////////////////
//
// SurfaceMatcher.h
//
//   Match Two Mesh Surfaces
//
// by Song Peng ( songpeng@ustc.edu.cn )
// 
// 11/May/2015
// 
///////////////////////////////////////////////////////////////


#ifndef _GLUI_EXAMPLE
#define _GLUI_EXAMPLE


// Version  
#define PROGRAM_NAME			   "SurfaceMatcher"
#define VERSION_MAJOR			   "0"
#define VERSION_MINOR			   "8"

// About this Program
#define PROGRAM_ABOUT			   "Program by Song Peng\n \
- songpeng@ustc.edu.cn\n \
\n \
11 May. 2015 -"
#define PROGRAM_ABOUT_ONE_LINE	   "Program by Song Peng"	// for status bar

// Window Color
#define GLUI_BGCOLOR			   200, 200, 200
#define GLUI_FGCOLOR			    20,  20,  20
#define STATUSBAR_BGCOLOR		   200, 200, 200
#define STATUSBAR_FGCOLOR		    20,  20,  20

// window size and title
#define WIN_W				       1060
#define WIN_H				        960

// Z range
#define _Z_NEAR			            1.0f
#define _Z_FAR			         5000.0f

// default fovy
#define _DEFAULT_FOVY			   30.0f

// acceleration factor on SHIFT key
#define _SHIFT_ACCELERATION		    5.0f

// if the action between motion and mouse's up is within this limit (in seconds),
// may set idle rotation (depends on other condition as well)
#define _IDLE_ROTATION_TIME_LIMIT	0.2f

// minmax of fovy
#define FOVY_MIN			        1.0f
#define FOVY_MAX			      150.0f

// Help Menu
#ifdef _WIN32
#define _HELP_MESSAGE "  Mouse Control\n\n \
left        : rotate the world\n \
+ ALT : drag active point\n \
middle   : xy-translation\n \
+ ALT : changing fovy\n \
left+mid : xz-translation\n \
right       : scale the world\n \
\n \
+ shift : acceleration key\n \
+ ctrl  : deceleration key"
#else
#define _HELP_MESSAGE "Mouse Control\n\n \
left     : rotate the world\n \
+ ALT : rotate at viewpoint\n \
middle   : xy-translation\n \
+ ALT : changing fovy\n \
left+mid : xz-translation\n \
\n \
+ shift : acceleration key\n \
+ ctrl  : deceleration key"
#endif



#endif
