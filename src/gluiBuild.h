///////////////////////////////////////////////////////////////
//
// gluiBuild.h
//
// - building the GLUI interface
//
// by Song Peng ( songpeng@ustc.edu.cn )
//
// 25/Apr/2015
//
///////////////////////////////////////////////////////////////


#ifndef _GLUI_BUILD
#define _GLUI_BUILD


#include "glui.h"

#define _MAX_ROLLOUT	             32

// ID for GLUI callback
//#define ID_IMPORT                  10
#define ID_QUIT					     11
#define ID_HELP					     12
#define ID_ABOUT				     13
#define ID_ROTATE_PANELS		     14
#define ID_CLOSE_PANELS			     15

// Show/Hide for Scans
#define ID_SHOW_SCAN                 20
#define ID_SHOW_SCAN_WIRE            21
#define ID_SHOW_SCAN_BBOX            22
#define ID_SHOW_FEATURE_POINTS       23
#define ID_SHOW_SAMPLE_POINTS        24
#define ID_SHOW_EVALUATE_POINTS      25
#define ID_SHOW_ALIGNED_SCAN         26
#define ID_SHOW_PICK_BBOX            27

// Show/Hide for Voxelizer
#define ID_SHOW_KEY_POINT            30   
#define ID_SHOW_LRF_SHAPE            31
#define ID_SHOW_LRF_SPHERE           32
#define ID_SHOW_LRF                  33
#define ID_SHOW_LOCAL_SHAPE          34
#define ID_SHOW_LOCAL_SHAPE_WIRE     35
#define ID_SHOW_LOCAL_BBOx           36
#define ID_SHOW_LOCAL_GRID           37
#define ID_SHOW_CANDI_MATCHES        38  

// Show/Hide for 3D Model
#define ID_SHOW_MODEL                40
#define ID_SHOW_MODEL_WIRE           41   
#define ID_SHOW_MODEL_BBOX           42

// Show/Hide for 3D Scene
#define ID_SHOW_AXES		         50
#define ID_SHOW_RAY                  51

// ID for application
#define ID_READ_SCANS                60
#define ID_ALIGN_SCANS               61
#define ID_RESET_SCENE               62
#define ID_FUNCTION_TEST             63
#define ID_SAVE_MATRICES             64
#define ID_READ_MATRICES             65
#define ID_SAVE_MODEL                66
#define ID_READ_MODEL                67

// ID for basic panel
#define ID_FOVY					     70
#define ID_CHANGE_BGCOLOR		     71
#define ID_DRAW_MODE                 72
#define ID_CONTROL_MODE              73

// ID for light panel
#define  ID_LIGHT_AMBIENT_RED        80
#define  ID_LIGHT_AMBIENT_GREEN      81
#define  ID_LIGHT_AMBIENT_BLUE       82
#define  ID_LIGHT_DIFFUSE_RED        83
#define  ID_LIGHT_DIFFUSE_GREEN      84
#define  ID_LIGHT_DIFFUSE_BLUE       85
#define  ID_LIGHT_SPECULAR_RED       86
#define  ID_LIGHT_SPECULAR_GREEN     87
#define  ID_LIGHT_SPECULAR_BLUE      88
#define  ID_LIGHT_ROT                90
#define  ID_RESET_LIGHT_DIR          91
#define  ID_SHOW_LIGHTDIR            92


// Build the GLUI interface
extern void initGLUI();

// rotate the rollout panels in main panel
extern void closeAllRollouts();
extern void rotatePanelsOffset(int value);
extern int  rotatePanelsTo(char key);

// Build the GLUI interface
extern void initGLUI();
static void buildInterface();

// Add Panels
static void addMatchPanel(GLUI_Panel *panel);
static void addRenderPanel(GLUI_Panel *panel);
static void addLightPanel(GLUI_Panel *panel);

// Callback function
extern void callbackGLUI(int id);
void ReadScanFiles();
void WriteAlignMatrices_MATS();
void ReadAlignMatrices_MATS();
void ReadModelFile();
void WriteModelFile();
bool GetOpenFileName(HWND hWnd, LPSTR szFile, int StringSize);
bool GetSaveFileName(HWND hWnd, LPSTR szFile, int StringSize);

// Auxiliary functions
void myExit(int exitCode);
void myHelp();
void myAbout();
extern void resetStatusBar(const char * fmt, ... );

#endif
