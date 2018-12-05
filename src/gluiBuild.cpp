///////////////////////////////////////////////////////////////
//
// gluiBuild.cpp
//
// - building the GLUI interface
//
// by Song Peng ( songpeng@ustc.edu.cn )
//
// 11/May/2015
//
///////////////////////////////////////////////////////////////


#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#ifdef WIN32
#include <windows.h>
#include <float.h>
#else
#include <limits.h>
#endif


#include "gluiBuild.h"
#include "SurfaceMatcher.h"
#include "HelpDefines.h"
#include "Scan.h"
#include "Matcher.h"


////////////////////////////////////////////////////////
// (1) Global Variables


#define _UNKNOWN	"---"


//////////////////////////////////////////////
// (A) External from driver program

// window
extern int winW,winH;
extern int mainWindowID;
extern int backgndColor[3];
extern char programName[_MAX_STR_SIZE];
extern char winTitle[_MAX_STR_SIZE];
extern char aboutStr[_MAX_STR_SIZE];
extern int  preMouseX;
extern int  preMouseY;

// external state variables
extern float currFovy;
extern int   showStatusBar;

// Show/Hide for Scans
extern int   showScan;
extern int   showScanWire;
extern int   showScanBBox;
extern int   showFeatPts;
extern int   showSampPts;
extern int   showEvalPts;
extern int   showAlignScan;

extern int   showPickBBox;

// Show/Hide for Voxelizer
extern int   showKeyPt;
extern int   showLRFShape;
extern int   showLRFSphere;
extern int   showLRF;
extern int   showLShape;
extern int   showLShapeWire;
extern int   showLBBox;
extern int   showLGrid;
extern int   showHighCandis;

// Show/Hide for 3D Model
extern int   showModel;
extern int   showModelWire;
extern int   showModelBBox;

// Show/Hide for 3D Scene
extern int   showAxes;
extern int   showRay;

// Mesh model
extern Matcher myMatcher;
extern int   startScanID;
extern int   endScanID;
extern int   pickScanID;
extern int   PSRLevel;
extern int   seedNum;

// Various mode
extern int   drawMode;
extern int   ctrlMode;

// lighting
extern GLfloat lightAmbient[4];
extern GLfloat lightDiffuse[4];
extern GLfloat lightSpecular[4];
extern GLfloat lightPosition[4];
extern float light_rot_mat[16] ;
extern int   showLightDir;


//////////////////////////////////////////////
// (B) Internal

vec scaleVec = vec(1,1,1);
float rotAngle = 0.0;
vec rotAxis = vec(1,0,0);

// using GLUI
GLUI *glui          = NULL;
GLUI *gluiStatusBar = NULL;

// list of rollout
GLUI_Rollout *rolloutList[_MAX_ROLLOUT];
int  rolloutStart;
int  nRolloutList  = 0;
int  rolloutOffset = 0;
char rolloutChar[_MAX_ROLLOUT];

// Status Bar
GLUI_StaticText *statusBar;
char             statusStr[_MAX_STR_SIZE];


//**************************************************************************************//
//                            Rotating the RollOuts in Main Panel
//**************************************************************************************//

void closeAllRollouts()
{
    int i;

    for (i=0; i<nRolloutList; i++)
	rolloutList[i]->close();
}


static void resetPanels()
{
    int offset,i;

    // remove all the rollouts
    while ( glui->remove_control(glui->get_main_panel(),rolloutStart,false) ) ;

    // current starting point
    offset = rolloutOffset;

    // add all the rollouts
    for (i=0; i<nRolloutList; i++) {

	glui->add_separator(false,0);

	glui->add_control(glui->get_main_panel(),rolloutList[offset],true);
	offset = (offset+1) % nRolloutList;
    }
}


void rotatePanelsOffset(int value)
{
    // any rollout?
    if (!nRolloutList)
	return;

    // rollout offset
    rolloutOffset += value;

    // make sure it is positive
    if (rolloutOffset >= nRolloutList)
        rolloutOffset -= ((rolloutOffset-nRolloutList) / nRolloutList + 1) * nRolloutList;
    if (rolloutOffset < 0)
        rolloutOffset += ((-rolloutOffset) / nRolloutList + 1) * nRolloutList;

    // rotate the rollout panels
    resetPanels();
}


int rotatePanelsTo(char key)
{
    GLUI_Rollout *tmpControl,*tmp2Control;

    int index,i,from,tmpCH,tmp2CH;


    // any rollout?
    if (!nRolloutList)
	return _FALSE;


    // 1. find the index

    // to lowercase
    if ('A' <= key && key <= 'Z')
	key = key-'A'+'a';

    // index to the selected stuff
    index = -1;
    for (i=0; i<nRolloutList; i++) 
	if (rolloutChar[i] == key)
	    index = i;

    // cannot find the key?
    if (index == -1)
	return _FALSE;


    // 2. close the rollout if it is currently open

    if (rolloutList[index]->is_open) {
	rolloutList[index]->close();
	return _TRUE;
    }


    // 3. rearrange rolloutList and rolloutChar consistently

    tmpControl = rolloutList[index];
    tmpCH      = rolloutChar[index];
    from       = rolloutOffset;

    if (from != index) {

	from--;

	do {

	    from = (from+1) % nRolloutList;

	    // swap(tmpControl,rolloutList[from]);
	    tmp2Control       = tmpControl;
	    tmpControl        = rolloutList[from];
	    rolloutList[from] = tmp2Control;

	    // swap(tmpCH,rolloutList[from]);
	    tmp2CH            = tmpCH;
	    tmpCH             = rolloutChar[from];
	    rolloutChar[from] = tmp2CH;

	} while (from != index);
    }


    // 4. reset the panel

    resetPanels();


    // 5. open this panel

    rolloutList[rolloutOffset]->open();


    return _TRUE;
}




//**************************************************************************************//
//                              Main Interface Builder
//**************************************************************************************//

void initGLUI()
{
	int tx,ty,tw,th;

	// create the GLUI panel
	glui = GLUI_Master.create_glui_subwindow(mainWindowID,GLUI_SUBWINDOW_RIGHT);
	glui->bkgd_color.set(GLUI_BGCOLOR);
	glui->fogd_color.set(GLUI_FGCOLOR);

	// create the status bar
	gluiStatusBar = GLUI_Master.create_glui_subwindow(mainWindowID,GLUI_SUBWINDOW_BOTTOM);
	gluiStatusBar->bkgd_color.set(STATUSBAR_BGCOLOR);
	gluiStatusBar->fogd_color.set(STATUSBAR_FGCOLOR);

	// create the controls/widgets
	buildInterface();

	// who is the main window
	glui->set_main_gfx_window( mainWindowID );
	gluiStatusBar->set_main_gfx_window( mainWindowID );

	// reset the window size
	GLUI_Master.get_viewport_area( &tx, &ty, &tw, &th );
	//winW += (winW - tw);
	glutReshapeWindow(winW*2-tw,winH*2-th);
}


static void buildInterface()
{
	///////////////////////////////////////////
	// (1) This is a control panel

	glui->add_separator(false,2);
	glui->add_separator();

	GLUI_StaticText *infoText
		= glui->add_statictext( "CONTROL PANEL" );
	infoText->set_alignment( GLUI_ALIGN_CENTER );

	glui->add_separator();


	// Quit / Help / About

	glui->add_separator(false,-5);

	GLUI_Panel *lastPanel = glui->add_panel( "", GLUI_PANEL_NONE );
	lastPanel->set_w( 220, false );

	//glui->add_column_to_panel(lastPanel,false);
	//GLUI_Button *importButton = glui->add_button_to_panel(lastPanel,"Import",ID_IMPORT,callbackGLUI);
	glui->add_column_to_panel(lastPanel,false);
	GLUI_Button *quitButton = glui->add_button_to_panel(lastPanel,"Quit",ID_QUIT,callbackGLUI);
	glui->add_column_to_panel(lastPanel,false);
	GLUI_Button *helpButton = glui->add_button_to_panel(lastPanel,"Help",ID_HELP,callbackGLUI);
	glui->add_column_to_panel(lastPanel,false);
	GLUI_Button *aboutButton = glui->add_button_to_panel(lastPanel,"About",ID_ABOUT,callbackGLUI);
	glui->add_column_to_panel(lastPanel,false);
	GLUI_Button *menuButton_down = glui->add_button_to_panel(lastPanel,"-",ID_CLOSE_PANELS,callbackGLUI);
	glui->add_column_to_panel(lastPanel,false);
	GLUI_Button *menuButton_up = glui->add_button_to_panel(lastPanel,"m",ID_ROTATE_PANELS,callbackGLUI);

	//importButton->set_w(42, 42);
	quitButton->set_w(42,42);
	helpButton->set_w(42,42);
	aboutButton->set_w(42,42);
	menuButton_up->set_w(12,12);
	menuButton_down->set_w(12,12);

	quitButton->set_alignment(GLUI_ALIGN_CENTER);
	helpButton->set_alignment(GLUI_ALIGN_CENTER);
	aboutButton->set_alignment(GLUI_ALIGN_CENTER);
	menuButton_up->set_alignment(GLUI_ALIGN_CENTER);
	menuButton_down->set_alignment(GLUI_ALIGN_CENTER);

	glui->add_separator(false,-14);

	glui->add_separator(false,0);

	rolloutStart = (glui->get_main_panel())->get_num_childs();


	///////////////////////////////////////////
	// (2) Matcher Panel

	glui->add_separator(false,0);

	GLUI_Rollout *rollout1 = glui->add_rollout("Matcher Control", true);
	rollout1->set_w( 250, true );

	addMatchPanel(rollout1);

	// add to rollout list
	rolloutChar[nRolloutList] = 'b';
	rolloutList[nRolloutList++] = (GLUI_Rollout *) rollout1;


	///////////////////////////////////////////
	// (3) Render Panel

	glui->add_separator(false,0);

	GLUI_Rollout *rollout2 = glui->add_rollout("Render Control", true);
	rollout2->set_w( 250, true );

	addRenderPanel(rollout2);

	// add to rollout list
	rolloutChar[nRolloutList] = 'a';
	rolloutList[nRolloutList++] = (GLUI_Rollout *) rollout2;


	///////////////////////////////////////////
	// (4) Lighting Panel

	glui->add_separator(false,0);

	GLUI_Rollout *rollout3 = glui->add_rollout("Lighting Control",false);
	rollout3->set_w( 250, true );

	addLightPanel(rollout3);

	// add to rollout list
	rolloutChar[nRolloutList] = 'a';
	rolloutList[nRolloutList++] = (GLUI_Rollout *) rollout3;


	///////////////////////////////////////////
	// (4) Status Bar

	statusBar = gluiStatusBar->add_statictext(statusStr);

	statusBar->set_w( winW-10 );
	statusBar->set_alignment(GLUI_ALIGN_LEFT);
}




//**************************************************************************************//
//                                    Add Panels
//**************************************************************************************//

static void addMatchPanel(GLUI_Panel *panel)
{
    glui->add_separator_to_panel(panel, false, 6);

	//////////////////////////////////////////////
	// Edit box

	GLUI_Spinner *startScanIDSpinner = glui->add_spinner_to_panel(panel, "Start ScanID:   ", GLUI_SPINNER_INT, &(startScanID) );
	GLUI_Spinner *endScanIDSpinner   = glui->add_spinner_to_panel(panel, "End  ScanID:   ", GLUI_SPINNER_INT, &(endScanID) );
	GLUI_Spinner *PSRLevelSpinner    = glui->add_spinner_to_panel(panel, "PSR   Level:    ", GLUI_SPINNER_INT, &(PSRLevel) );

	GLUI_Spinner *seedNum_spinner = glui->add_spinner_to_panel(panel, "Seed Number: ",GLUI_SPINNER_INT, &(seedNum) );

	startScanIDSpinner->set_alignment(GLUI_ALIGN_CENTER);
	startScanIDSpinner->set_w( 100 );
	startScanIDSpinner->set_float_limits( 0, 100, GLUI_LIMIT_WRAP );

	endScanIDSpinner->set_alignment(GLUI_ALIGN_CENTER);
	endScanIDSpinner->set_w( 100 );
	endScanIDSpinner->set_float_limits( 0, 100, GLUI_LIMIT_WRAP );

	PSRLevelSpinner->set_alignment(GLUI_ALIGN_CENTER);
	PSRLevelSpinner->set_w( 100 );
	PSRLevelSpinner->set_float_limits( 4, 12, GLUI_LIMIT_WRAP );

	seedNum_spinner->set_alignment(GLUI_ALIGN_CENTER);
	seedNum_spinner->set_w( 100 );
	seedNum_spinner->set_int_limits( 0, 10000000, GLUI_LIMIT_WRAP );


    //////////////////////////////////////////////
    // Button

    glui->add_separator_to_panel(panel, false, -10);

	glui->add_separator_to_panel(panel, false, 4);
	GLUI_Panel *subPanel = glui->add_panel_to_panel(panel, "", GLUI_PANEL_NONE);

	GLUI_Button *readScansButton = glui->add_button_to_panel(subPanel,   "Read Scans",   ID_READ_SCANS,    callbackGLUI);
	glui->add_separator_to_panel(subPanel,false, 4); 
	GLUI_Button *resetScansButton = glui->add_button_to_panel(subPanel,  "Reset Scene",  ID_RESET_SCENE,   callbackGLUI);
	glui->add_separator_to_panel(subPanel,false, 10);
	GLUI_Button *saveMatsButton = glui->add_button_to_panel(subPanel,    "Save .MATS",    ID_SAVE_MATRICES, callbackGLUI);
	glui->add_separator_to_panel(subPanel,false, 4);
	GLUI_Button *saveModelButton = glui->add_button_to_panel(subPanel,   "Save Model",    ID_SAVE_MODEL,    callbackGLUI);
	glui->add_separator_to_panel(subPanel,false, 4);

	glui->add_column_to_panel(subPanel,false);
	GLUI_Button *alignScansButton = glui->add_button_to_panel(subPanel,  "Align Scans",  ID_ALIGN_SCANS,   callbackGLUI);
	glui->add_separator_to_panel(subPanel,false, 4);
	GLUI_Button *funcTestButton = glui->add_button_to_panel(subPanel,    "Func Test",    ID_FUNCTION_TEST, callbackGLUI);
	glui->add_separator_to_panel(subPanel,false, 10);
	GLUI_Button *readMatsButton = glui->add_button_to_panel(subPanel,    "Read .MATS",   ID_READ_MATRICES, callbackGLUI);
	glui->add_separator_to_panel(subPanel,false, 4);
	GLUI_Button *readModelButton = glui->add_button_to_panel(subPanel,   "Read Model",   ID_READ_MODEL,    callbackGLUI);
	glui->add_separator_to_panel(subPanel,false, 4);

	readScansButton->set_alignment(GLUI_ALIGN_LEFT);
	resetScansButton->set_alignment(GLUI_ALIGN_LEFT);
	saveMatsButton->set_alignment(GLUI_ALIGN_LEFT);
	saveModelButton->set_alignment(GLUI_ALIGN_LEFT);
	alignScansButton->set_alignment(GLUI_ALIGN_LEFT);
	funcTestButton->set_alignment(GLUI_ALIGN_LEFT);
	readMatsButton->set_alignment(GLUI_ALIGN_LEFT);
	readModelButton->set_alignment(GLUI_ALIGN_LEFT);
}


static void addRenderPanel(GLUI_Panel *panel)
{
	glui->add_separator_to_panel(panel, false, 6);

	//////////////////////////////////////////////
	// Change FOVY

	glui->add_separator_to_panel(panel, false, -2);
	GLUI_Slider *fovy_slider = glui->add_slider_to_panel(panel, "Fovy:", GLUI_SLIDER_FLOAT, &currFovy, (float) 0.0f, 179.9f, ID_FOVY, callbackGLUI );
	fovy_slider->set_speed( 1.0f );
	fovy_slider->set_w(160,true);
	fovy_slider->set_alignment( GLUI_ALIGN_CENTER );


	//////////////////////////////////////////////
	// Change Background Colors

	glui->add_separator_to_panel(panel, false, 2);

	GLUI_Rollout *subRollout = glui->add_rollout_to_panel(panel, "Background Color", GLUI_PANEL_NONE);
	subRollout->set_w(160,true);

	GLUI_Spinner *backgndR_spinner =
		glui->add_spinner_to_panel( subRollout, "  backgnd red", GLUI_SPINNER_INT, &(backgndColor[0]), ID_CHANGE_BGCOLOR, callbackGLUI );
	backgndR_spinner->set_alignment( GLUI_ALIGN_CENTER );
	backgndR_spinner->set_int_limits( 0, 255 );

	GLUI_Spinner *backgndG_spinner =
		glui->add_spinner_to_panel( subRollout, "backgnd green", GLUI_SPINNER_INT, &(backgndColor[1]), ID_CHANGE_BGCOLOR, callbackGLUI );
	backgndG_spinner->set_alignment( GLUI_ALIGN_CENTER );
	backgndG_spinner->set_int_limits( 0, 255 );

	GLUI_Spinner *backgndB_spinner =
		glui->add_spinner_to_panel( subRollout, " backgnd blue", GLUI_SPINNER_INT, &(backgndColor[2]), ID_CHANGE_BGCOLOR, callbackGLUI );
	backgndB_spinner->set_alignment( GLUI_ALIGN_CENTER );
	backgndB_spinner->set_int_limits( 0, 255 );

	glui->add_separator_to_panel(panel, false, 2);


	//////////////////////////////////////////////
	// Subpanel 1: shwo/hide check box (scans + voxlizer)

    glui->add_separator_to_panel(panel, true,  -2);
	glui->add_separator_to_panel(panel, false,  1);

	GLUI_Panel * subPanel1 = glui->add_panel_to_panel(panel, "", GLUI_PANEL_NONE);
	subPanel1->set_alignment(GLUI_ALIGN_CENTER);

	// Scans flags
	GLUI_Panel * subPanel1Left = glui->add_panel_to_panel( subPanel1 , "Scans");
	glui->add_checkbox_to_panel( subPanel1Left, "show scan ",    &showScan,       ID_SHOW_SCAN,            callbackGLUI );
	glui->add_checkbox_to_panel( subPanel1Left, "show wire ",    &showScanWire,   ID_SHOW_SCAN_WIRE,       callbackGLUI );
	glui->add_checkbox_to_panel( subPanel1Left, "show pick",     &showPickBBox,   ID_SHOW_PICK_BBOX,       callbackGLUI );
	//glui->add_checkbox_to_panel( subPanel1Left, "show bbox",     &showScanBBox,   ID_SHOW_SCAN_BBOX,       callbackGLUI );
	glui->add_checkbox_to_panel( subPanel1Left, "show featPts",  &showFeatPts,    ID_SHOW_FEATURE_POINTS,  callbackGLUI );
	glui->add_checkbox_to_panel( subPanel1Left, "show sampPts",  &showSampPts,    ID_SHOW_SAMPLE_POINTS,   callbackGLUI );
	glui->add_checkbox_to_panel( subPanel1Left, "show evalPts",  &showEvalPts,    ID_SHOW_EVALUATE_POINTS, callbackGLUI );
	glui->add_checkbox_to_panel( subPanel1Left, "show align",    &showAlignScan,  ID_SHOW_ALIGNED_SCAN,    callbackGLUI );

	// Voxelizer flags
	glui->add_column_to_panel(subPanel1, false);
	GLUI_Panel * subPanel1Right = glui->add_panel_to_panel( subPanel1 , "Voxelizer");
	glui->add_checkbox_to_panel( subPanel1Right, "show keyPt",   &showKeyPt,      ID_SHOW_KEY_POINT,         callbackGLUI );
	glui->add_checkbox_to_panel( subPanel1Right, "show SShape",  &showLRFShape,   ID_SHOW_LRF_SHAPE,         callbackGLUI );
	glui->add_checkbox_to_panel( subPanel1Right, "show Sphere",  &showLRFSphere,  ID_SHOW_LRF_SPHERE,        callbackGLUI );
	glui->add_checkbox_to_panel( subPanel1Right, "show LRF",     &showLRF,        ID_SHOW_LRF,               callbackGLUI );
	glui->add_checkbox_to_panel( subPanel1Right, "show LShape",  &showLShape,     ID_SHOW_LOCAL_SHAPE,       callbackGLUI );
	glui->add_checkbox_to_panel( subPanel1Right, "show LSWire",  &showLShapeWire, ID_SHOW_LOCAL_SHAPE_WIRE,  callbackGLUI );
	glui->add_checkbox_to_panel( subPanel1Right, "show LBBox",   &showLBBox,      ID_SHOW_LOCAL_BBOx,        callbackGLUI );
	glui->add_checkbox_to_panel( subPanel1Right, "show LGrid",   &showLGrid,      ID_SHOW_LOCAL_GRID,        callbackGLUI );
	glui->add_checkbox_to_panel( subPanel1Right, "show candis",  &showHighCandis, ID_SHOW_CANDI_MATCHES,     callbackGLUI );


	//////////////////////////////////////////////
	// Subpanel 2: hide/draw check box (model + scene)

	//glui->add_separator_to_panel(panel, true,  -2);
	glui->add_separator_to_panel(panel, false,  -12);

	GLUI_Panel * subPanel2 = glui->add_panel_to_panel(panel, "", GLUI_PANEL_NONE);
	subPanel2->set_alignment(GLUI_ALIGN_CENTER);

	// 3D model flags
	GLUI_Panel * subPanel2Left = glui->add_panel_to_panel( subPanel2 , "Model");
	glui->add_checkbox_to_panel( subPanel2Left,  "show model ",  &showModel,       ID_SHOW_MODEL,       callbackGLUI );
	glui->add_checkbox_to_panel( subPanel2Left,  "show wire ",   &showModelWire,   ID_SHOW_MODEL_WIRE,  callbackGLUI );
	glui->add_checkbox_to_panel( subPanel2Left,  "show bbox",    &showModelBBox,   ID_SHOW_MODEL_BBOX,  callbackGLUI );

	// 3D scene flags
	glui->add_column_to_panel(subPanel2, false);
	GLUI_Panel * subPanel2Right = glui->add_panel_to_panel( subPanel2 , "Scene");
	glui->add_checkbox_to_panel( subPanel2Right, "show axes ",   &showAxes,        ID_SHOW_AXES,         callbackGLUI );
	glui->add_checkbox_to_panel( subPanel2Right, "show ray ",    &showRay,        ID_SHOW_RAY,               callbackGLUI );


	//////////////////////////////////////////////
	// Subpanel 3: draw and control mode

	glui->add_separator_to_panel(panel, true,   2);
	glui->add_separator_to_panel(panel, false, -8);

	GLUI_Panel * subPanel3 = glui->add_panel_to_panel(panel, "", GLUI_PANEL_NONE);
	subPanel3->set_alignment(GLUI_ALIGN_CENTER);

	// Draw mode
	GLUI_Panel * subPanel3Left = glui->add_panel_to_panel( subPanel3 , "Draw mode");
	GLUI_RadioGroup *drawRadioGroup;
	drawRadioGroup = glui->add_radiogroup_to_panel(subPanel3Left, &drawMode, ID_DRAW_MODE, callbackGLUI);
	glui->add_radiobutton_to_group( drawRadioGroup, "Point" );
	glui->add_radiobutton_to_group( drawRadioGroup, "Wireframe" );
	glui->add_radiobutton_to_group( drawRadioGroup, "Flat" );
	glui->add_radiobutton_to_group( drawRadioGroup, "Smooth" );

	// Control mode
	glui->add_column_to_panel(subPanel3, false);
	GLUI_Panel * subPanel3Right = glui->add_panel_to_panel( subPanel3 , "Control mode");
	GLUI_RadioGroup *projectRadioGroup;
	projectRadioGroup = glui->add_radiogroup_to_panel(subPanel3Right, &ctrlMode, ID_CONTROL_MODE, callbackGLUI);
	glui->add_radiobutton_to_group( projectRadioGroup, "Scan" );
	glui->add_radiobutton_to_group( projectRadioGroup, "World" );

	//glui->add_separator_to_panel(basicAdjPanel,false,2);
}

static void addLightPanel(GLUI_Panel *panel)
{
	GLUI_StaticText * infoText ;
	GLUI_Button     * button   ;

	glui->add_separator_to_panel(panel,false,1);

	//////////////////////////////////////////////
	// Ambient color

	GLUI_Rollout *ambSubRollout = glui->add_rollout_to_panel(panel,"Ambient Color",GLUI_PANEL_NONE);
	ambSubRollout->set_w(180,true);

	GLUI_Spinner *ambR_spinner = glui->add_spinner_to_panel(ambSubRollout, "Ambient R: ",GLUI_SPINNER_FLOAT, &(lightAmbient[0]) );
	ambR_spinner->set_alignment(GLUI_ALIGN_CENTER);
	ambR_spinner->set_int_limits( 0, 1.0, GLUI_LIMIT_WRAP );

	GLUI_Spinner *ambG_spinner = glui->add_spinner_to_panel(ambSubRollout, "Ambient G: ",GLUI_SPINNER_FLOAT, &(lightAmbient[1]) );
	ambG_spinner->set_alignment(GLUI_ALIGN_CENTER);
	ambG_spinner->set_int_limits( 0, 1.0, GLUI_LIMIT_WRAP );

	GLUI_Spinner *ambB_spinner = glui->add_spinner_to_panel(ambSubRollout, "Ambient B: ",GLUI_SPINNER_FLOAT, &(lightAmbient[2]) );
	ambB_spinner->set_alignment(GLUI_ALIGN_CENTER);
	ambB_spinner->set_int_limits( 0, 1.0, GLUI_LIMIT_WRAP );

	//////////////////////////////////////////////
	// Diffuse color

	GLUI_Rollout *diffSubRollout = glui->add_rollout_to_panel(panel,"Diffuse Color",GLUI_PANEL_NONE);
	diffSubRollout->set_w(180,true);

	GLUI_Spinner *diffR_spinner = glui->add_spinner_to_panel(diffSubRollout, "Diffuse R: ",GLUI_SPINNER_FLOAT, &(lightDiffuse[0]) );
	diffR_spinner->set_alignment(GLUI_ALIGN_CENTER);
	diffR_spinner->set_int_limits( 0, 1.0, GLUI_LIMIT_WRAP );

	GLUI_Spinner *diffG_spinner = glui->add_spinner_to_panel(diffSubRollout, "Diffuse G: ",GLUI_SPINNER_FLOAT, &(lightDiffuse[1]) );
	diffG_spinner->set_alignment(GLUI_ALIGN_CENTER);
	diffG_spinner->set_int_limits( 0, 1.0, GLUI_LIMIT_WRAP );

	GLUI_Spinner *diffB_spinner = glui->add_spinner_to_panel(diffSubRollout, "Diffuse B: ",GLUI_SPINNER_FLOAT, &(lightDiffuse[2]) );
	diffB_spinner->set_alignment(GLUI_ALIGN_CENTER);
	diffB_spinner->set_int_limits( 0, 1.0, GLUI_LIMIT_WRAP );

	//////////////////////////////////////////////
	// Specular color

	GLUI_Rollout *specSubRollout = glui->add_rollout_to_panel(panel,"Specular Color",GLUI_PANEL_NONE);
	specSubRollout->set_w(180,true);

	GLUI_Spinner *specR_spinner = glui->add_spinner_to_panel(specSubRollout, "Specular R: ",GLUI_SPINNER_FLOAT, &(lightSpecular[0]) );
	specR_spinner->set_alignment(GLUI_ALIGN_CENTER);
	specR_spinner->set_int_limits( 0, 1.0, GLUI_LIMIT_WRAP );

	GLUI_Spinner *specG_spinner = glui->add_spinner_to_panel(specSubRollout, "Specular G: ",GLUI_SPINNER_FLOAT, &(lightSpecular[1]) );
	specG_spinner->set_alignment(GLUI_ALIGN_CENTER);
	specG_spinner->set_int_limits( 0, 1.0, GLUI_LIMIT_WRAP );

	GLUI_Spinner *specB_spinner = glui->add_spinner_to_panel(specSubRollout, "Specular B: ",GLUI_SPINNER_FLOAT, &(lightSpecular[2]) );
	specB_spinner->set_alignment(GLUI_ALIGN_CENTER);
	specB_spinner->set_int_limits( 0, 1.0, GLUI_LIMIT_WRAP );

	//////////////////////////////////////////////
	// Light direction

	GLUI_Panel * subPanel = glui->add_panel_to_panel( panel , "" , GLUI_PANEL_NONE );
	subPanel->set_alignment(GLUI_ALIGN_CENTER);

	glui->add_separator_to_panel(subPanel, false, 6);
	GLUI_Button *resetButton = glui->add_button_to_panel(subPanel, "Reset Light", ID_RESET_LIGHT_DIR, callbackGLUI);
	resetButton->set_w(80,true);
	resetButton->set_alignment(GLUI_ALIGN_LEFT);
	glui->add_separator_to_panel(subPanel, false,6 );
	glui->add_checkbox_to_panel( subPanel,"show light dir", &showLightDir, ID_SHOW_LIGHTDIR, callbackGLUI );

	glui->add_column_to_panel(subPanel,false);
	GLUI_Rotation * rotater = glui->add_rotation_to_panel( subPanel, "Rotation", light_rot_mat , ID_LIGHT_ROT,callbackGLUI  );
}




//**************************************************************************************//
//                                   GLUI Callback
//**************************************************************************************//

void callbackGLUI(int id)
{
    void myExit(int exitCode);
    void myHelp();
    void myAbout();
    void resetStatusBar(const char * fmt, ... );
    void resetProj();

    int winID;


    // Make sure the view is mainWindow
    winID = glutGetWindow();
    if (winID != mainWindowID)
	glutSetWindow(mainWindowID);

    switch(id)
	{
    ////////////////////////////////////////////////////////////////////
    // QUIT / HELP / ABOUT
    ////////////////////////////////////////////////////////////////////

	//case ID_IMPORT:
	//	glutPostRedisplay();
	//	break;

	case ID_QUIT :
		myExit(EXIT_SUCCESS);
		break;

    case ID_HELP :
		myHelp();
		break;

    case ID_ABOUT :
		myAbout();
		break;

    case ID_ROTATE_PANELS :
		rotatePanelsOffset(1);
		resetStatusBar("cycle control panel.");
		break;

    case ID_CLOSE_PANELS :
		closeAllRollouts();
		resetStatusBar("all panels closed.");
		break;


	////////////////////////////////////////////////////////////////////
    // Matcher Adjustment Panel
    ////////////////////////////////////////////////////////////////////

	case ID_READ_SCANS:
		ReadScanFiles();
		preMouseX = -1; // Prevent idle rotation
		preMouseY = -1; 
		glutPostRedisplay();
		break;

	case ID_ALIGN_SCANS:
		myMatcher.AlignAllScans(startScanID, endScanID);
		break;

	case ID_RESET_SCENE:
		myMatcher.ResetScene();
		//myMatcher.Function_RefineAlign();
		break;

	case ID_FUNCTION_TEST:
		 srand( seedNum );
		//myMatcher.Function_Sample(startScanID, endScanID);
		//myMatcher.Function_Match(startScanID, startScanID+1);
		//myMatcher.Function_Evaluate(startScanID, endScanID);
		//myMatcher.Function_MeshClipping();
		//myMatcher.Function_PairwiseAlign();
		//myMatcher.Function_Experiment(startScanID, endScanID);
		 myMatcher.Function_RefineAlign();
		break;

	case ID_SAVE_MATRICES:
		WriteAlignMatrices_MATS();
		break;

	case ID_READ_MATRICES:
		ReadAlignMatrices_MATS();
		glutPostRedisplay();
		break;

	case ID_SAVE_MODEL:
		WriteModelFile();
		break;

	case ID_READ_MODEL:
		ReadModelFile();
		glutPostRedisplay();
		break;


	////////////////////////////////////////////////////////////////////
	// Render Adjustment Panel
	////////////////////////////////////////////////////////////////////

    case ID_FOVY :
		//resetStatusBar("changing fovy.");
		resetProj();
		glutPostRedisplay();
		break;

    case ID_CHANGE_BGCOLOR :
		//resetStatusBar("Background Color changed");
		glClearColor( backgndColor[0]/255.0f, backgndColor[1]/255.0f, backgndColor[2]/255.0f, 0.0f );
		glutPostRedisplay();
		break;


	////////////////////////////////////////////////////////////////////
	// Light Panel
	////////////////////////////////////////////////////////////////////

	case ID_RESET_LIGHT_DIR :
		light_rot_mat[0] = 1.0f ; light_rot_mat[4] = 0.0f ; light_rot_mat[8]  = 0.0f ; light_rot_mat[12] = 0.0f ; 
		light_rot_mat[1] = 0.0f ; light_rot_mat[5] = 1.0f ; light_rot_mat[9]  = 0.0f ; light_rot_mat[13] = 0.0f ; 
		light_rot_mat[2] = 0.0f ; light_rot_mat[6] = 0.0f ; light_rot_mat[10] = 1.0f ; light_rot_mat[14] = 0.0f ; 
		light_rot_mat[3] = 0.0f ; light_rot_mat[7] = 0.0f ; light_rot_mat[11] = 0.0f ; light_rot_mat[15] = 1.0f ; 
		resetStatusBar("Reset lighting direction");
		glutPostRedisplay();
		break;

	case ID_SHOW_LIGHTDIR :
		if ( showLightDir )
			resetStatusBar("Show lighting direction");
		else
			resetStatusBar("Hide lighting direction");
		glutPostRedisplay();
		break;

    }

    if (winID != mainWindowID)
	glutSetWindow(winID);
}




//**************************************************************************************//
//                                 Open and Save File
//**************************************************************************************//

void ReadScanFiles()
{
	HWND hWnd=GetForegroundWindow();
	char objFileName[256];

	if( GetOpenFileName(hWnd, objFileName, 256 ) )
	{
		myMatcher.LoadScans( objFileName );
		resetStatusBar(objFileName);
	}
}

void ReadAlignMatrices_MATS()
{
	HWND hWnd=GetForegroundWindow();
	char matsFileName[FILE_NAME_LENGTH];

	if( GetOpenFileName(hWnd, matsFileName, FILE_NAME_LENGTH ) )
	{
		myMatcher.ReadAlignMatrices( matsFileName );
		resetStatusBar( matsFileName );
	}
}

void WriteAlignMatrices_MATS()
{
	HWND hWnd=GetForegroundWindow();
	char matsFileName[FILE_NAME_LENGTH];

	if( GetSaveFileName(hWnd, matsFileName, FILE_NAME_LENGTH ) )
	{
		myMatcher.WriteAlignMatrices( matsFileName );
		resetStatusBar( matsFileName );
	}
}

void ReadModelFile()
{
	HWND hWnd=GetForegroundWindow();
	char objFileName[256];

	if( GetOpenFileName(hWnd, objFileName, 256 ) )
	{
		myMatcher.ReadModel( objFileName );
		resetStatusBar(objFileName);
	}
}

void WriteModelFile()
{
	HWND hWnd=GetForegroundWindow();
	char objFileName[FILE_NAME_LENGTH];

	if( GetSaveFileName(hWnd, objFileName, FILE_NAME_LENGTH ) )
	{
		myMatcher.PoissonReconstruction(objFileName, PSRLevel);
		resetStatusBar( objFileName );
	}
}

bool GetOpenFileName(HWND hWnd, LPSTR szFile, int StringSize)
{
	OPENFILENAME   ofn;
	char   szFileTitle[256];
	strcpy(szFileTitle, "Open File");
	szFile[0]=0;

	memset(&ofn,0,sizeof(OPENFILENAME));
	ofn.lStructSize     =   sizeof(OPENFILENAME);
	ofn.hwndOwner   =   hWnd;
	ofn.lpstrFilter  = "All\0*.*\0.obj\0*.obj\0";
	ofn.nFilterIndex    =   1;
	ofn.lpstrFile          =   szFile;
	ofn.nMaxFile        =   StringSize;
	ofn.lpstrTitle         =   szFileTitle;
	ofn.Flags              =   OFN_FILEMUSTEXIST;

	if(GetOpenFileName(&ofn)!=TRUE)
	{
		DWORD   Errval;
		char   Errstr[50]= "Common   Dialog   Error:   ";
		char   buf[5];
		Errval=CommDlgExtendedError();
		if(Errval!=0)
		{
			wsprintf(buf, "%ld ",Errval);
			strcat(Errstr,buf);
			MessageBox(NULL,Errstr, "Warning ",MB_OK|MB_ICONSTOP);
		}
		return false;
	}

	return   true;
} 

bool GetSaveFileName(HWND hWnd, LPSTR szFile, int StringSize)
{
	OPENFILENAME   ofn;
	char   szFileTitle[256];
	strcpy(szFileTitle, "Save File");
	szFile[0]=0;

	memset(&ofn,0,sizeof(OPENFILENAME));
	ofn.lStructSize     =   sizeof(OPENFILENAME);
	ofn.hwndOwner   =   hWnd;
	ofn.lpstrFilter  = "All\0*.*\0.obj\0*.obj\0";
	ofn.nFilterIndex    =   1;
	ofn.lpstrFile          =   szFile;
	ofn.nMaxFile        =   StringSize;
	ofn.lpstrTitle         =   szFileTitle;
	ofn.Flags              =   OFN_FILEMUSTEXIST;

	if(GetSaveFileName(&ofn)!=TRUE)
	{
		DWORD   Errval;
		char   Errstr[50]= "Common   Dialog   Error:   ";
		char   buf[5];
		Errval=CommDlgExtendedError();
		if(Errval!=0)
		{
			wsprintf(buf, "%ld ",Errval);
			strcat(Errstr,buf);
			MessageBox(NULL,Errstr, "Warning ",MB_OK|MB_ICONSTOP);
		}
		return false;
	}

	return   true;
} 




//**************************************************************************************//
//                                  Auxiliary Functions
//**************************************************************************************//

void resetStatusBar(const char * fmt, ... )
{
	va_list argp;
	char    str[_MAX_STR_SIZE];


	// 1. string to be put onto the status bar

	va_start(argp, fmt);
	vsprintf(str, fmt, argp);
	va_end(argp);


	// 2. put it on the status bar if it is not on the status bar now

	if (strcmp(str,statusStr)) {
		strcpy(statusStr,str);
		if (statusBar) statusBar->set_text(statusStr);
	}
}


#ifdef WIN32
DWORD WINAPI
openHelpMessageBox()
{
	MessageBox(NULL,_HELP_MESSAGE,programName,MB_ICONINFORMATION|MB_SETFOREGROUND|MB_TASKMODAL);
	return 0;
}

DWORD WINAPI
openAboutMessageBox()
{
	MessageBox(NULL,aboutStr,programName,MB_SETFOREGROUND|MB_TASKMODAL);
	return 0;
}
#endif


void myExit(int exitCode)
{
	///////////////////////////////////////////
	// really quit?

#ifdef WIN32

	int r = MessageBox(NULL,"Really Quit?",programName,MB_ICONEXCLAMATION|MB_YESNO|MB_DEFBUTTON2);
	if (r == IDNO)
		return;

#endif

	exit(exitCode);
}


void myHelp()
{
	// feedback
#ifdef _WIN32

	HANDLE hThread;
	DWORD  IDThread;

	// Create a thread for the message box window
	hThread = CreateThread( NULL,						// no security attributes 
		0,						// use default stack size 
		(LPTHREAD_START_ROUTINE) openHelpMessageBox,	// thread function 
		NULL,						// no thread function argument 
		0,						// use default creation flags 
		&IDThread );					// returns thread identifier 

	// Check the return value for success.
	if (hThread == NULL) fprintf(stderr,"Error: creation menu.\n");

#else

	fprintf(stderr,_HELP_MESSAGE);

#endif
}


void myAbout()
{
#ifdef _WIN32

	HANDLE hThread;
	DWORD  IDThread;

	// Create a thread for the message box window
	hThread = CreateThread( NULL,						// no security attributes 
		0,						// use default stack size 
		(LPTHREAD_START_ROUTINE) openAboutMessageBox,	// thread function 
		NULL,						// no thread function argument 
		0,						// use default creation flags 
		&IDThread );					// returns thread identifier 

	// Check the return value for success.
	if (hThread == NULL) fprintf(stderr,"Error: creation menu.\n");

#else

	// feedback
	resetStatusBar(PROGRAM_ABOUT_ONE_LINE);

#endif
}