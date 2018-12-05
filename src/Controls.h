///////////////////////////////////////////////////////////////
//
// Controls.h
//
//   Defined Variables to Control the Program
//
// by Song Peng ( songpeng@ustc.edu.cn )
// 
// 11/May/2015
//
///////////////////////////////////////////////////////////////


#ifndef _CONTROLS
#define _CONTROLS


//**************************************************************************************//
//                                Control Parameters 
//**************************************************************************************//

// Input Mesh Model
#define INPUT_HAS_NORMAL              1  

// Mesh Vertex Sampling 
#define SAMPLE_POINT_NUM              2500    // Sample points for description
#define SAMPLE_POINT_DIST             0.017   
//#define FEATURE_POINT_NUM           100     // Sample points for matching (feature)      
//#define FEATURE_POINT_DIST          0.1
#define FEATURE_POINT_NUM             200     // Sample points for matching (feature)      
#define FEATURE_POINT_DIST            0.06
//#define FEATURE_POINT_NUM             600     // Sample points for matching (for the overlap experiment)      
//#define FEATURE_POINT_DIST            0.04

// Local Voxelization
#define LOCAL_VOLUME_DIMEN           Vector3i(8,8,8)
//#define LOCAL_VOLUME_DIMEN           Vector3i(10,10,10)
//#define LOCAL_VOLUME_DIMEN             Vector3i(12,12,12)
//#define LOCAL_SHAPE_RADIUS             0.08    // TODO: need to tune it (e.g., 0.05)
#define LOCAL_SHAPE_RADIUS             0.10    // TODO: need to tune it (e.g., 0.05)

// Feature Selected for Each Voxel
#define FEATURE_TYPE_ID                1

// Scan Registration
#define SCAN_ALIGN_SCORE_THRES         0.50   
#define INIT_ALIGN_SCORE_THRES         0.10
#define DESCRIPTOR_DIST_THRES          6.0
#define DESCRIPTOR_TOP_MATCH_NUM       5
#define KDTREE_MAX_DIST_TIMES          3




//**************************************************************************************//
//                                   Descriptor Comparison
//**************************************************************************************//

// Select Descriptor (by uncommenting it)
//#define USE_DESCRIPTOR_OURS              1
//#define USE_DESCRIPTOR_3DSC              2
//#define USE_DESCRIPTOR_SHOT              3
#define USE_DESCRIPTOR_SGC              4



//**************************************************************************************//
//                            Defined Values for Application 
//**************************************************************************************//

// ID of Features or Their Combinations
#define FEATURE_SURFACE_AREA              1
#define FEATURE_SURFACE_CRNTER            2
#define FEATURE_SURFACE_NORMAL            3
#define FEATURE_SURFACE_CURVATURE         4
#define FEATURE_SHAPE_VOLUME              5
#define FEATURE_SHAPE_CENTER              6
#define FEATURE_COMB_AREA_SURFCEN         7
#define FEATURE_COMB_AREA_NORMAL          8
#define FEATURE_COMB_VOLUME_VOLCEN        9
//#define FEATURE_COMB_SURFCEN_NORMAL       10
//#define FEATURE_COMB_AREA_SURFCEN_NORM    11

// Manipulation Mode
#define MANIPU_SCAN                       0
#define MANIPU_WORLD                      1




//**************************************************************************************//
//                              Default Variable Values
//**************************************************************************************//

// Descriptor Related
//#define SHOW_LRF_ALIGNED_SCAN          1


// Voxel Related
#define  DEFAULT_VOXEL_VALUE           0.0
#define  VOXEL_EDGE_SAMPLE_NUM         5     // For building point grid

// Voxel State
#define  VOXEL_OUT_MESH                0
#define  VOXEL_CROSS_MESH              1
//#define  VOXEL_IN_MESH               2

// Spherical Bin State
#define  BIN_OUT_MESH                  0
#define  BIN_CROSS_MESH                1

// Sample point state
#define  POINT_UNKNOWN                -1
#define  POINT_OUT_MESH                0
#define  POINT_IN_MESH                 1

// Initial World Pose
#define INIT_WORLD_POSITION            Vector3f(0.0,  0.0,  -6.5)
#define INIT_WORLD_ROT_AXIS            Vector3f(0.0,  0.0,  1.0)
#define INIT_WORLD_ROT_ANGLE           90.0

// Initial World Axes Pose
#define INIT_WORLD_AXES_POSITION       Vector3f(0.0,  0.0, -4.5)
#define INIT_WORLD_AXES_ROT_AXIS       Vector3f(0.0,  0.0,  1.0)
#define INIT_WORLD_AXES_ROT_ANGLE      90.0

// Initial Scan Pose
#define INITIAL_SCAN_SHIFT_VECTOR              Vector3f(0.0,  1.0, 0.0)
#define ALIGNED_SCAN_SHIFT_VECTOR             Vector3f(0.0, -1.0, 0.0)


#endif
