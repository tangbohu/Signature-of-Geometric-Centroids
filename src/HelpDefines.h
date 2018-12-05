///////////////////////////////////////////////////////////////
//
// HelpDefines.h
//
//   Common Definitions
//
// by Song Peng ( songpeng@ustc.edu.cn )
// 
// 02/Apr/2015
//
///////////////////////////////////////////////////////////////


#ifndef _HELP_DEFINES_H
#define _HELP_DEFINES_H


/////////////////////////////////////////////////////////////////////
// General Definition
/////////////////////////////////////////////////////////////////////

// File and buffer
#define FILE_NAME_LENGTH        1024
#define BUFFER_LENGTH           256

#define _MAX_STR_SIZE		    512
#define _MAX_LINE_SIZE		    512
#define _MAX_PATH_SIZE		    512

#define _TRUE			        1
#define _FALSE			        0

#define SQUARE_ROOT_TWO         1.414213562373f
#define SQUARE_ROOT_THREE       1.732050807569f

// Max and Min Integer/Float
#define  MIN_INT               -10000000
#define  MAX_INT                10000000
#define  MIN_FLOAT             -10000000.0
#define  MAX_FLOAT              10000000.0

// Float Precision 
#define FLOAT_LARGE_ERROR       0.00001
#define FLOAT_SMALL_ERROR       0.0000001

#ifndef M_PI
#define M_PI			        3.1415926535897932384626433832795
#endif

#ifndef _EPSILON
#define _EPSILON		         1e-7
#endif


/////////////////////////////////////////////////////////////////////
// General Formulas
/////////////////////////////////////////////////////////////////////

#ifndef _MAX
#define _MAX(a,b)		      ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef _MIN
#define _MIN(a,b)		      ( ((a) < (b)) ? (a) : (b) )
#endif

#define _ToRadian(X)		  ((X)/180.0*M_PI)
#define _ToDegree(X)		  ((X)*180.0/M_PI)

#ifndef _IN_BETWEEN
#define _IN_BETWEEN(v,a,b)	  ( ((v) >= (a)) && ((v) <= (b)) )
#endif


#endif
