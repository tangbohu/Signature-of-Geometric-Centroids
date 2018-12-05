///////////////////////////////////////////////////////////////
//
// gluiMouse.h
//
// - handle mouse event and idle callback
//
// by Song Peng ( song0083@ntu.edu.sg )
//
// 24/Apr/2015
//
//
///////////////////////////////////////////////////////////////


#ifndef _ND_MOUSE_H
#define _ND_MOUSE_H


////////////////////////////////////////////
// Mouse Button and Motion Callback

extern void mouse  ( int button, int state, int x, int y ) ;
extern void motion ( int x, int y ) ;


////////////////////////////////////////////
// Idle Callback

extern void idle();
extern void stopIdleMotion();



#endif
