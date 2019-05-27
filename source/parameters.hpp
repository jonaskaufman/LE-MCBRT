#ifndef PARAMETERS_H
#define PARAMETERS_H

//#define N 1000      /// pixels per side
#define D 10.0     /// source distance
#define MEAN 0.0	/// source angle mean
#define SIGMA 1.0   /// source angle std dev in radians
#define E0 10.0     /// initial primary ray energy
#define A 1.0       /// primary ray interaction probability scaling
#define F 0.5       /// primary ray interaction energy fraction
#define G 0.1       /// secondary ray deposition constant (probably not necessary)
#define KS 4		// number of secondary rays to spawn after interaction
#define MIN_ENERGY  0.00001 // minimum energy a secondary ray can have before dying
							// useful to avoid comparing to 0 which can cause error

/// Global definitions (maybe make these constructor arguments?)
//#define PI 3.14159265358979323846 /// use standard def?
 /// make global definitions file??
 // I'm treating this as a global definitions file

 // Label the pixel edges
enum class PIXEL_EDGE {TOP, BOTTOM, LEFT, RIGHT};


/* Conditional print macros for debugging */
#define DEBUG(debugCode, action) {\
	if (CUR_DEBUG == (debugCode)  || CUR_DEBUG == DEBUG_ALL){\
		action;\
	}\
}




/* Debug codes */
#define CUR_DEBUG DEBUG_ALL
#define DEBUG_ALL -1

#define NO_DEBUG 0
#define DB_ARGPASS 1 // passing arguments to run_serial
#define DB_SIMCONST 2 // building the simulation
#define DB_INITPRIM 3 // initializing primary rays
#define DB_TRACE 4    // tracing ray from pixel to pixel
#define DB_INTERACT 5  // probabilitisticly determining an interaction point
#define DB_SECONDARY 6 // spawn and evolve secondary rays to completion

#endif