

//#define N 1000      /// pixels per side
#define D 100.0     /// source distance
#define SIGMA 1.0   /// source angle std dev in radians
#define E0 10.0     /// initial primary ray energy
#define A 1.0       /// primary ray interaction probability scaling
#define F 0.5       /// primary ray interaction energy fraction
#define G 0.1       /// secondary ray deposition constant (probably not necessary)

/// Global definitions (maybe make these constructor arguments?)
#define PI 3.14159265358979323846 /// use standard def?
 /// make global definitions file??
 // I'm treating this as a global definitions file

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
#define DB_ARGPASS 1
#define DB_SIMCONST 2
