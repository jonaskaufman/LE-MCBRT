#ifndef PARAMETERS_H
#define PARAMETERS_H

/// Simulation parameters
#define PARAM_D 10.0    /// source distance
#define PARAM_MEAN 0.0  /// source angle mean
#define PARAM_SIGMA 0.5 /// source angle std dev in radians
#define PARAM_E0 1000.0 /// initial primary ray energy
#define PARAM_A 3.0     /// primary ray interaction probability scaling
#define PARAM_F 0.25    /// primary ray interaction energy fraction
#define PARAM_G 1.0     /// secondary ray deposition constant
#define PARAM_KS 10     /// number of secondary rays to spawn after interaction
#define PARAM_MINERGY                                                                                                  \
    0.000000001 /// minimum energy a secondary ray can have before dying
                /// useful to avoid comparing to 0 which can cause error

/// Label the pixel edges
enum class PIXEL_EDGE
{
    TOP,
    BOTTOM,
    LEFT,
    RIGHT
};

/* Conditional print macros for debugging */
#define DEBUG(debugCode, action)                                                                                       \
    {                                                                                                                  \
        if (CUR_DEBUG == (debugCode) || CUR_DEBUG == DEBUG_ALL)                                                        \
        {                                                                                                              \
            action;                                                                                                    \
        }                                                                                                              \
    }

/* Debug codes */
#define CUR_DEBUG DB_EVOLVE_SEC
#define DEBUG_ALL -1

#define NO_DEBUG 0
#define DB_ARGPASS 1    // passing arguments to run_serial
#define DB_SIMCONST 2   // building the simulation
#define DB_INITPRIM 3   // initializing primary rays
#define DB_TRACE 4      // tracing ray from pixel to pixel
#define DB_INTERACT 5   // probabilitisticly determining an interaction point
#define DB_SECONDARY 6  // spawn and evolve secondary rays to completion
#define DB_INITSECOND 7 // spawn secondary rays
#define DB_EVOLVE_PRI 8
#define DB_EVOLVE_SEC 9
#endif
