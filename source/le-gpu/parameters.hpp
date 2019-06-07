#ifndef PARAMETERS_H
#define PARAMETERS_H

/// Simulation parameters
#define PARAM_D 1.0     /// source distance as relative to grid size
#define PARAM_MEAN 0.0  /// source angle mean
#define PARAM_SIGMA 0.1 /// source angle std dev in radians
#define PARAM_E0 100.0  /// initial primary ray energy
#define PARAM_A 2.0     /// primary ray interaction probability scaling
#define PARAM_F 0.05    /// primary ray interaction energy deposit fraction
#define PARAM_G 1.0     /// secondary ray deposition constant
#define PARAM_KS 10     /// number of secondary rays to spawn after interaction
#define PARAM_MINERGY                                                                                                  \
    0.000000001 /// minimum energy a secondary ray can have before dying
                /// useful to avoid comparing to 0 which can cause error
#define PARAM_EPSILON 0.00000001 // for comparing equality of angle and 2*pi

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
        if (CUR_DEBUG == debugCode || CUR_DEBUG == DB_ALL)                                                             \
        {                                                                                                              \
            action;                                                                                                    \
        }                                                                                                              \
        else if (CUR_DEBUG == DB_GENERAL && (debugCode == DB_INIT_PRI || debugCode == DB_INIT_SEC))                    \
        {                                                                                                              \
            action;                                                                                                    \
        }                                                                                                              \
    }

/* Debug codes */
#define CUR_DEBUG DB_GPU

#define DB_ALL -1
#define DB_NONE 0

#define DB_GENERAL 1    // general progress
#define DB_INIT_PRI 2   // initialize primary ray
#define DB_INIT_SEC 3   // initialize secondary ray
#define DB_EVOLVE_PRI 4 // evolve primary ray
#define DB_EVOLVE_SEC 5 // evolve secondary ray
#define DB_TRACE 6      // trace ray from pixel to pixel
#define DB_GPU 7        // base GPU version

#endif

