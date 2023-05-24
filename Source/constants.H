#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

#define numrigidbodies 0
#define XDIR 0
#define YDIR 1
#define ZDIR 2

#define fourbythree 1.3333333333
#define half 0.5
#define PI 3.141592654
#define three 3.0
#define two 2.0
#define one 1.0
#define four 4.0
#define twobyfive 0.4
#define fourth 0.25
#define eighth 0.125
#define ZERO 0.0
#define TINYVAL 1e-20
#define BIGVAL 1e20

#define NCOMP_TENSOR 6
#define NCOMP_FULLTENSOR 9

#define NUM_STATES 20
#define MASS_INDEX 0
#define VELX_INDEX 1
#define VELY_INDEX 2
#define VELZ_INDEX 3
#define FRCX_INDEX 4
#define FRCY_INDEX 5
#define FRCZ_INDEX 6
#define DELTA_VELX_INDEX 7		//Difference in velocity between two successive time steps Vx(t+dt)-Vx(t)
#define DELTA_VELY_INDEX 8
#define DELTA_VELZ_INDEX 9
#define MASS_OLD_INDEX 10
#define VELX_RIGID_INDEX 11
#define VELY_RIGID_INDEX 12
#define VELZ_RIGID_INDEX 13
#define MASS_RIGID_INDEX 14
#define STRESS_INDEX 15
#define RIGID_BODY_ID 16
#define NORMALX 17
#define NORMALY 18
#define NORMALZ 19



#define XX 0
#define XY 1
#define XZ 2
#define YY 3
#define YZ 4
#define ZZ 5

#define BC_PERIODIC 0
#define BC_NOSLIPWALL 1
#define BC_SLIPWALL 2
#define BC_PARTIALSLIPWALL 3
#define BC_OUTFLOW 4

#endif
