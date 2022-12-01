#ifndef _P_VERLET_VELOCITY_
#define _P_VERLET_VELOCITY_

extern void P_run_verlet_velocity_energy_log(int step);
extern void P_run_verlet_velocity_motion_log(int step);
extern void Verlet_Velocity_Output(bool system,bool motion,int step);

#endif