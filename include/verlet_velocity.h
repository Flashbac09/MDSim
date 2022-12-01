#ifndef _VERLET_VELOCITY_
#define _VERLET_VELOCITY_

extern void run_verlet_velocity(double time, bool process_notice, bool MotionOut);
extern void R2FP();
extern double Ek2T(int im);
extern double Ep();

extern void NVT_Berendson(double T, double T0);
extern void NVT_Andersen(double step);

extern void random_speed();

extern void RDF(int step);
extern void RDF_out();
#endif