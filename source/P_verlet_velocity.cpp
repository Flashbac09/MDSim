#include "../include/GeneralConstruction.h"
#include "../include/verlet_velocity.h"
#include "../include/P_verlet_velocity.h"



void P_run_verlet_velocity_energy_log(int step)
{
    ofstream fout;
    if (step == 0)
        fout.open("../product/run.log");
    else
        fout.open("../product/run.log", ios::app);
    if (step == 0)
    {
        fout << "Step    Time(ps)    SystemEnergy(eV)    KineticEnergy(eV)    PotentialEnergy(eV)    Temperature(K)" << endl;
    }
    int precision = 12;
    fout << std::left << setw(10) << step << std::left << setw(10) << setiosflags(ios::fixed) << setprecision(3) << p.dt * step
         << setiosflags(ios::fixed) << setprecision(precision) << Ep() + Ek2T(1) << "     "
         << setiosflags(ios::fixed) << setprecision(precision) << Ek2T(1) << "     "
         << setiosflags(ios::fixed) << setprecision(precision) << Ep() << "     "
         << setiosflags(ios::fixed) << setprecision(precision) << p.temperature << endl;
    fout.close();
    return;
}

void P_run_verlet_velocity_motion_log(int step)
{
    // position
    ofstream fout;
    if (step == 0)
        fout.open("../product/position.txt");
    else
        fout.open("../product/position.txt", ios::app);
    fout << "Step: " << step << "   "
         << "Time: " << (double)(p.dt * step) << endl;
    int precision = 12;
    for (int i = 1; i <= p.atomsum; i++)
    {
        fout << setiosflags(ios::fixed) << setprecision(precision) << p.x[i] << "    " << setiosflags(ios::fixed) << setprecision(precision) << p.y[i] << "    "
             << setiosflags(ios::fixed) << setprecision(precision) << p.z[i] << endl;
    }
    fout << endl;
    fout.close();
    // velocity

    if (step == 0)
        fout.open("../product/velocity.txt");
    else
        fout.open("../product/velocity.txt", ios::app);
    fout << "Step: " << step << "   "
         << "Time: " << (double)(p.dt * step) << endl;
    for (int i = 1; i <= p.atomsum; i++)
    {
        fout << setiosflags(ios::fixed) << setprecision(precision) << p.u[i] << "    " << setiosflags(ios::fixed) << setprecision(precision) << p.v[i] << "    "
             << setiosflags(ios::fixed) << setprecision(precision) << p.w[i] << endl;
    }
    fout << endl;
    fout.close();
    // force
    if (step == 0)
        fout.open("../product/force.txt");
    else
        fout.open("../product/force.txt", ios::app);
    fout << "Step: " << step << "   "
         << "Time: " << (double)(p.dt * step) << endl;
    for (int i = 1; i <= p.atomsum; i++)
    {
        fout << setiosflags(ios::fixed) << setprecision(precision) << p.fx[i] << "    " << setiosflags(ios::fixed) << setprecision(precision) << p.fy[i] << "    "
             << setiosflags(ios::fixed) << setprecision(precision) << p.fz[i] << endl;
    }
    fout << endl;
    fout.close();
    return;
}

void Verlet_Velocity_Output(bool system, bool motion, int step)
{
    if (system == false && motion == false)
        return;
    else if (system == true && motion == false)
    {
        P_run_verlet_velocity_energy_log(step);
        return;
    }
    else if (system == false && motion == true)
    {
        P_run_verlet_velocity_motion_log(step);
        return;
    }
    else
    {
        P_run_verlet_velocity_energy_log(step);
        P_run_verlet_velocity_motion_log(step);
        return;
    }
}