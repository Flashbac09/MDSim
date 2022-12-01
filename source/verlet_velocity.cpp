#include "../include/GeneralConstruction.h"
#include "../include/verlet_list.h"
#include "../include/P_verlet_list.h"
#include "../include/verlet_velocity.h"
#include "../include/P_verlet_velocity.h"
static double u1[Atom_limit], v1[Atom_limit], w1[Atom_limit];
void run_verlet_velocity(double time, bool process_notice, bool MotionOut)
{
    int i = 0, promark = 0, timepart = 10;
    bool SystemOut = true;
    double t = time;
    auto start = chrono::system_clock::now();

    while (i <= p.nsteps)
    {

        if (i % p.nstep_search == 0)
        {
            verlet_list(t);
        }
        p.temperature = Ek2T(2);
        if (i % p.nstep_output == 0)
        {
            Verlet_Velocity_Output(SystemOut, MotionOut, i);
        }
        if (i != 0)
        {
            // NVT adjustment:
            if (p.thermostat == "Berendson")
            {
                NVT_Berendson(p.temperature, p.NVT_Temperature);
            }
            else if (p.thermostat == "Andersen")
            {
                NVT_Andersen(i);
            }
        }
        if (i >= p.rdf_start && i <= p.rdf_stop && i % p.rdf_interval == 0 && p.cal_rdf == 1)
        {
            RDF(i);
        }
        // t->t+dt:
        i += 1;
        t += p.dt;
        for (int j = 1; j <= p.atomsum; j++) // ev:e.kg.m.m/s.s   f:du/dr*dr/dx--ev/Am    r:Am     f= e.kg.10^10.m/(s.s)
        {
            p.x[j] += p.u[j] * p.dt + 1000 * p.fx[j] * p.dt * p.dt * Avogadro_NA * ElementaryCharge / (2 * p.molar_mass); // f*(dt*dt)/2molarmass=e.Am.10^-4.NA/(2)  =e.NA.Am
            p.y[j] += p.v[j] * p.dt + 1000 * p.fy[j] * p.dt * p.dt * Avogadro_NA * ElementaryCharge / (2 * p.molar_mass); // ev*ps*ps*mol/(Am*Kg)=(e*NA) (kg*m*m*ps*ps)/(s*s*AM*kg)=(e*NA)(Am*10^-4)
            p.z[j] += p.w[j] * p.dt + 1000 * p.fz[j] * p.dt * p.dt * Avogadro_NA * ElementaryCharge / (2 * p.molar_mass);
            // p.u[j] += p.dt * 1000 * p.fx[j] * Avogadro_NA * ElementaryCharge / (2 * p.molar_mass);
            // p.v[j] += p.dt * 1000 * p.fy[j] * Avogadro_NA * ElementaryCharge / (2 * p.molar_mass);
            // p.w[j] += p.dt * 1000 * p.fz[j] * Avogadro_NA * ElementaryCharge / (2 * p.molar_mass);
            // p.x[j] += p.dt * p.u[j];
            // p.y[j] += p.dt * p.v[j];
            // p.z[j] += p.dt * p.w[j];
        }
        memmove(u1, p.fx, sizeof(double) * (p.atomsum + 10));
        memmove(v1, p.fy, sizeof(double) * (p.atomsum + 10));
        memmove(w1, p.fz, sizeof(double) * (p.atomsum + 10));//为了处理Berendson热库的速度信息迭代，不得不使用最老实的VerletVelocity公式
        // R2FP():位移到势能和力场

        memset(p.LJ_Pot, 0, (p.atomsum + 10) * sizeof(double));
        memset(p.fx, 0, (p.atomsum + 10) * sizeof(double));
        memset(p.fy, 0, (p.atomsum + 10) * sizeof(double));
        memset(p.fz, 0, (p.atomsum + 10) * sizeof(double));
        for (int i = 1; i <= p.atomsum; i++)
        {
            for (int k = 1; k <= p.nlist[i]; k++)
            {
                int j = p.list[i][k];
                if (j > i)
                {
                    double rx = 0, ry = 0, rz = 0, dist = 0, pot = 0, unitforce = 0;
                    rx = p.x[j] - p.x[i];
                    while (rx >= p.Xborder / 2)
                    {
                        rx -= p.Xborder;
                    }
                    while (rx <= -p.Xborder / 2)
                    {
                        rx += p.Xborder;
                    }
                    if (fabs(rx) > (p.r_cut + p.neighbor_r))
                        continue;

                    ry = p.y[j] - p.y[i];
                    while (ry >= p.Yborder / 2)
                    {
                        ry -= p.Yborder;
                    }
                    while (ry <= -p.Yborder / 2)
                    {
                        ry += p.Yborder;
                    }
                    if (fabs(ry) > (p.r_cut + p.neighbor_r))
                        continue;

                    rz = p.z[j] - p.z[i];
                    while (rz >= p.Zborder / 2)
                    {
                        rz -= p.Zborder;
                    }
                    while (rz <= -p.Zborder / 2)
                    {
                        rz += p.Zborder;
                    }
                    if (fabs(rz) > (p.r_cut + p.neighbor_r))
                        continue;
                    else
                    {
                        dist = sqrt(rx * rx + ry * ry + rz * rz);
                        if (dist > (p.r_cut + p.neighbor_r))
                            continue;
                        else if (dist <= p.r_cut)
                        {

                            pot = 4 * p.epsilon * (pow(p.sigma / dist, 12) - pow(p.sigma / dist, 6)) - 4 * p.epsilon * (pow(p.sigma / p.r_cut, 12) - pow(p.sigma / p.r_cut, 6));
                            p.LJ_Pot[i] += pot;
                            p.LJ_Pot[j] += pot;
                            unitforce = (-24) * p.epsilon * ((2 * pow(p.sigma / dist, 12) - pow(p.sigma / dist, 6)) / (dist * dist));
                            p.fx[i] += unitforce * rx;
                            p.fy[i] += unitforce * ry;
                            p.fz[i] += unitforce * rz;
                            p.fx[j] += -unitforce * rx;
                            p.fy[j] += -unitforce * ry;
                            p.fz[j] += -unitforce * rz;
                        }
                    }
                }
            }
        }
        for (int j = 1; j <= p.atomsum; j++)
        {
            p.u[j] += 1000 * (p.fx[j] + u1[j]) * p.dt * Avogadro_NA * ElementaryCharge / (2 * p.molar_mass);
            p.v[j] += 1000 * (p.fy[j] + v1[j]) * p.dt * Avogadro_NA * ElementaryCharge / (2 * p.molar_mass);
            p.w[j] += 1000 * (p.fz[j] + w1[j]) * p.dt * Avogadro_NA * ElementaryCharge / (2 * p.molar_mass);
        }

        // Time Notice:
        if (process_notice == true)
        {
            if (i % ((int)(p.nsteps / timepart)) == 0 && promark <= timepart)
            {

                promark += 1;
                cout << "Process: " << setiosflags(ios::fixed) << setprecision(2) << (double)promark / timepart * 100 << "%" << endl;
            }
        }
    }

    auto end = chrono::system_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(end - start);
    cout << "End: MDsim (" << setiosflags(ios::fixed) << setprecision(6) << (double)chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000000 << "s)" << endl;
}

double Ek2T(int im)
{
    double eksum = 0;
    for (int i = 1; i <= p.atomsum; i++)
    {
        eksum += 0.5 * p.molar_mass * (p.u[i] * p.u[i] + p.v[i] * p.v[i] + p.w[i] * p.w[i]);
    }
    // m-/na v^2-*10000 __Joule

    if (im == 1)
        return eksum / (1000 * ElementaryCharge * Avogadro_NA); // joule---eV,kg-g,Am-m,ps-s
    else
        return 10 * eksum / (Avogadro_NA * 1.5 * p.atomsum * Boltzmann_K);
}
double Ep()
{
    double epsum = 0;
    for (int i = 1; i <= p.atomsum; i++)
    {
        epsum += p.LJ_Pot[i];
    }
    return 0.5 * epsum;
}

void NVT_Berendson(double T, double T0)
{
    double lambda, dT;
    dT = (p.dt / p.tau) * (T0 - T);
    lambda = sqrt((T + dT) / T);
    for (int i = 1; i <= p.atomsum; i++)
    {
        p.u[i] *= lambda;
        p.v[i] *= lambda;
        p.w[i] *= lambda;
    }
}

void NVT_Andersen(double step)
{
    double rand1, sum;
    default_random_engine g1((time(0) - 100) * step);
    uniform_real_distribution<double> uniform_real_dis(0.0, 1.0);
    default_random_engine g2((time(0) + 100) * step);
    normal_distribution<double> normal_dis(0, sqrt(Boltzmann_K * Avogadro_NA * p.NVT_Temperature / (10 * p.molar_mass)));

    for (int i = 1; i <= p.atomsum; i++)
    {

        rand1 = uniform_real_dis(g1);
        if (rand1 < 1.0 / p.nraise)
        {

            p.u[i] = normal_dis(g2);
            p.v[i] = normal_dis(g2);
            p.w[i] = normal_dis(g2);
        }
    }
}

void random_speed()
{
    default_random_engine g1((time(0) - 100));
    uniform_real_distribution<double> uniform_real_dis(-0.5, 0.5);
    for (int i = 1; i <= p.atomsum; i++)
    {
        p.u[i] = uniform_real_dis(g1);
        p.v[i] = uniform_real_dis(g1);
        p.w[i] = uniform_real_dis(g1);
    }
    double vx = 0, vy = 0, vz = 0;
    for (int i = 1; i <= p.atomsum; i++)
    {
        vx += p.u[i];
        vy += p.v[i];
        vz += p.w[i];
    }
    vx /= p.atomsum;
    vy /= p.atomsum;
    vz /= p.atomsum;
    for (int i = 1; i <= p.atomsum; i++)
    {
        p.u[i] -= vx;
        p.v[i] -= vy;
        p.w[i] -= vz;
    }
    double ek = 0, factor = 0;
    for (int i = 1; i <= p.atomsum; i++)
    {
        ek += 0.5 * p.molar_mass * (p.u[i] * p.u[i] + p.v[i] * p.v[i] + p.w[i] * p.w[i]);
    }
    factor = sqrt(1.5 * Boltzmann_K * Avogadro_NA * p.atomsum * p.initial_temperature / (ek * 10));
    for (int i = 1; i <= p.atomsum; i++)
    {
        p.u[i] *= factor;
        p.v[i] *= factor;
        p.w[i] *= factor;
    }
    p.temperature = Ek2T(2);
}

void RDF(int step)
{
    double dis = 0;
    int j = 0;
    memset(p.n, 0, sizeof(int) * Atom_limit * 1000);
    for (int i = 1; i <= p.atomsum; i++)
    {
        for (int w = 1; w <= p.nlist[i]; w++)
        {
            j = p.list[i][w];
            if (j > i)
            {
                dis = get_distance(i, j);
                if (dis <= p.rdf_rcut)
                {
                    p.n[i][(int)(dis / p.rdf_dr)] += 1;
                    p.n[j][(int)(dis / p.rdf_dr)] += 1;
                }
            }
        }
    }
    for (int i = 0; i < (int)(p.rdf_rcut / p.rdf_dr); i++)
    {
        for (int j = 1; j <= p.atomsum; j++)
        {
            p.g_rdf[step][i] += p.n[j][i];
        }
        if (i != 0)
            p.g_rdf[step][i] /= (4 * PAI * i * i * p.rdf_dr * p.rdf_dr * p.rdf_dr * p.atomsum * p.atomsum / (p.Xborder * p.Yborder * p.Zborder));
    }
    return;
}

void RDF_out()
{
    for (int j = 0; j < p.rdf_rcut / p.rdf_dr; j++)
    {
        for (int i = p.rdf_start; i <= p.rdf_stop; i++)
        {
            p.gfinal[j] += p.g_rdf[i][j];
        }
        p.gfinal[j] /= (p.rdf_stop - p.rdf_start + 1);
    }
    ofstream fout;
    fout.open("../product/rdf1.txt");
    for (int j = 0; j < p.rdf_rcut / p.rdf_dr; j++)
    {
        fout << std::left << setw(10) << setiosflags(ios::fixed) << setprecision(3) << (j + 0.5) * p.rdf_dr
             << setiosflags(ios::fixed) << setprecision(12) << p.gfinal[j] << endl;
    }
    fout.close();
}