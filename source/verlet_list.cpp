#include "../include/GeneralConstruction.h"
#include "../include/verlet_list.h"
void verlet_list(double t)
{
    memset(p.nlist, 0, sizeof(int) * (p.atomsum + 10));
    memset(p.list, 0, sizeof(int) * Atom_limit * Neighbor_limit);
    // memset(p.LJ_Pot, 0, (p.atomsum + 1) * sizeof(double));
    // memset(p.fx, 0, (p.atomsum + 1) * sizeof(double));
    // memset(p.fy, 0, (p.atomsum + 1) * sizeof(double));
    // memset(p.fz, 0, (p.atomsum + 1) * sizeof(double));
    for (int i = 1; i <= p.atomsum; i++)
    {
        for (int j = i + 1; j <= p.atomsum; j++)
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
                else
                {
                    p.nlist[i] += 1;
                    p.nlist[j] += 1;
                    p.list[i][p.nlist[i]] = j;
                    p.list[j][p.nlist[j]] = i;
                    if (dist <= p.r_cut && t == 0)
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
    // cout<<"Done: Verlet List Update, time = "<<time<<" s."<<endl;
    return;
}