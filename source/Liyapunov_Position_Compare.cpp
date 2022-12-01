#include "../include/GeneralConstruction.h"
#include "../include/initial_input.h"
#include "../include/geo_input.h"
#include "../include/verlet_velocity.h"
#include "../include/Liyapunov_Position_Compare.h"

void Liyapunov_Position_Compare()
{
    double r1[Atom_limit][5] = {0}, r2[Atom_limit][5] = {0}, compare[Atom_limit] = {0};
    ifstream fin1, fin2;
    string temp1, temp2;

    fin1.open("../product/position_liyapunov1.txt", ios::in);
    fin2.open("../product/position.txt", ios::in);
    ofstream fout,fout2;
    fout.open("../product/Liyapunov_compare_x");
    fout2.open("../product/Liyapunov_compare_x1");
    for (int i = 0; i <= p.nsteps; i += p.nstep_output)
    {
        // 1
        getline(fin1, temp1);
        for (int j = 1; j <= p.atomsum; j++)
        {
            fin1 >> r1[j][0] >> r1[j][1] >> r1[j][2];
        }

        getline(fin1, temp1);
        getline(fin1, temp1);

        // 2
        getline(fin2, temp2);
        for (int j = 1; j <= p.atomsum; j++)
        {
            fin2 >> r2[j][0] >> r2[j][1] >> r2[j][2];
        }
        getline(fin2, temp2);
        getline(fin2, temp2);
        // delta sum
        for (int j = 1; j <= p.atomsum; j++)
        {
            compare[i / p.nstep_output] += pow(r1[j][0] - r2[j][0], 2) + pow(r1[j][1] - r2[j][1], 2) + pow(r1[j][2] - r2[j][2], 2);
        }
        fout<< std::left << setw(10) <<  i << setiosflags(ios::fixed) << setprecision(12) << r1[1][0] << endl;
        fout2<< std::left << setw(10) <<  i << setiosflags(ios::fixed) << setprecision(12) << r2[1][0] << endl;
    }
    fin1.close();
    fin2.close();
    fout.close();
    fout2.close();

    fout.open("../product/Liyapunov_Position_Compare.log");
    fout << "Time(ps)    Sum of DistDiff^2(Am)" << endl;
    for (int i = 0; i <= p.nsteps; i += p.nstep_output)
    {
        if (compare[i / p.nstep_output] >= 1000000000)
        {
            cout << "Liyapunov Test End: " << i * p.dt << "s" << endl;
            fout << "Liyapunov Test End: " << i * p.dt << "s" << endl;
            break;
        }
        fout << std::left << setw(10) << setiosflags(ios::fixed) << setprecision(3) << p.dt * i << setiosflags(ios::fixed) << setprecision(15) << compare[i / p.nstep_output] << endl;
    }
    fout.close();
    return;
}

void Liyapunov_Test()
{  
    memset(p.nlist, 0, sizeof(int) * (p.atomsum + 10));
    memset(p.list, 0, sizeof(int) * Atom_limit * Neighbor_limit);
     memset(p.LJ_Pot, 0, (p.atomsum + 1) * sizeof(double));
     memset(p.fx, 0, (p.atomsum + 1) * sizeof(double));
     memset(p.fy, 0, (p.atomsum + 1) * sizeof(double));
     memset(p.fz, 0, (p.atomsum + 1) * sizeof(double));
    if (access("../product/run_liyapunov1.log", 0) == 0)
        remove("../product/run_liyapunov1.log");
    if (access("../product/position_liyapunov1.txt", 0) == 0)
        remove("../product/position_liyapunov1.txt");
    if(access("../product/run.log",0)==0)
    rename("../product/run.log", "../product/run_liyapunov1.log");
    else 
    {
        cout<<"Liyapunov Test Error: Can't find run.log"<<endl;
        return;
    }
    if(access("../product/position.txt",0)==0)
    rename("../product/position.txt", "../product/position_liyapunov1.txt");
    else
    {
        cout<<"Liyapunov Test Error: Can't find position.txt"<<endl;
        return;
    }
    if (initial_input() && geo_input())
    {
        p.x[1] += 0.0000000001; // 0.0000000001 case
        run_verlet_velocity(0, true, true);
        Liyapunov_Position_Compare();
    }
    return;
}