#include "../include/GeneralConstruction.h"
#include "../include/P_LJ_Potential.h"

void all_LJPotential()
{
    ofstream fout;
    fout.open("../product/energy.txt");
    int precision = 13;
    for (int i = 1; i <= p.atomsum; i++)
    {
        fout << std::left << setw(8) << i << setiosflags(ios::fixed) << setprecision(precision) << p.LJ_Pot[i] << endl;
    }
    fout.close();
    cout << "Done: Output energy.txt" << endl;
    return;
}

void all_LJPforce()
{
    ofstream fout;
    fout.open("../product/force.txt");
    int precision = 13;
    for (int i = 1; i <= p.atomsum; i++)
    {
        fout << std::left << setw(8) << i
             << setiosflags(ios::fixed) << setprecision(precision) << p.fx[i] << "    "
             << setiosflags(ios::fixed) << setprecision(precision) << p.fy[i] << "    "
             << setiosflags(ios::fixed) << setprecision(precision) << p.fz[i] << endl;
    }
    fout.close();
    cout << "Done: Output force.txt" << endl;
}