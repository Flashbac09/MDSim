#include "../include/GeneralConstruction.h"
#include "../include/P_verlet_list.h"

void all_verlet_list_extended()
{

    ofstream fout;
    string folderpath = {"../product/all_verlet_list_extended"};
    /*
    if (access(folderpath.c_str(), 0) != 0) //查看是否存在这个文件夹
    {
        mkdir(folderpath.c_str());
    }
    */
    for (int i = 1; i <= p.atomsum; i++)
    {
        char name[80] = {"../product/all_verlet_list_extended/verlet_list_extended@"};
        fout.open(strcat(name, to_string(i).data())); //拼接起一个路径+名字
        int precision = 13;                           //可以修改小数位或有效数字位
        fout << p.nlist[i] << endl;
        fout << setiosflags(ios::fixed) << setprecision(precision) << p.x[i] << "   " << setiosflags(ios::fixed) << setprecision(precision) << p.y[i] << "   " << setiosflags(ios::fixed) << setprecision(precision) << p.z[i] << "  " << endl;
        for (int j = 1; j <= p.nlist[i]; j++)
        {
            int m = p.list[i][j];
            fout << std::left << setw(8) << m << setiosflags(ios::fixed) << setprecision(precision) << p.x[m]
                 << "    " << setiosflags(ios::fixed) << setprecision(precision) << p.y[m]
                 << "    " << setiosflags(ios::fixed) << setprecision(precision) << p.z[m]
                 << "        " << setiosflags(ios::fixed) << setprecision(precision) << get_distance(i,m) << endl;
        }
        fout.close();
    }
    cout << "Done: all verlet list(extended), " << p.atomsum << " atoms' verlet list in File all_verlet_list_extended." << endl;
}

void all_verlet_list_standard()
{

    ofstream fout;
    string folderpath = {"../product/all_verlet_list_standard"};
    /*
    if (access(folderpath.c_str(), 0) != 0) //查看是否存在这个文件夹
    {
        mkdir(folderpath.c_str());
    }
    */
    for (int i = 1; i <= p.atomsum; i++)
    {
        char name[80] = {"../product/all_verlet_list_standard/verlet_list_standard@"};
        fout.open(strcat(name, to_string(i).data())); //拼接起一个路径+名字
        int precision = 13;                           //可以修改小数位或有效数字位
        // fout<<p.nlist[i]<<endl;
        for (int j = 1; j <= p.nlist[i]; j++)
        {
            int m = p.list[i][j];
            fout << std::left << setw(8) << m << setiosflags(ios::fixed) << setprecision(precision) << p.x[m]
                 << "    " << setiosflags(ios::fixed) << setprecision(precision) << p.y[m]
                 << "    " << setiosflags(ios::fixed) << setprecision(precision) << p.z[m] << endl;
        }
        fout.close();
    }
    cout << "Done: all verlet list(standard), " << p.atomsum << " atoms' verlet list in File all_verlet_list_standard." << endl;
}

void nlist()
{
    ofstream fout;
    fout.open("../product/nlist");
    for (int i = 1; i <= p.atomsum; i++)
    {
        fout << std::left << setw(8) << i << p.nlist[i] << endl;
    }
    fout.close();
    cout << "Done: Ouput nlist" << endl;
}

double get_distance(int i, int j)
{

    double rx = 0, ry = 0, rz = 0, dist = 0;
    rx = p.x[j] - p.x[i];
    while (rx >= p.Xborder / 2)
    {
        rx -= p.Xborder;
    }
    while (rx <= -p.Xborder / 2)
    {
        rx += p.Xborder;
    }
    ry = p.y[j] - p.y[i];
    while (ry >= p.Yborder / 2)
    {
        ry -= p.Yborder;
    }
    while (ry <= -p.Yborder / 2)
    {
        ry += p.Yborder;
    }
    rz = p.z[j] - p.z[i];
    while (rz >= p.Zborder / 2)
    {
        rz -= p.Zborder;
    }
    while (rz <= -p.Zborder / 2)
    {
        rz += p.Zborder;
    }
    dist = sqrt(rx * rx + ry * ry + rz * rz);
    return dist;
}