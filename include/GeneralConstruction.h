#ifndef _INIT_
#define _INIT_

// Pure C Libs
/*
extern"C"
{
    #include<stdio.h>
    #include<stdlib.h>
    #include<string.h>
    #include<math.h>
}
*/

// C++ Libs
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <string>
#include <cmath>
#include <algorithm>

#include <random>
#include <ctime>

#include <fstream> //file input output
#include <io.h>
#include <direct.h> //mkdir
#include <chrono>   //clock,c++11
#include<unistd.h>//linux

//Linux: /usr/include/
#include<sys/stat.h>
#include<sys/types.h>
#include<dir.h>


#include <unistd.h> //linux
using namespace std;

const int Neighbor_limit = 400;
const int Atom_limit = 1200;

const double Boltzmann_K = 1.38064852;
const double Avogadro_NA = 6.02214086;
const double ElementaryCharge = 1.602176634;
const double PAI=3.1415926535897932;

class AtomFile
{
public:
    string id;
    string geo_dir;
    int atomsum;
    double r_cut;         //势截断半径
    double neighbor_r;    //额外考虑近邻半径
    double neighbor_n;    //近邻表原子个数上限
    double epsilon;       // LJ势参数
    double sigma;         // LJ势参数
    double x[Atom_limit]; // position
    double y[Atom_limit];
    double z[Atom_limit];
    double u[Atom_limit]; // velocity
    double v[Atom_limit];
    double w[Atom_limit];
    int nlist[Atom_limit];                //第i个原子的近邻原子数
    int list[Atom_limit][Neighbor_limit]; //第i个原子的第j个近邻原子的编号
    // int list_sorted[Atom_limit][Atom_limit];
    // double distance[Atom_limit][Atom_limit];
    // double rx[Atom_limit][Atom_limit];
    // double ry[Atom_limit][Atom_limit];
    // double rz[Atom_limit][Atom_limit];
    double LJ_Pot[Atom_limit]; // LJ-Potential
    double fx[Atom_limit];     // force
    double fy[Atom_limit];
    double fz[Atom_limit];
    double Xborder;
    double Yborder;
    double Zborder;
    int nstep_output;   // for output
    int nstep_search;   // update verlet list
    double read_v;      // 1:read 2:random
    double dt;          // ps
    double ek;          // kinetic energy
    double ep;          // potential energy
    double temperature; // temperature
    double molar_mass;
    double nsteps;          // number of sim steps
    double tau;             // for NVT-Berendson thermostat
    double nraise;          // for NVT-Andersen thermostat
    double NVT_Temperature; // NVT
    double initial_temperature;
    int cal_rdf;
    double rdf_rcut;
    double rdf_dr;
    int rdf_start;
    int rdf_interval;
    int rdf_stop;
    double g_rdf[10100][1000];
    int n[Atom_limit][1000];//second number:RDF shell number
    double gfinal[1000];
    string thermostat;
    ifstream fin; // use it to read atomfile info
    ofstream fout;
    string filenotes; //read files
};
extern class AtomFile p;
#endif