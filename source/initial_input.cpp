#include "../include/GeneralConstruction.h"
#include "../include/initial_input.h"
int initial_input()
{

    p.fin.open("../initial_input/INPUT", ios::in);
    if (!p.fin.is_open())
    {
        cout << "Failed to Read File INPUT" << endl;
        p.fin.close();
        return 0;
    }
    else
    {
        while (!(p.fin).eof())
        {
            p.fin >> p.filenotes;
            if (p.filenotes.compare("id") == 0)
            {
                p.fin >> p.id;
            }
            if (p.filenotes.compare("natoms") == 0)
            {
                p.fin >> p.atomsum;
            }
            if (p.filenotes.compare("r_cut") == 0)
            {
                p.fin >> p.r_cut;
            }
            if (p.filenotes.compare("neighbor_r") == 0)
            {
                p.fin >> p.neighbor_r;
            }
            if (p.filenotes.compare("geo_dir") == 0)
            {
                p.fin >> p.geo_dir;
            }
            if (p.filenotes.compare("neighbor_n") == 0)
            {
                p.fin >> p.neighbor_n;
            }
            if (p.filenotes.compare("epsilon") == 0)
            {
                p.fin >> p.epsilon;
            }
            if (p.filenotes.compare("sigma") == 0)
            {
                p.fin >> p.sigma;
            }
            if (p.filenotes.compare("nstep_output") == 0)
            {
                p.fin >> p.nstep_output;
            }
            if (p.filenotes.compare("nstep_search") == 0)
            {
                p.fin >> p.nstep_search;
            }
            if (p.filenotes.compare("read_v") == 0)
            {
                p.fin >> p.read_v;
            }
            if (p.filenotes.compare("dt") == 0)
            {
                p.fin >> p.dt;
            }
            if (p.filenotes.compare("molar_mass") == 0)
            {
                p.fin >> p.molar_mass;
            }
            if (p.filenotes.compare("nsteps") == 0)
            {
                p.fin >> p.nsteps;
            }
             if (p.filenotes.compare("tau") == 0)
            {
                p.fin >> p.tau;
            }
             if (p.filenotes.compare("nraise") == 0)
            {
                p.fin >> p.nraise;
            }
            if (p.filenotes.compare("thermostat") == 0)
            {
                p.fin >> p.thermostat;
            }
            if (p.filenotes.compare("NVT_Temperature") == 0)
            {
                p.fin >> p.NVT_Temperature;
            }
             if (p.filenotes.compare("initial_temperature") == 0)
            {
                p.fin >> p.initial_temperature;
            }
            if (p.filenotes.compare("cal_rdf") == 0)
            {
                p.fin >> p.cal_rdf;
            }
            if (p.filenotes.compare("rdf_rcut") == 0)
            {
                p.fin >> p.rdf_rcut;
            }
            if (p.filenotes.compare("rdf_dr") == 0)
            {
                p.fin >> p.rdf_dr;
            }
            if (p.filenotes.compare("rdf_start") == 0)
            {
                p.fin >> p.rdf_start;
            }
            if (p.filenotes.compare("rdf_interval") == 0)
            {
                p.fin >> p.rdf_interval;
            }
            if (p.filenotes.compare("rdf_stop") == 0)
            {
                p.fin >> p.rdf_stop;
            }
        }
        cout << "Done: Read in INPUT" << endl;
        p.fin.close();
        return 1;
    }
}