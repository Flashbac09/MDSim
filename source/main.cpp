// General:
#include "../include/GeneralConstruction.h"
// Input:
#include "../include/initial_input.h"
#include "../include/geo_input.h"
// Work:
#include "../include/verlet_list.h"
#include "../include/verlet_velocity.h"
#include "../include/Liyapunov_Position_Compare.h"
// Output:
#include "../include/P_verlet_list.h"
#include "../include/P_LJ_Potential.h"
#include "../include/P_verlet_velocity.h"
// Main:
#include "../include/main.h"
AtomFile p;

int main()
{

    if (initial_input() && geo_input())
    {
        if(p.read_v==2)random_speed();
        run_verlet_velocity(0, true,false);
        if(p.cal_rdf==1)RDF_out();
        // Ensembles and thermostats(if there exists) are controlled by INPUT file.
        // 3 Parameters: time,process_notice,MotionOut
    }
  // Liyapunov_Test();

    // system("pause");
    return 0;
}
