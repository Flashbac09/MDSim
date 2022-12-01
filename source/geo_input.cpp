#include "../include/GeneralConstruction.h"
#include "../include/geo_input.h"

int geo_input()
{
    p.fin.open(p.geo_dir, ios::in);
    if (!p.fin.is_open())
    {
        cout << "Failed to Read in File geo.in." << endl;
        p.fin.close();
        return 0;
    }
    else
    {
        int t = 0; // temp calc atomsum
        int posimark = 0;
        double temp1, temp2;
        while (getline(p.fin, p.filenotes)) // Read cell parameters
        {
            if ((p.filenotes).compare("%CELL_PARAMETER") == 0)
            {
                p.fin >> p.Xborder >> temp1 >> temp2; //目前的距离计算方式只对长方体晶胞有效，对平行六面体无效.
                p.fin >> temp1 >> p.Yborder >> temp2;
                p.fin >> temp1 >> temp2 >> p.Zborder;
                break;
            }
        }
        while (getline(p.fin, p.filenotes))
        {
            if (p.filenotes.compare("%ATOMIC_POSTION") == 0)
            {
                break;
            }
            // CONFUSED about this word : POSTION
        }
        while (p.fin >> p.filenotes) // Read Position
        {
            if (strcmp((p.filenotes).data(), "%ATOMIC_VELOCITY") == 0) //只是与上面不同的方法，关于basic string的各种成员函数
            {
                posimark = 1;
                cout << "Done: Read in Position from " << p.geo_dir << endl;
                break;
            }
            t += 1;
            p.fin >> p.x[t];
            p.fin >> p.y[t];
            p.fin >> p.z[t];
            if(t>=p.atomsum&&p.geo_dir.compare("../initial_input/geo2.in")==0)return 1;
        }
        t = 0;
        if (posimark == 1)
        {
            while (t <= p.atomsum)
            {
                t += 1;
                p.fin >> p.filenotes;
                p.fin >> p.u[t];
                p.fin >> p.v[t];
                p.fin >> p.w[t];
            }
            cout << "Done: Read in Velocity from " << p.geo_dir << endl;
        }
        else
            return 0;
        p.fin.close();
        return 1;
    }
}