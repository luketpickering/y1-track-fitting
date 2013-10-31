#include <bitset>
#include <iomanip>
#include "read.hpp"
#include <cmath>
#include "fit.hpp"
#include <limits>

// ---------------------START USING STATEMENTS-----------------------
using std::cout;
using std::cerr;
using std::endl;
using std::bitset;
using std::setw;
using std::numeric_limits;
// ----------------------END USING STATEMENTS------------------------

// --------------------------START GLOBALS---------------------------
const bool quiet = true;

// ---------------------------END GLOBALS----------------------------

// ---------------------------UTIL FUNCS-----------------------------

void print_soln(double** circs, double** eqns, int eqn_to_print){
    cout << "Proffered Soln" << endl << endl;
    
    for(size_t i = 0; i < 8; ++i)
    {
        cout << setw(4) << circs[i][0] << " " << setw(4) << circs[i][1]
             << " " << setw(4) << circs[i][2] << endl;
    }
    
    for(size_t i = 0; i < 4; ++i)
    {
        if(eqn_to_print & (1 << i)){
            cout << setw(4) << eqns[i][0] << " " << setw(4) << eqns[i][1] << endl;
        }
    }
    
    cout << endl << endl;
}

bitset<sizeof(unsigned char)*8> bsc(unsigned char a){
    return bitset<sizeof(unsigned char)*8>(a);
}
const int twidth = 4;
inline void print_ro(unsigned short *wro){
#ifndef VERBOSE
    return;
#endif
    cout << setw(twidth) << wro[0] << setw(twidth) << wro[1]
         << setw(twidth) << wro[2] << endl;
}
inline void print_ro_header(){
#ifndef VERBOSE
    return;
#endif
    cout << "Event Read Out" << endl
         << setw(twidth) << "X" << setw(twidth) << "Y"
         << setw(twidth) << "TDC" << endl;
}

int main(int argv, char* argc[]){

    bool file_ok = stream_init(argc[1]);
    unsigned short wr[24];
    unsigned short* event[8] = {wr, wr+3, wr+6, wr+9, wr+12, wr+15, wr+18, wr+21};
    double c_array[24];
    double* circs[8] = {c_array, c_array+3, c_array+6, c_array+9, c_array+12, c_array+15, 
        c_array+18, c_array+21};
    double eqn1[12];
    double* eqns[] = {eqn1,eqn1+3,eqn1+6,eqn1+9}; 
        
    int num_events = 0;
    int tstepcnt = 0;
    
    int cnt = 0;
    int stepcnt = 0;
    double sum_phi = 0.0;
    int btwn[2] = {0,7};
    int ambig_ev = 0;
    double phi_step = 1.0;
    double sum_alph = 0.0;

    
    while(get_event(event)){
        phi = 1.0;
        double e0y = 0.0;
        double sumy = 0.0;
        for(size_t i = 0; i < 8; i++)
        {
            circs[i][0] = event[i][0]*10000.0;
            circs[i][1] = event[i][1]*10000.0 + (int(i%2))*5000.0;
            circs[i][2] = event[i][2]*tdc_to_ns;
            if (i == 0)
                e0y = circs[i][1];
            sumy += circs[i][1] - e0y;
        }
        if(abs(sumy)>27.0*5000.0){
            ambig_ev++;
            continue;
        }
    
        double phi_init = max_poss_phi(event);
        phi_step = 0.1*phi_init;
        circs_change_phi(circs, 8, phi_init);
        double c_phi = phi_init;
        bool found_ln = false;
        int num_lines_miss = 0;
        int llm = 0;
        stepcnt = 0;
        double c2c_uv[3];
        double nll[4];
        double mnll[4] = {numeric_limits<float>::max()};
        double glob_minll = numeric_limits<float>::max();
        bool glob_minll_changed;
        int glob_min_ln;
        
        
        unit_vect_2p(circs[btwn[0]],circs[btwn[1]],c2c_uv);
        
        while (c_phi > 0.0 && !found_ln){
            stepcnt++;
            circs_change_phi(circs, 8, c_phi - phi_step);
            c_phi -= phi_step;
            glob_minll_changed = false;
            
            nll_array((const double**)circs, c2c_uv, btwn, eqns, nll);
            for(size_t i = 0; i < 4; ++i)
            {
                if(nll[i] < mnll[i]){
                    nll[i] = mnll[i];
                }
                if(nll[i] < glob_minll){
                    glob_minll = nll[i];
                    glob_min_ln = i;
                    glob_minll_changed = true;
                }
            }

            if(!glob_minll_changed){
                found_ln = true;
                cnt += 1; 
                sum_phi += c_phi; found_ln = true; tstepcnt += stepcnt;
                sum_alph += atan(eqns[glob_min_ln][0]);
                
                //cout << glob_min_ln <<  endl;
                cout << atan(eqns[glob_min_ln][0])*(180.0/3.1415) << "\n";
            }
        }

    }
    
    cout << flush;
    cerr << sum_phi << " " << cnt << " avg phi = " << (sum_phi/double(cnt)) << endl;
    cerr << "avg# steps = " << (float(tstepcnt)/float(cnt)) << " avg angl = " << (sum_alph/double(cnt))*(180.0/3.1415) << endl;
    return 0;
}

