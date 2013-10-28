#include <bitset>
#include <iomanip>
#include "read.hpp"
#include <cmath>
#include "fit.hpp"

// ---------------------START USING STATEMENTS-----------------------
using std::cout;
using std::endl;
using std::bitset;
// ----------------------END USING STATEMENTS------------------------

// --------------------------START GLOBALS---------------------------
const bool quiet = true;

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
    unsigned short event[24];
    int num_events = 0;
    
    const int max_cout = 1025;
    const int bin_size = 5;
    const int num_bins = max_cout/bin_size;
    
    int hist[num_bins] = {0};
    
    
    vector< vector<double> > circs(2);
    vector< vector<double> > all_circs(8,vector<double>(3,0.0));
    
    vector< vector<double> > tuv(4,vector<double>(2,0.0));
    vector< vector<double> > ep(4,vector<double>(2,0.0));
    
    int cnt = 0;
    int stepcnt = 0;
    double sum_phi = 0.0;
    int btwn[2] = {0,7};
    while(get_event(event)){
        phi = 1.0;
        for(size_t i = 0; i < 8; i++)
        {
            all_circs[i][0] = event[i*3]*10000.0;
            all_circs[i][1] = event[i*3 + 1]*10000.0 + (int(i%2))*5000.0;
            all_circs[i][2] = event[i*3 + 2];
        
           // cout << all_circs[i][0] << " " << all_circs[i][1] << " " 
             //   << all_circs[i][2] << endl;
        }
    
        double phi_init = max_poss_phi(event)+0.5;
    
        double c_phi = circs_change_phi(all_circs,phi_init);
        bool found_ln = false;
        stepcnt = 0;
        while (c_phi > 0.0 && !found_ln){
            stepcnt++;
            c_phi = circs_change_phi(all_circs,c_phi-0.5);
            circs[0] = all_circs[0];
            circs[1] = all_circs[7];
            find_com_tang_and_ep(circs,tuv, ep);
            for(size_t i = 0; i < 4; ++i)
            {
                Eqn ln_eq = v_pt_line(tuv[i],ep[i]);
                bool miss = line_miss_all_others(all_circs, ln_eq, btwn);
                //pvect(circs[0]);pvect(circs[1]);
                /*cout << "Phi:" << c_phi << endl 
                    << "Edge Points: " << ep[i][0] << " " << ep[i][1] << endl
                        << "UV = " << tuv[i][0] << " " << tuv[i][1] << endl
                            << "Eqn: y = " << ln_eq[0] << "x + " << ln_eq[1] << endl 
                                << "Line misses all others? " << miss << endl << endl;
                */if (miss){cnt += 1; sum_phi += c_phi; found_ln = true; break;}
            }
        }
        if (!(cnt%10000)) { 
            cout << "Read event:" << cnt+1 << endl
                 << "\tsteps required:" << stepcnt
                 << "\tphi_init = " << phi_init
                 << "\tphi calc = " << c_phi << endl << endl;
	    return 0;
        }
    }
    cout << sum_phi << " " << cnt << " avg phi = " << (sum_phi/double(cnt)) << endl;
    cout << "avg# steps = " << (stepcnt/cnt) << endl;
#ifdef HIST    
    while (get_event(event)) {
#ifdef VERBOSE
        if (! (num_events % 100000) ) {
            cout << "EVENT: " << num_events << endl;
            for (int i = 0; i < 8; i++) {
                print_ro( ( event+(3*i) ) );
            }
        }
#endif
        for (int i = 0; i < 8; i++) {
            unsigned short tdc = (event+(3*i))[2];
            int bin = floor(tdc/bin_size);
            if ( (bin <0) || (bin > (num_bins - 1)) ){
                cout << bin << endl;
                return -1;
            }
            hist[bin] += 1;
        }
        num_events++;
    }
    
    cout << "Found " << num_events << " events." << endl;
    
    int count = 0;
    for (int i = 0; i < num_bins; i++){
        cout << i << ", " << hist[i] << endl;
        count += hist[i];
        
    }
#endif        
    return 0;
}

