#include <bitset>
#include <iomanip>
#include "read.hpp"
#include <cmath>
#include "fit.hpp"

// ---------------------START USING STATEMENTS-----------------------
using std::cout;
using std::endl;
using std::bitset;
using std::setw;
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
    unsigned short wr[24];
    unsigned short* event[8] = {wr, wr+3, wr+6, wr+9, wr+12, wr+15, wr+18, wr+21};
    double c_array[24];
    double* circs[8] = {c_array, c_array+3, c_array+6, c_array+9, c_array+12, c_array+15, 
        c_array+18, c_array+21};
        
    int num_events = 0;
    
    const int max_cout = 1025;
    const int bin_size = 5;
    const int num_bins = max_cout/bin_size;
    
    int hist[num_bins] = {0};
    
    int cnt = 0;
    int stepcnt = 0;
    double sum_phi = 0.0;
    int btwn[2] = {0,7};
    
    double phi_step = 0.5;
    
    while(get_event(event)){
        phi = 1.0;
        for(size_t i = 0; i < 8; i++)
        {
            circs[i][0] = event[i][0]*10000.0;
            circs[i][1] = event[i][1]*10000.0 + (int(i%2))*5000.0;
            circs[i][2] = event[i][2]*1.0;
        
           // cout << all_circs[i][0] << " " << all_circs[i][1] << " " 
             //   << all_circs[i][2] << endl;
        }
    
        double phi_init = max_poss_phi(event) + phi_step;
        circs_change_phi(circs, 8, phi_init);
        double c_phi = phi_init;
        bool found_ln = false;
        int num_lines_miss = 0;
        stepcnt = 0;
        double c2c_uv[3];
        
        unit_vect_2p(circs[btwn[0]],circs[btwn[1]],c2c_uv);
        
        while (c_phi > 0.0 && !found_ln){
            stepcnt++;
            circs_change_phi(circs, 8, c_phi - phi_step);
            c_phi -= phi_step;
            
            num_lines_miss = check_miss((const double**)circs, c2c_uv, btwn);
            //pvect(circs[0]);pvect(circs[1]);
                /*cout << "Phi:" << c_phi << endl 
                    << "Edge Points: " << ep[i][0] << " " << ep[i][1] << endl
                        << "UV = " << tuv[i][0] << " " << tuv[i][1] << endl
                            << "Eqn: y = " << ln_eq[0] << "x + " << ln_eq[1] << endl 
                                << "Line misses all others? " << miss << endl << endl;
                */
            if (num_lines_miss){cnt += 1; sum_phi += c_phi; found_ln = true; break;}
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

