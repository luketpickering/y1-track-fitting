#include <iomanip>
#include "read.hpp"
#include <cmath>
#include "fit.hpp"

// ---------------------START USING STATEMENTS-----------------------
using std::cout;
using std::cerr;
using std::endl;
using std::bitset;
using std::setw;

// ----------------------END USING STATEMENTS------------------------

// --------------------------START GLOBALS---------------------------
const int twidth = 4;
// ---------------------------END GLOBALS----------------------------

// ---------------------------UTIL FUNCS-----------------------------

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
        
    int num_events = 0;
    int tstepcnt = 0;
    
    int cnt = 0;
    double sum_phi = 0.0;
    int ambig_ev = 0;
    double sum_alph = 0.0;
    double c_phi;
    double c_alph;
    int c_ln;
    int stpcnt;
    int no_min = 0;
    
    while(get_event(event)){
        stpcnt = min_nll_for_ev((const unsigned short**)event,c_phi,c_alph,c_ln);
        
        if(c_ln == -1){
            ambig_ev++;
        }else if (c_ln == -2){
            no_min++;
        }else {
            sum_phi += c_phi;
            sum_alph += c_alph;
            cnt++;
            tstepcnt += stpcnt;
        }
        
    }
    
    cout << flush;
    cerr << sum_phi << " " << cnt << " ambig_steps = " << ambig_ev << " no_min = " << no_min << endl 
         << "avg phi = " << (sum_phi/double(cnt)) << endl;
    cerr << "avg# steps = " << (float(tstepcnt)/float(cnt)) << " avg angl = " 
         << (sum_alph/double(cnt))*(180.0/3.1415) << endl;
    return 0;
}

