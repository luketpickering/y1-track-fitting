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
    
    cout << count << endl;
    
    return 0;
}

