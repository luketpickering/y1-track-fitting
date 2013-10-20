#include <bitset>
#include <iomanip>
#include "read.hpp"

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
    cout << setw(twidth) << wro[0] << setw(twidth) << wro[1]
         << setw(twidth) << wro[2] << endl;
}
inline void print_ro_header(){
    if (quiet)
        return;
    cout << "Event Read Out" << endl
         << setw(twidth) << "X" << setw(twidth) << "Y"
         << setw(twidth) << "TDC" << endl;
}


int main(int argv, char* argc[]){

    bool file_ok = stream_init(argc[1]);
    unsigned short event[24];
    int num_events = 0;
    while (get_event(event)) {
        //for (int i = 0; i < 8; i++) {
          //print_ro( ( event+(3*i) ) );
        //}
        num_events++;
    }
    
    cout << "Found " << num_events << " events." << endl;
    
    return 0;
}

