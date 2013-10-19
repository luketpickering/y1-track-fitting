#include <iostream>
#include <fstream>
#include <bitset>
#include <assert.h>
#include <iomanip>

// ---------------------START USING STATEMENTS-----------------------
using std::cout;
using std::endl;
using std::cerr;
using std::ifstream;
using std::ios;
using std::flush;
using std::bitset;
using std::hex;
using std::setw;
// ----------------------END USING STATEMENTS------------------------

// --------------------------START GLOBALS---------------------------
const bool quiet = true;
ifstream ifs;

bitset<sizeof(unsigned char)*8> bsc(unsigned char a){
    return bitset<sizeof(unsigned char)*8>(a);
}

inline unsigned short LRTs( unsigned short a, unsigned short l, unsigned short r){
    return (((unsigned short)(a << l)) >> (r+l));
}
//Slow const method
/**/
const unsigned short xl = 13;
const unsigned short xr = 0;
const unsigned short yl = 10;
const unsigned short yr = 3;
const unsigned short tl = 0;
const unsigned short tr = 6;

inline unsigned char sget_x(unsigned short wro){
    return LRTs(wro, xl, xr);
}
inline unsigned char sget_y(unsigned short wro){
    return LRTs(wro,yl,yr);
}
inline unsigned short sget_tdc(unsigned short wro){
    return LRTs(wro,tl,tr);
}
/**/

inline unsigned char get_x(unsigned short wro){
    return LRTs(wro, 13, 0);
}
inline unsigned char get_y(unsigned short wro){
    return LRTs(wro,10,3);
}
inline unsigned short get_tdc(unsigned short wro){
    return LRTs(wro,0,6);
}

inline void get_readout(unsigned short wro, unsigned short *out){
    out[0] = get_x(wro);
    out[1] = get_y(wro);
    out[2] = get_tdc(wro);
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
void test_suite(){

    if (!quiet)
        cout << "-----------Running tests-----------\n\n" << flush;

    assert( sizeof(unsigned short) == 2);
    
    //Corresponds to 52-3-2
    unsigned short test = 0xD1A;
    if (!quiet){
    cout << "Read out data retrieval sanity check" << endl
         << "\tASSERT( " << (int)get_x(test) << " " << (int)get_y(test)
         << " " << (int)get_tdc(test) << " ) == ( 2 3 52 )\n" << endl;
    }
    
    assert( (int)get_x(test) == 2 );
    assert( (int)get_y(test) == 3 );
    assert( (int)get_tdc(test) == 52);
    
    if(!quiet)
        cout << "\n-------------Succeded--------------\n\n" << endl;

}

int main(int argv, char* argc[]){

    if (argv != 2){
        cerr << "Wrong number of args supplied,"
                " please just supply a filename" << endl;
        return -1;
    }
    
    test_suite();
    
    std::ifstream ifs(argc[1], ios::in | ios::binary);
    if( ! ifs.good() ){
        cerr << "File is no good." << endl;
        return -1;
    }    

    unsigned short event_ro[8];
    unsigned short var_ro[3];
    
    ifs.read((char *)event_ro, 16);
    print_ro_header();
    
    for (int i = 0; i < 8; i += 1) {
        get_readout(event_ro[i], var_ro);
        print_ro(var_ro);
    }
    
    return 0;
}

