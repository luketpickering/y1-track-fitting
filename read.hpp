#include <iostream>
#include <fstream>
#include <bitset>
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
ifstream ifs;

bitset<sizeof(unsigned short)*8> bss(unsigned short *a){
    return bitset<sizeof(unsigned short)*8>(*a);
}

inline unsigned short LRTs( unsigned short a, unsigned short l, unsigned short r){
    return (((unsigned short)(a << l)) >> (r+l));
}

inline unsigned char get_x(unsigned short wro){
    return LRTs(wro, 13, 0);
}
inline unsigned char get_y(unsigned short wro){
    return LRTs(wro,10,3);
}
inline unsigned short get_tdc(unsigned short wro){
    return LRTs(wro,0,6);
}

inline void get_readout_expl(unsigned short wro, unsigned short *out){
    out[0] = get_x(wro);
    out[1] = get_y(wro);
    out[2] = get_tdc(wro);
}
inline void get_readout_fast(unsigned short wro, unsigned short *out){
	out[0] = wro & 7;
	out[1] = (wro >> 3) & 7;
	out[2] = (wro >> 6) & 1023;
}

inline void parse_event(unsigned char* data, unsigned short* hit){
    unsigned short holder;
    holder = (data[0] << 0) | (data[1] << 8);
#ifndef VERBOSE
    get_readout_fast(holder, hit);
#else
	get_readout_expl(holder, hit);
#endif
    return;
}

bool stream_init(char *fn){

    ifs.open(fn, ios::in | ios::binary);
    if( ! ifs.good() ){
        cerr << "File is no good." << endl;
    }
    
    return ifs.good();

}

bool get_event(unsigned short* hits){
    static unsigned char temp[16];
    if(!ifs.good() || ifs.eof())
        return false;
    
    ifs.read((char*)temp , 16);
    
    for (char ctr = 0; ctr < 8; ctr += 1)
        parse_event((temp + (ctr*2)), hits + (ctr*3) );
    return !ifs.eof();
}

