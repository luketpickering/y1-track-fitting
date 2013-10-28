#include <iostream>
#include <fstream>

// ---------------------START USING STATEMENTS-----------------------
using std::endl;
using std::cerr;
using std::ifstream;
using std::ios;
// ----------------------END USING STATEMENTS------------------------

// --------------------------START GLOBALS---------------------------
static ifstream ifs;
// ---------------------------END GLOBALS----------------------------

// ------------------------START PARSE FUNCS-------------------------
inline static void get_readout(unsigned short wro, unsigned short *out){
	out[0] = wro & 7;
	out[1] = (wro >> 3) & 7;
	out[2] = (wro >> 6) & 1023;
	return;
}

inline static void parse_event(unsigned char* data, unsigned short* hit){
    unsigned short holder;
    holder = (data[0] << 0) | (data[1] << 8);
    get_readout(holder, hit);
    return;
}
// -------------------------END PARSE FUNCS--------------------------

// -------------------------START READ FUNCS-------------------------
bool stream_init(char *fn){
    ifs.open(fn, ios::in | ios::binary);
    if( ! ifs.good() ){
        cerr << "File is no good." << endl;
    }   
    return ifs.good();
}

bool get_event(unsigned short** hits){
    unsigned char temp[16];
    if(!ifs.good() || ifs.eof())
        return false;
    ifs.read((char*)temp , 16);
    for (char ctr = 0; ctr < 8; ctr += 1)
        parse_event((temp + (ctr*2)), hits[ctr]);
    return !ifs.eof();
}
// --------------------------END READ FUNCS--------------------------

