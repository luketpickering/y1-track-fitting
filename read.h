#include <stdio.h>
#include <stdlib.h>
#ifdef WSNOW
#include "snow.h"
#endif

static FILE* f_ptr;

#ifndef BOOLTDEF
#define BOOLTDEF
typedef int bool;
#define false 0;
#define true 1;
#endif

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
#ifdef WSNOW
    init_with_snow();
#endif
    
    f_ptr = fopen(fn, "rb");
    if( ferror(f_ptr) ){
        perror("File is no good:");
    }   
    return ferror(f_ptr);
}
bool get_event(unsigned short** hits){
    unsigned char temp[16];
    
    if(ferror(f_ptr) || feof(f_ptr))
        return false;
    fread((void*)temp,sizeof(unsigned char), 16, f_ptr);
#ifdef WSNOW
    read_with_snow(temp,16,50000);
#endif
    for (char ctr = 0; ctr < 8; ctr += 1)
        parse_event((temp + (ctr*2)), hits[ctr]);
    return (!feof(f_ptr));
}
// --------------------------END READ FUNCS--------------------------

