#include <stdio.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#ifdef WSNOW
#include "snow.h"
#endif

static FILE* f_ptr;
static int fd;
static unsigned char* bdata;
static struct stat stats;

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
    
    fd = open(fn, O_RDONLY);
    if( fd == -1 ){
        perror("File is no good:");
        return false;
    }   
    
    fstat(fd,&stats);
    fprintf(stderr, "File size = %lld\n", stats.st_size);
    bdata = mmap(NULL, stats.st_size,PROT_READ,MAP_SHARED,fd,0);
    if(bdata == MAP_FAILED){
        perror("Failed to mmap the file.");
        return false;
    }
    return true;
}
bool stream_close(){
    munmap(bdata, stats.st_size);
    close(fd);
    return 1;
}
static int offset = 0;
bool get_event(unsigned short** hits){
#ifdef WSNOW
    read_with_snow(bdata + offset,16,15000);
#endif
    for (char ctr = 0; ctr < 8; ctr += 1)
        parse_event((bdata + offset + (ctr*2)), hits[ctr]);
    offset += 16;
    return (offset != stats.st_size); //(!feof(f_ptr));
}
// --------------------------END READ FUNCS--------------------------

