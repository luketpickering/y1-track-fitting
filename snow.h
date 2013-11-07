#include <sys/ioctl.h>
#include <stdio.h>

#ifndef BOOLTDEF
#define BOOLTDEF
typedef int bool;
#define false 0;
#define true 1;
#endif

static unsigned char* snow;
static int sno_ctr = 0;
static int osc = 0;
static int cols = 80;
static int true_rows = 24;
static int rows;
static int row_len;
static int len;
static int ectr = 0;
static bool filled = false;
static bool do_snow = false;
static FILE* out_f;

static void sno_shift_print(){
    fprintf(out_f,"\n");
    for(int i = rows-1; i >= 0; --i)
        for(size_t j = 0; j < row_len; ++j)
            snow[(i+1)*row_len + j] = snow[i*row_len+j];

    for(size_t i = row_len; i < len; ++i)
    {
        if((i <= 2*row_len-1) || (i > (len-row_len-1))){
            for(size_t j = 0; j < 8; ++j)
                fprintf(out_f, "|");
            continue;
         }
        for(size_t j = 0; j < 8; ++j)
        {
            if(!(i%row_len) && j < 2){
                fprintf(out_f, "|");
                continue;
            }
            else if(!((i+1)%row_len) && j > 5){
                fprintf(out_f, "|");
                continue;
            }
            else if((snow[i] >> j)&1){
                fprintf(out_f, "*");
            }else{
                fprintf(out_f, " ");
            }
        }
    }
    fflush(out_f);
}
void init_with_snow(){
    
    struct winsize w;
    ioctl(0,TIOCGWINSZ, &w);
    
    true_rows = w.ws_row;
    cols = w.ws_col;    
#ifndef SNOW_COUT
    out_f = stderr;
#else
    out_f = stdout;
#endif
        
    rows = true_rows+1;
    len = cols*rows/(sizeof(unsigned char)*8);
    row_len = cols/(sizeof(unsigned char)*8);
    if(cols%(sizeof(unsigned char)*8)){
        fprintf (out_f, "Please resize your terminal to a round multiple of 8 or this will be ugly.\n"
                        "Press [Enter] to continue . . ." );
        fflush ( out_f );
        getchar();
    }
    if(row_len > 0 && len > 0){do_snow = true;}
    snow = malloc(sizeof(unsigned char)*len);
}
void read_with_snow(const unsigned char* dat, int arsize, int freq){
    if(!snow){return;}
    ectr++;
    if(!filled || !(ectr%freq)){
        for(size_t i = osc; i < (arsize+osc) ; ++i){
            snow[sno_ctr] = dat[i%arsize];
            sno_ctr++;

            if(!filled){
                if (sno_ctr == len){
                    filled = true;
                    sno_ctr = 0;
                }else{
                    continue;
                }
            }
            if(sno_ctr == row_len){
                sno_shift_print();
                osc = i % arsize;
                sno_ctr = 0;
                break;
            }
        }
    }
}
void end_snow(){
    free(snow);
}