#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "read.h"
#include "fit.h"

#ifndef OPT_BLOB
#define OPT_BLOB
typedef struct _opt{
    double step_size;
    int min_steps;
    char* input_fn;
} opts;
#endif

// ---------------------------UTIL FUNCS-----------------------------

int parse_opts(int argc, char** argv, opts* out_opt_blob){
    int c;
    opterr = 0;
    int reqopts = 0;
    int *bla = 0;
    while((c = getopt(argc, argv, "f:s:m:")) != -1){
        switch(c){
            case 'f':
                asprintf(&(out_opt_blob->input_fn),"%s", optarg);
                reqopts ^= 1;
                break;
            case 's':
                out_opt_blob->step_size = atof(optarg);
                reqopts ^= 2;
                
                break;
            case 'm':
                out_opt_blob->min_steps = atoi(optarg);
                reqopts ^= 4;
                
                break;
            case '?':
                fprintf(stderr, "Invalid Options.\n Require -f \'filename\' -s <step size> "
                                "-m <min number of steps>.\n");
                *bla = 1;
                return -1;
            default:
            *bla = 1;
            
                return -1;       
        }
    }
    if(reqopts != 7){
        fprintf(stderr, "Requires all options.\n Require -f \'filename\' -s <step size> "
                        "-m <min number of steps>.\n");
        *bla = 1;
        
        return -1;
    }
}

int main(int argc, char* argv[]){
    
    opts opt_b;
    int valid_opts = parse_opts(argc, argv, &opt_b);
    
    printf("%i %f %s\n ", opt_b.min_steps, opt_b.step_size, opt_b.input_fn);
    
    if(valid_opts == -1){
        fprintf(stderr, "Invalid Ops. Bugging Out.\n");
        free(opt_b.input_fn);
        stream_close();
        return -1;
    }

    if (!stream_init(opt_b.input_fn)){
        fprintf(stderr, "Failed to read in file.\n");
        return -1;
    }
    
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
    int stpcnt;
    
    //stats
    double running_mean_phi = 52.0;
    double sq_sum_phi = 0.0;
    double running_2nd_mmt = 0.0;
    double stdv_phi = 0.0;
    
    while(get_event(event)){
        if(stdv_phi < 1e-6){
            stdv_phi = 52.0/5.0;
            fprintf(stderr, "used default stdv on step:%d\n",cnt);
        }
        stpcnt = min_nll_for_ev((const unsigned short**)event, running_mean_phi, 
            stdv_phi, &c_phi, &c_alph, &opt_b);
        
        if(stpcnt == -1){
            ambig_ev++;
        }else {
            sum_phi += c_phi; 
            sq_sum_phi += (c_phi*c_phi);  
            sum_alph += c_alph;
            
            cnt++;
            tstepcnt += stpcnt;
            
            running_mean_phi = sum_phi/((double)cnt);
            running_2nd_mmt = sq_sum_phi/((double)cnt);
            stdv_phi = sqrt(running_2nd_mmt - running_mean_phi*running_mean_phi);
            fprintf(stdout, "%f\t%f\n",c_phi, c_alph);
        }        
    }
    fprintf(stderr,"OUT\n\tAvg V_d = %f +/- %f\n\tavg #steps = %f\n\tavg angle = %f\n", 
            running_mean_phi, stdv_phi, ((float)tstepcnt/(float)cnt), 
            (sum_alph/(double)cnt)*(180.0/3.1415));
    stream_close();
    return 0;
}

