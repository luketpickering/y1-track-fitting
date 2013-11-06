#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef OPT_BLOB
#define OPT_BLOB
typedef struct _opt{
    double step_size;
    int min_steps;
    char* input_fn;
} opts;
#endif

const double tdc_to_ns = 0.5;

static inline void circs_change_phi(double** circs, 
                                    double new_phi, 
                                    double old_phi){
    for(size_t i = 0; i < 8; ++i)
        circs[i][2] *= (new_phi/old_phi);
    return;
}

static inline double max_poss_phi(const unsigned short** event){
    double longest_TDC = 0.0;
    for(size_t i = 0; i < 8; ++i){
        longest_TDC = (event[i][2] > longest_TDC) ? 
            event[i][2] : longest_TDC;
    }
    //assuming that the longest count corresponds to 1/2cm
    return (5000.0/(longest_TDC*tdc_to_ns));
}

static inline double p_dist_track(const double* pt, 
                                  const double* eq){
    return abs(pt[1] - eq[0]*pt[0] - eq[1])/eq[2];
}

static inline double neg_log_lhood(const double** circs, 
                                   const double* eq, 
                                   const int ln_btw[]){
    double nll = 0.0;
    for(size_t i = 0; i < 8; i++){
        if( i == ln_btw[0] || i == ln_btw[1])
            continue;
        double leng = p_dist_track(circs[i], eq);
        nll += (leng-circs[i][2])*(leng-circs[i][2]);
    }
    return nll;
}
       
#define PC_ooD 0
#define PC_X 1
#define PC_Y 2
static inline void calc_com_tan_precalcs(const double **extr_circs, 
                                         double* opc){
    double dx,dy;
    dx = extr_circs[1][0] - extr_circs[0][0];
    dy = extr_circs[1][1] - extr_circs[0][1];
    opc[PC_ooD] = 1.0/sqrt(dx*dx + dy*dy);
    opc[PC_X]= dx*opc[PC_ooD];
    opc[PC_Y] = dy*opc[PC_ooD];
}
//k = + 1 bote
//k = - 1 tope
//r2n, k = +1 boti 
//r2n, k = -1 topi
#define BOTE 0
#define BOTI 1
#define TOPE 2
#define TOPI 3
static inline void get_com_tang_eqns(const double **extr_circs, 
                                     const double* pre_calc, 
                                     double* out_eqn, 
                                     int ln_tc){

    double R,a,b,c,a1,a2,b1,b2,rtr2;
    int k,nr;
    nr = (ln_tc&1)? -1.0 : 1.0;
    k = (ln_tc&2)? -1.0 : 1.0;
    R = nr*(extr_circs[1][2] - nr*extr_circs[0][2])*pre_calc[PC_ooD];
    rtr2 = sqrt(1.0-R*R);
    a1 = R*pre_calc[PC_X];
    a2 = pre_calc[PC_Y]*rtr2;
    b1 = R*pre_calc[PC_Y];
    b2 = pre_calc[PC_X]*rtr2;
    a = a1 - k*a2;
    b = b1 + k*b2;
    c = extr_circs[0][2] - (a*extr_circs[0][0] + b*extr_circs[0][1]);
    out_eqn[0] = -1.0*a/b;
    out_eqn[1] = -1.0*c/b;
    out_eqn[2] = sqrt(1.0+out_eqn[0]*out_eqn[0]);
}

static inline void nll_array(const double** circs, 
                             const int* tans_betw,
                             double** eqn, 
                             double* nll, 
                             const double* pre_calcs){

    const double* extremal_circs[2] = {circs[tans_betw[0]], circs[tans_betw[1]]};

    for(size_t i = 0; i < 4; ++i){
        get_com_tang_eqns((const double**) extremal_circs, (const double*)pre_calcs, eqn[i], i);
        nll[i] = neg_log_lhood(circs, eqn[i],tans_betw);
    }
    return;
}

static inline int parse_ev(const unsigned short** ev, 
                           const double phi_init, 
                           double **out_circs){
    int e0y = 0;
    int sumy = 0;
    
    for(size_t i = 0; i < 8; i++)
    {
        out_circs[i][0] = ev[i][0]*10000.0;
        out_circs[i][1] = ev[i][1]*10000.0 + 5000.0 * (double)(i%2) ;
        out_circs[i][2] = ev[i][2]*tdc_to_ns*phi_init;
        if (i == 0)
            e0y = ev[i][1]*2;
        sumy += ev[i][1]*2 + (ev[i][0]%2) - e0y;
    }
    if(abs(sumy)>27){
        return 0;
    }
    return 1;
}
int min_nll_for_ev(const unsigned short** ev, 
                   const double init_phi, 
                   const double phi_sig, 
                   double* out_min_phi, 
                   double* out_min_alph, 
                   opts* opt_b){
        
    double c_array[24];
    double* circs[8] = {c_array, c_array+3, c_array+6, c_array+9, c_array+12, c_array+15, c_array+18, c_array+21};

    double eqn1[12];
    double* eqns[] = {eqn1,eqn1+3,eqn1+6,eqn1+9};
    int lns_to_prop = 15;
    int btwn[2] = {0,7};
    double max_phi = max_poss_phi(ev);
    
    const double phi_resolution = opt_b->step_size;
    int step_dirn = 1;
    double search_dist = 5.0;
    
    double range[2] = {init_phi - search_dist*phi_sig, init_phi + search_dist*phi_sig};
    range[0] = (range[0] < phi_resolution ) ? phi_resolution : range[0];
    range[1] = (range[1] > max_phi) ? max_phi : range[1];
    
    double step_size = step_dirn*phi_resolution;
    double start_phi = (step_dirn == -1) ? range[1] : range[0];
    
    if(!parse_ev(ev, start_phi, circs)){
        return -1;
    }
    
    int steps = (range[1] - range[0])/phi_resolution;

    int min_steps = opt_b->min_steps;
    if(steps < min_steps){
        steps = min_steps;
        step_size = step_dirn*(range[1] - range[0])/((double)min_steps);
    }
    double glob_min = DBL_MAX;
    double glob_min_phi;
    double glob_min_talpha;
    double c_lls[4];
    double c_phi = start_phi;
    
    double* extremal_circs[2] = {circs[btwn[0]],circs[btwn[1]]};
    double pre_calcs[3];
    calc_com_tan_precalcs((const double**) extremal_circs, pre_calcs);
    
    
    for(size_t i = 0; i < steps; ++i){   
        nll_array((const double**)circs, btwn, eqns, c_lls, 
                    (const double*)pre_calcs);
        for(size_t j = 0; j < 4; ++j){
            if(c_lls[j] < glob_min){
                glob_min = c_lls[j];
                glob_min_phi = c_phi;
                glob_min_talpha = eqns[j][0];
            }
        }
        circs_change_phi(circs, c_phi + step_size, c_phi);
        c_phi += step_size;
    }
    *out_min_phi = glob_min_phi;
    *out_min_alph = atan(glob_min_talpha);
    return steps;  
}
