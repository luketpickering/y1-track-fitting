#include <cmath>
#include <string>
#include <iomanip>
#include <limits>
#include <math.h>
using std::numeric_limits;

using std::cout; using std::flush; 
using std::endl; using std::string;
using std::abs; using std::setw;
using std::acos; using std::ceil;
using std::sqrt;


const double tdc_to_ns = 0.5;

void print_soln(const double** circs, const double** eqns, int eqn_to_print){
    cout << "Proffered Soln" << endl << endl;
    
    for(size_t i = 0; i < 8; ++i)
    {
        cout << setw(4) << circs[i][0] << " " << setw(4) << circs[i][1]
            << " " << setw(4) << circs[i][2] << endl;
    }
    
    for(size_t i = 0; i < 4; ++i)
    {
        if(eqn_to_print & (1 << i)){
            cout << setw(4) << eqns[i][0] << " " << setw(4) << eqns[i][1] << endl;
        }
    }
    
    cout << endl << endl;
}

inline void circs_change_phi(double** circs, int num_circs, double new_phi, double old_phi){
    for(size_t i = 0; i < num_circs; ++i)
        circs[i][2] *= (new_phi/old_phi);
    return;
}
inline double max_poss_phi(const unsigned short** event){
    double longest_TDC = 0.0;
    for(size_t i = 0; i < 8; ++i){
        longest_TDC = (event[i][2] > longest_TDC) ? 
            event[i][2] : longest_TDC;
    }
    //assuming that the longest count corresponds to 1/2cm
    return (5000.0/(longest_TDC*tdc_to_ns));
}

inline double p_dist_track(const double* pt, const double* eq){
    return abs(pt[1] - eq[0]*pt[0] - eq[1])/eq[2];
}

inline double neg_log_lhood(const double** circs, const double* eq, const int ln_btw[]){
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
inline void calc_com_tan_precalcs(const double **extr_circs, double* opc){
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
inline void get_com_tang_eqns(const double **extr_circs, const double* pre_calc, double* out_eqn, int ln_tc){

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


inline void nll_array(const double** circs, const int* tans_betw,
double** eqn, double* nll){

    const double* extremal_circs[2];
    double pre_calcs[3];
  
    //ensure extremal_circs[0] is always lower in x
    if(circs[tans_betw[0]][0] < circs[tans_betw[1]][0]){ 
        extremal_circs[0] = circs[tans_betw[0]];
        extremal_circs[1] = circs[tans_betw[1]];
    } else {
        extremal_circs[1] = circs[tans_betw[0]];
        extremal_circs[0] = circs[tans_betw[1]];
    }
  
    calc_com_tan_precalcs(extremal_circs, pre_calcs);
    for(size_t i = 0; i < 4; ++i){
        get_com_tang_eqns(extremal_circs, (const double*)pre_calcs, eqn[i], i);
        //cout << eqn[i][0] << " " << eqn[i][1] << " " << eqn[i][2] << endl;
        nll[i] = neg_log_lhood(circs, eqn[i],tans_betw);
    }
    //print_soln(circs,(const double**)eqn,15);
    //int *b = 0;
    //*b = 1;
    return;
}

inline bool parse_event(const unsigned short** ev, const double& phi_init, double **out_circs){
    int e0y = 0;
    int sumy = 0;
    
    for(size_t i = 0; i < 8; i++)
    {
        out_circs[i][0] = ev[i][0]*10000.0;
        out_circs[i][1] = ev[i][1]*10000.0 + (int(i%2))*5000.0;
        out_circs[i][2] = ev[i][2]*tdc_to_ns*phi_init;
        if (i == 0)
            e0y = ev[i][1]*2;
        sumy += ev[i][1]*2 + (ev[i][0]%2) - e0y;
    }
    if(abs(sumy)>27){
        return false;
    }
    return true;
}

int min_ind(const double* arr, int size){
    int mini;
    for(size_t i = 0; i < size; ++i)
    {
        mini = arr[i] < arr[mini] ? i : mini;
    }
    return mini;
}

int min_nll_for_ev(const unsigned short** ev,const double& init_phi, const double& phi_sig, double& out_min_phi, 
    double& out_min_alph, int& out_min_ln){
        
    double c_array[24];
    double* circs[8] = {c_array, c_array+3, c_array+6, c_array+9, c_array+12, c_array+15, c_array+18, c_array+21};

    double eqn1[12];
    double* eqns[] = {eqn1,eqn1+3,eqn1+6,eqn1+9};
    int lns_to_prop = 15;
    int btwn[2] = {0,7};
    double max_phi = max_poss_phi(ev);
    
    const double phi_resolution = 0.01;
    double step_size = phi_resolution;
    double search_dist = 5.0;
    
    double range[2] = {init_phi - search_dist*phi_sig, init_phi + search_dist*phi_sig};
    range[0] = (range[0] < phi_resolution ) ? phi_resolution : range[0];
    range[1] = (range[1] > max_phi) ? max_phi : range[1];
    
    if(!parse_event(ev, range[0], circs)){
        out_min_ln = -1;
        return 0;
    }
    
    int steps = (range[1] - range[0])/phi_resolution;

    if(steps < 10){
        steps = 10;
        step_size = (range[1] - range[0])/10.0;
        //cerr << "too few steps." << endl;
    }

    /*cout << "step deets: " << phi_sig << " " << init_phi << " " << max_phi << " " << steps 
            << "\n\t\t" << range[0] << "=" << init_phi - search_dist*phi_sig << "<->" << range[1] <<
                    "=" << init_phi + search_dist*phi_sig <<  endl;*/
    
    double max_dbl = numeric_limits<double>::max();
    
    int glob_min_ln = -1;
    double glob_min = max_dbl;
    double glob_min_phi;
    double min_llh[] = {max_dbl,max_dbl,max_dbl,max_dbl};
    double min_phi[4];
    double min_eqn[8];
    double* min_eqns[] = {min_eqn,min_eqn+2,min_eqn+2,min_eqn+2};
    double c_lls[4];
    double c_phi = range[0];
      
    for(size_t i = 0; i < steps; ++i)
    {        
        nll_array((const double**)circs, btwn, eqns, c_lls);
        //cout << c_phi << " " << c_lls[0] << " " << c_lls[1] << " " << c_lls[2] << " " << c_lls[3] << endl;
        for(size_t j = 0; j < 4; ++j)
        {
            if(c_lls[j] < min_llh[j]){
                if(!isnormal(c_lls[j]))
                    cerr << "is not normal - lhood " << c_lls[j]  << endl;
                if(!isnormal(eqns[j][0]))
                    cerr << "is not normal - grad " << eqns[j][0] << endl;
                if(!isnormal(eqns[j][1]))
                    cerr << "is not normal - int " << eqns[j][1] << endl;
                //cout << "cll:" << c_lls[j] << " " << c_phi << endl;
                min_llh[j] = c_lls[j];
                min_phi[j] = c_phi;
                min_eqns[j][0] = eqns[j][0];
                min_eqns[j][1] = eqns[j][1];
            }
            if(c_lls[j] < glob_min){
                glob_min_ln = j;
                glob_min = c_lls[j];
                glob_min_phi = c_phi;
                //cout << glob_min_phi << endl;
            }
        }
        circs_change_phi(circs, 8, c_phi + step_size, c_phi);
        c_phi += step_size;
    }
    if(glob_min_ln == -1){
        cout << "BUGGING OUT:" << range[0] << " " << init_phi << "+/-" << phi_sig << " " << range[1] << " " << max_phi << " steps: " << steps << endl;
        int *bla = 0;
        *bla = 1;
    }
    
    out_min_ln = glob_min_ln;
    out_min_phi = glob_min_phi;
    out_min_alph = atan(min_eqns[glob_min_ln][0]);
    return steps;  
    
}
