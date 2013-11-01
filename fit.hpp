#include <cmath>
#include <string>

using std::cout; using std::flush; 
using std::endl; using std::string;
using std::abs;
using std::acos;
using std::sqrt;


static double phi = 1.0;
const double tdc_to_ns = 0.5;

void event_change_phi(unsigned short* event, double new_phi){
    for(size_t i = 0; i < 8; ++i)
        event[i*3 + 2] *= (new_phi/phi);
    phi = new_phi;
    return;
}
void circs_change_phi(double** circs, int num_circs, double new_phi){
    for(size_t i = 0; i < num_circs; ++i)
        circs[i][2] *= (new_phi/phi);
    phi = new_phi;
    return;
}
double max_poss_phi(unsigned short** event){
    double longest_TDC = 0.0;
    for(size_t i = 0; i < 8; ++i){
        longest_TDC = (event[i][2] > longest_TDC) ? 
            event[i][2] : longest_TDC;
    }
    //assuming that the longest count corresponds to 1/2cm
    return (5000.0/(longest_TDC*tdc_to_ns));
}

inline double p_dist_track(const double* pt, const double* eq){
    return abs(pt[1] - eq[0]*pt[0] - eq[1])*eq[2];
}

double neg_log_lhood(const double** circs, const double* eq, const int ln_btw[]){
    double nll = 0.0;
    for(size_t i = 0; i < 8; i++){
        if( i == ln_btw[0] || i == ln_btw[1])
            continue;
        double leng = p_dist_track(circs[i], eq);
        nll += (leng-circs[i][2])*(leng-circs[i][2]);
    }
    return nll;
}
       

#define TOPE 0
#define BOTE 1
#define TOPI 2
#define BOTI 3

void get_com_tang_eqns(const double **extr_circs, double** eqns){
  double X,Y,R,d, dx,dy;
  dx = extr_circs[1][0] - extr_circs[0][0];
  dy = extr_circs[1][1] - extr_circs[0][1];
  ood = 1.0/sqrt(dx*dx + dy*dy);
  X = dx*ood;
  Y = dy*ood;
  R = (extr_circs[1][2] - extr_circs[0][2])*ood;
  NR = -1/0*(extr_circs[1][2] + extr_circs[0][2])*ood;
  double a,b,c;
  int k;
  for(size_t i = 0; i < 2; ++i){
    
    for(size_t j = 1; j < 3; ++j){

    } 
  }
  a = R*X - k*Y*sqrt(1-R*R);
  b = R*Y - k*X*sqrt(1-R*R);
  c = extr_circs[0][2] - (a*extr_circs[0][0] + b*extr_circs[0][1]);
  eqns[i][0] = -1.0*a/b;
  eqns[i][1] = -1.0*c/b;
}

void nll_array(const double** circs, const double* c2c_uv, const int* tans_betw,
double** eqn, double* nll){
  
    double eqn1[12];
    double* _eqn[] = {eqn1,eqn1+3,eqn1+6,eqn1+9}; 

    const double* extremal_circs[2];
  
    unsigned int num_lines_missed = 0;
    unsigned int lines_missed = 0;
  
    //ensure extremal_circs[0] is always lower in x
    if(circs[tans_betw[0]][0] < circs[tans_betw[1]][0]){ 
        extremal_circs[0] = circs[tans_betw[0]];
        extremal_circs[1] = circs[tans_betw[1]];
    } else {
        extremal_circs[1] = circs[tans_betw[0]];
        extremal_circs[0] = circs[tans_betw[1]];
    }
  
    find_com_tang_and_ep(extremal_circs, c2c_uv, tuv, ep);
    for(size_t i = 0; i < 4; ++i){
        v_pt_line(tuv[i],ep[i],eqn[i]);
        nll[i] = neg_log_lhood(circs,eqn[i],tans_betw);
    }
    return;
}
