#include <cmath>
#include <string>

using std:: cout; using std::flush; 
using std::endl; using std::string;
using std::abs; using std::cos; using std::sin;
using std::acos; using std::asin;

#define ANGL_T_ANGL(a) (a == (1 << 0))
#define ANGL_T_COS(a) (a == (1 << 1))
#define ANGL_T_SIN(a) (a == (1 << 2))
#define SIN_ANG_FLAG (1<<2)

static const double tdc_to_ns = 1.0;
static double phi = 1.0;
static const double pi = asin(1.0)*2.0;

void event_change_phi(unsigned short* event, double new_phi){
  for(size_t i = 0; i < 8; ++i)
    event[i*3 + 2] *= (new_phi/phi);
  phi = new_phi;
  return;
}
void circs_change_phi(double** circs, int num_circs, double* new_phi){
  for(size_t i = 0; i < num_circs; ++i)
    circs[i][2] *= (*new_phi/phi);
  phi = *new_phi;
  return;
}
double max_poss_phi(unsigned short* event){
  double longest_TDC = 0.0;
  for(size_t i = 0; i < 8; ++i){
    longest_TDC = (event[i*3 + 2] > longest_TDC) ? 
      event[i*3 + 2] : longest_TDC;
  }
  //assuming that the longest count corresponds to 1/2cm
  return (5000.0/(longest_TDC*tdc_to_ns));
}

inline double swap_cosa_sina(const double &triga){
  return sqrt(1.0-triga*triga);
}

void vadd_2D(const double* v1, const double* v2, double* vout){
  vout[0] = v1[0] + v2[0];
  vout[1] = v1[1] + v2[1];
  return;
}
void vscale_2D(const double* v1, const double* sf, double* vout){
  vout[0] = v1[0]*(*sf);
  vout[1] = v1[1]*(*sf);
  return;
}

void vrot_2D(const double* v, const double &angl_info, 
	     const unsigned char angular_type, double* vout){
  vector<double> rtnv(2,0.0);
  double cosa, sina;
  if (ANGL_T_ANGL(angular_type)){
    cosa = cos(angl_info);
    sina = sin(angl_info);
  }
  else if (ANGL_T_COS(angular_type)){
    cosa = angl_info;
    sina = swap_cosa_sina(angl_info);
  }
  else if (ANGL_T_SIN(angular_type)){
    cosa = swap_cosa_sina(angl_info);
    sina = angl_info;
  } else {
    int failboat = 1/0;
    
  }
  vout[0] = v[0]*cosa - v[1]*sina;
  vout[1] = v[0]*sina + v[1]*cosa;
  return;
}

inline double sina_ext_ctang(const double &len, const Circle &C1, const Circle &C2){
#ifdef DEBUG
  if ((C2[2] - C1[2])<0.0){ int failboat = 1/0;}
#endif
  return (C2[2] - C1[2])/len;
}
inline double sina_int_ctang(const double &len, const Circle &C1, const Circle &C2){
  return (C2[2] + C1[2])/len;
}

inline double p_dist_track(const vector<double>& pt, const Eqn& eq){
  return abs(pt[1] - eq[0]*pt[0] - eq[1])/sqrt(1.0+eq[0]*eq[0]);
}

void unit_vect_2p(const vector<double>& P1, const vector<double>& P2,
		  double& out_len, double *out_dxdy, double* vout){
  out_dxdy[0] = P2[0] - P1[0];
  out_dxdy[1] = P2[1] - P1[1];
  out_len = sqrt(out_dxdy[0]*out_dxdy[0] + out_dxdy[1]*out_dxdy[1]);
  vout[0] = (out_dxdy[0]/out_len);
  vout[1] = out_dxdy[1]/out_len;
  return;
} 

bool line_miss_all_others(const double** circs, const double* eq,
			  const int ln_btw[]){
  for(size_t i = 0; i < 8; i++){
    if( i == ln_btw[0] || i == ln_btw[1])
      continue;
    double leng = p_dist_track(circs[i], eq);
    if (leng < circs[i][2])
      return false;
  }
  return true;
}

void v_pt_line(const double* V1, const double* pt, double* eqn_out){
  eqn_out[0] = (V1[1]/V1[0]);
  eqn_out[1] = pt[1] - eq[0]*pt[0];
  eqn_out[2] = sqrt(1.0+eq[0]*eq[0]);
  return;
}       

#define TOPE 0
#define BOTE 1
#define TOPI 2
#define BOTI 3
void pvect(const double* vect, int vsize, const string& head="Vect:"){
  cout << head <<  " "  << flush;
  for(size_t i = 0; i < vsize; ++i)
    {
      cout << " x_"<<i<<": " << vect[i];
    }
  cout << endl;
  return;
}
void pavect(const double** vects, int asize, int vsize){
  for(size_t i = 0; i < asize; i++)
    pvect(vects[i],vsize);
  return;
}
void find_com_tang_and_ep(double** extremal_circs, const double* c2c_uv,
			  double** out_tuv, double** out_ep){

  static double* pert_angls[4];
  static double temp[2], rad_vect[2];
  vscale_2D(c2c_uv, extremal_circs[0][2], rad_vect);

  //Atempt to do it in one
  if(extremal_circs[0][2] < extremal_circ[1][2]){
    double a_ext = asin(sina_ext_ctang(len, extremal_circs[0], extremal_circs[1]));
    double a_int = asin(sina_int_ctang(len, extremal_circs[0], extremal_circs[1]));
    double a90pa_ext = 90.0*pi/180.0 + a_ext;
    double a90na_int = 90.0*pi/180 - a_int;

    pert_angls[0] = a_ext; pert_angls[1] = a_int;
    pert_angls[2] = a90pa_ext; pert_angls[3] = a90na_int;

  } else{
    double a_ext = acos(sina_ext_ctang(len, extremal_circs[1], extremal_circs[0]));
    double a_int = acos(sina_int_ctang(len, extremal_circs[1], extremal_circs[0]));
    double a90_mext = 90.0*pi/180.0 - a_ext;
    double a90_mint = 90.0*pi/180.0 - a_int;

    pert_angls[0] = -1.0*a90_mext; pert_angls[1] = -1.0*a90_mint;
    pert_angls[2] = a_ext; pert_angls[3] = a_int;
  }

  for(size_t i = 0; i < 2; i++){
    vrot_2D(c2c_uv, pert_angls[i], 1, out_tuv[i*2]);
    vrot_2D(c2c_uv, -1.0*pert_angls[i], 1, out_tuv[i*2 + 1]);
	    
    vrot_2D(rad_vect, pert_angls[i+2], 1, temp);
    vadd_2D(temp, extremal_circs[0], out_ep[i*2]);

    vrot_2D(rad_vect, -1.0*pert_angls[i+2], 1, temp);
    vadd_2D(temp, extremal_circs[0],out_ep[i*2 + 1]);
  } 

  return;
}

bool check_miss(const double** circs, const double* c2c_uv, const int* tans_betw){
  static double epTE[2],epBE[2],epTI[2],epBI[2];
  static double* ep[] = {epTE,epBE,epTI,epBI};
  static double tuvTE[2],tuvBE[2],tuvTI[2],tuvBI[2];
  static double* tuv[] = {tuvTE,tuvBE,tuvTI,tuvBI};

  static double* extremal_circs[2];

  //ensure extremal_circs[0] is always lower in x
  if(circs[tans_betw[0]][0] < extremal_circs[tans_betw[1]][0]){ 
    extremal_circs[0] = circs[tans_betw[0]];
    extremal_circs[1] = circs[tans_betw[1]];
  } else {
    extremal_circs[1] = circs[tans_betw[0]];
    extremal_circs[0] = circs[tans_betw[1]];
  }
}
    
    /*
    vrot_2D(c2c_uv, pert_angls[0], 1, out_tuv[TOPE]);
    vrot_2D(c2c_uv, -1.0*pert_angls[0], 1, out_tuv[BOTE]);
    vrot_2D(c2c_uv, pert_angls[1], 1, out_tuv[TOPI]);
    vrot_2D(c2c_uv, -1.0*pert_angls[1], 1, out_tuv[BOTI]);
    vrot_2D(rad_vect, pert_angls[2], 1, temp);
    vadd_2D(temp, circs[0], out_ep[TOPE]);

    vrot_2D(rad_vect, -1.0*pert_angls[2], 1, temp);
    vadd_2D(temp, circs[0],out_ep[BOTE]);

    vrot_2D(rad_vect, pert_angls[3], 1, temp);
    vadd_2D(temp,circs[0],out_ep[TOPI])

    vrot_2D(rad_vect, -1.0*pert_angls[3], 1, temp);
    vadd_2D(temp,circs[0],out_ep[BOTI])



    if(sort_in & 1){
        double a_ext = asin(sina_ext_ctang(len, circs[0], circs[1]));
	double a_int = asin(sina_int_ctang(len, circs[0], circs[1]));

        double a90pa_ext = 90.0*pi/180.0 + a_ext;
        double a90na_int = 90.0*pi/180 - a_int;

#define ALPH_EXT 0
#define ALPH_INT 1
#define _90_PL_A_EXT 2
#define _90_MIN_A_INT 3
	pert_angls[ALPH_EXT] = a_ext; pert_angls[ALPH_INT] = a_int;
	pert_angls[_90_PL_A_EXT] = a90pa_ext; pert_angls[_90_MIN_A_INT] = a90na_int;

    } else{
        double a_ext = acos(sina_ext_ctang(len, circs[1], circs[0]));
        double a_int = acos(sina_int_ctang(len, circs[1], circs[0]));
        double a90_mext = 90.0*pi/180.0 - a_ext;
        double a90_mint = 90.0*pi/180.0 - a_int;

#define ALPH_EXT 0
#define ALPH_INT 1
#define _90_MIN_A_EXT 2
#define _90_MIN_A_INT 3
	pert_angls[ALPH_EXT] = a_ext; pert_angls[ALPH_INT] = a_int;
	pert_angls[_90_MIN_A_EXT] = a90_mext; pert_angls[_90_MIN_A_INT] = a90_mint;
    }
        
    //if R of circs[0] is larger
      if(sort_in & 1){
            


        out_tuv[TOPE] = vrot_2D(c2c_uv, a_ext, 1);
        out_tuv[BOTE] = vrot_2D(c2c_uv, -1.0*a_ext, 1);
        out_tuv[TOPI] = vrot_2D(c2c_uv, -1.0*a_int, 1);
        out_tuv[BOTI] = vrot_2D(c2c_uv, a_int, 1);
            
        out_ep[TOPE] = vadd_2D(vscale_2D(vrot_2D(c2c_uv, a90pa_ext, 1), circs[0][2]),circs[0]);
        out_ep[BOTE] = vadd_2D(vscale_2D(vrot_2D(c2c_uv, -1.0*a90pa_ext, 1), circs[0][2]),circs[0]);
        out_ep[TOPI] = vadd_2D(vscale_2D(vrot_2D(c2c_uv, a90na_int, 1), circs[0][2]),circs[0]);
        out_ep[BOTI] = vadd_2D(vscale_2D(vrot_2D(c2c_uv, -1.0*a90na_int, 1), circs[0][2]),circs[0]);
    } else {


        out_tuv[TOPE] = vrot_2D(c2c_uv, -1.0*a90_mext, 1);
        out_tuv[BOTE] = vrot_2D(c2c_uv, a90_mext, 1);
        out_tuv[TOPI] = vrot_2D(c2c_uv, -1.0*a90_mint, 1);
        out_tuv[BOTI] = vrot_2D(c2c_uv, a90_mint, 1);

        out_ep[TOPE] = vadd_2D(vscale_2D(vrot_2D(c2c_uv, a_ext, 1), circs[0][2]),circs[0]);
        out_ep[BOTE] = vadd_2D(vscale_2D(vrot_2D(c2c_uv, -1.0*a_ext, 1), circs[0][2]),circs[0]);
        out_ep[TOPI] = vadd_2D(vscale_2D( vrot_2D(c2c_uv, a_int, 1), circs[0][2]),circs[0]);
        out_ep[BOTI] = vadd_2D(vscale_2D(vrot_2D(c2c_uv, -1.0*a_int, 1), circs[0][2]),circs[0]);
	} */
