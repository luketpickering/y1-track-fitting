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
    {
        event[i*3 + 2] *= (new_phi/phi);
    }
    phi = new_phi;
    return;
}
void circs_change_phi(double** circs, int num_circs, double* new_phi){
    for(size_t i = 0; i < num_circs; ++i)
    {
        circs[i][2] *= (*new_phi/phi);
    }
    phi = *new_phi;
    return;
}
double max_poss_phi(unsigned short* event){
    double longest_TDC = 0.0;
    for(size_t i = 0; i < 8; ++i)
    {
        longest_TDC = (event[i*3 + 2] > longest_TDC) ? 
                event[i*3 + 2] : longest_TDC;
    }
    //assuming that the longest count corresponds to 1/2cm
    return (5000.0/(longest_TDC*tdc_to_ns));
}

double swap_cosa_sina(const double &triga){
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
    else {
        if (ANGL_T_COS(angular_type)){
            cosa = angl_info;
            sina = swap_cosa_sina(angl_info);
        }
        else if (ANGL_T_SIN(angular_type)){
            cosa = swap_cosa_sina(angl_info);
            sina = angl_info;
        } else {
            int failboat = 1/0;
        }
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

#define first_sX_sR(a) (a == 0)
#define first_sX_lR(a) (a == 1)
#define first_lX_sR(a) (a == 2)
#define first_lX_lR(a) (a == 3)
char sort_circ_ind(const double** extremal_circs){
    return ((char(extremal_circs[0][0] > extremal_circs[1][0]) << 1) || 
        (char(extremal_circs[0][2] > extremal_circs[1][2]) << 0));
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
}
void find_com_tang_and_ep(double** circs,
double** out_tuv, double** out_ep){
    double* extremal_circs[] = {circs[0],circs[1]};
    char sort_in = sort_circ_ind(circs);
    //now circs[0] is always lower in x
    //bitflip to swap weather circs[0] was remebered as having larger R
    if(sort_in & 2){ 
        sort_in = ~sort_in; 
        circs[0].swap(circs[1]);
    }

    double dxdy[2];
    double len;
    vector<double> c2c_uv = unit_vect_2p(circs[0], circs[1], len, dxdy);
        
    //if R of circs[0] is larger
    if(sort_in & 1){
            
        double sin_ext = sina_ext_ctang(len, circs[0], circs[1]);
        double sin_int = sina_int_ctang(len, circs[0], circs[1]);
        double a_ext = asin(sin_ext);
        double a_int = asin(sin_int);


        //out_tuv[TOPE] = vrot_2D(c2c_uv, -1.0*sin_ext, SIN_ANG_FLAG);
        //out_tuv[BOTE] = vrot_2D(c2c_uv, sin_ext, SIN_ANG_FLAG);
        //out_tuv[TOPI] = vrot_2D(c2c_uv, sin_int, SIN_ANG_FLAG);
        //out_tuv[BOTI] = vrot_2D(c2c_uv, -1.0*sin_int, SIN_ANG_FLAG);

        out_tuv[TOPE] = vrot_2D(c2c_uv, a_ext, 1);
        out_tuv[BOTE] = vrot_2D(c2c_uv, -1.0*a_ext, 1);
        out_tuv[TOPI] = vrot_2D(c2c_uv, -1.0*a_int, 1);
        out_tuv[BOTI] = vrot_2D(c2c_uv, a_int, 1);
            
        //double sin_90p_ext = swap_cosa_sina(sin_ext);
        //double sin_90na_int = swap_cosa_sina(sin_int);
        double a90pa_ext = 90.0*pi/180.0 + a_ext;
        double a90na_int = 90.0*pi/180 - a_int;
            
        out_ep[TOPE] = vadd_2D(vscale_2D(
            vrot_2D(c2c_uv, a90pa_ext, 1), 
        circs[0][2]),circs[0]);
        out_ep[BOTE] = vadd_2D(vscale_2D(
            vrot_2D(c2c_uv, -1.0*a90pa_ext, 1), 
        circs[0][2]),circs[0]);
        out_ep[TOPI] = vadd_2D(vscale_2D(
            vrot_2D(c2c_uv, a90na_int, 1), 
        circs[0][2]),circs[0]);
        out_ep[BOTI] = vadd_2D(vscale_2D(
            vrot_2D(c2c_uv, -1.0*a90na_int, 1), 
        circs[0][2]),circs[0]);

        out_ep[TOPE] = vadd_2D(vscale_2D(
            vrot_2D(c2c_uv, a90pa_ext, 1), 
        circs[0][2]),vector<double>(2,10000000.0));
        out_ep[BOTE] = vscale_2D(
            vrot_2D(c2c_uv, -1.0*a90pa_ext, 1), 
        circs[0][2]);
        out_ep[TOPI] = vscale_2D(
            vrot_2D(c2c_uv, a90na_int, 1), 
        circs[0][2]);
        out_ep[BOTI] = vscale_2D(
            vrot_2D(c2c_uv, -1.0*a90na_int, 1), 
        circs[0][2]);
    } else {
        //cos(angle_at_larger_circle) = sin(angle_at_smaller_circle)
        //want sin to remove ambiguity over rotations about 0
        //double sin_90na_ext = sina_ext_ctang(len, circs[1], circs[0]);
        //double sin_90na_int = sina_int_ctang(len, circs[1], circs[0]);
        double a_ext = acos(sina_ext_ctang(len, circs[1], circs[0]));
        double a_int = acos(sina_int_ctang(len, circs[1], circs[0]));
        double a90_mext = 90.0*pi/180.0 - a_ext;
        double a90_mint = 90.0*pi/180.0 - a_int;

        out_tuv[TOPE] = vrot_2D(c2c_uv, -1.0*a90_mext, 1);
        out_tuv[BOTE] = vrot_2D(c2c_uv, a90_mext, 1);
        out_tuv[TOPI] = vrot_2D(c2c_uv, -1.0*a90_mint, 1);
        out_tuv[BOTI] = vrot_2D(c2c_uv, a90_mint, 1);
            
        //double sin_ext = swap_cosa_sina(sin_ext);
        //double sin_int = swap_cosa_sina(sin_int);

        out_ep[TOPE] = vadd_2D(vscale_2D(
            vrot_2D(c2c_uv, a_ext, 1), 
        circs[0][2]),circs[0]);
        out_ep[BOTE] = vadd_2D(vscale_2D(
            vrot_2D(c2c_uv, -1.0*a_ext, 1), 
        circs[0][2]),circs[0]);
        out_ep[TOPI] = vadd_2D(vscale_2D(
            vrot_2D(c2c_uv, a_int, 1), 
        circs[0][2]),circs[0]);
        out_ep[BOTI] = vadd_2D(vscale_2D(
            vrot_2D(c2c_uv, -1.0*a_int, 1), 
        circs[0][2]),circs[0]);

    }        
}
