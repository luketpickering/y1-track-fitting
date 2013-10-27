#include <vector>
#include <cmath>
#include <string>

using std::vector; using std:: cout; using std::flush; 
using std::endl; using std::string;

#define ANGL_T_ANGL(a) (a == (1 << 0))
#define ANGL_T_COS(a) (a == (1 << 1))
#define ANGL_T_SIN(a) (a == (1 << 2))
#define SIN_ANG_FLAG (1<<2)

const double tdc_to_ns = 0.5;
double phi = 38.0;

#ifdef DEBUG
class Circle{
public:
    vector<double> cont;
    vector<double>& operator[](const int ind){
        if (ind > 2){ int failboat = 1/0; }
        if ( cont[ind] < 0.0 ){ int failboat = 1/0; }
        return cont[ind];
    }
};
class Eqn{
public:
    vector<double> cont;
    vector<double>& operator[](const int ind){
        if (ind > 2){ int failboat = 1/0;}
        return cont[ind];
    }
};
#else
typedef vector<double> Circle;
typedef vector<double> Eqn;
#endif

double swap_cosa_sina(const double &triga){
    return sqrt(1.0-triga*triga);
}

vector<double> vadd_2D(const vector<double> &v1, const vector<double> &v2){
    vector<double> rtnv(2,v1[0] + v2[0]);
    rtnv[1] =  v1[1] + v2[0];
    return rtnv;
}
vector<double> vscale_2D(const vector<double> &v1, const double &sf){
    vector<double> rtnv(v1);
    rtnv[0] *= sf;
    rtnv[1] *= sf;
    return rtnv;
}

vector<double> vrot_2D(const vector<double> &v, const double &angl_info, 
    const unsigned char angular_type){
         vector<double> rtnv(2,0.0);
        if (ANGL_T_ANGL(angular_type)){
            rtnv[0] = v[0]*cos(angl_info) - v[1]*sin(angl_info);
            rtnv[1] = v[0]*sin(angl_info) + v[1]*cos(angl_info);
        }
        else {
            double cosa, sina;
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
            rtnv[0] = v[0]*cosa - v[1]*sina;
            rtnv[1] = v[0]*sina + v[1]*cosa;
        }
        return rtnv;
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
    return abs(pt[1] - eq[0]*pt[0] - eq[1])/eq[2];
}

vector<double> unit_vect_2p(const vector<double>& P1, const vector<double>& P2,
double& out_len, double *out_dxdy){
    out_dxdy[0] = P2[0] - P1[0];
    out_dxdy[1] = P2[1] - P1[1];
    out_len = sqrt(out_dxdy[0]*out_dxdy[0] + out_dxdy[1]*out_dxdy[1]);
    vector<double> rtnv(2,out_dxdy[0]/out_len);
    rtnv[1] = out_dxdy[1]/out_len;
    return rtnv;
} 

void micron_grid(vector< vector<double> >& grid){
    for(size_t i = 0; i < 8; i++)
        for(size_t j = 0; j < 8; j++)
            grid[i][j] = j*10000.0 + ( (j&1) ? 5000.0 : 0.0 );  
}

bool line_miss_all_others(const vector< Circle >& circs, const Eqn& eq,
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

Eqn v_pt_line(const vector<double>& V1, const vector<double>& pt){
    Eqn eq(2, V1[1]/V1[0]);
    eq[1] = pt[0] - eq[0]*pt[0];
    return eq;
}

double max_poss_phi(const vector< Circle >& circs){
    double max_tdc = 0;
    for(size_t i = 0; i < 8; ++i){
        if (circs[i][2] > max_tdc)
            max_tdc = circs[i][2];
    }
    //microns in 1/2cm
    return 5000.0/(max_tdc*tdc_to_ns);
}        

#define first_sX_sR(a) (a == 0)
#define first_sX_lR(a) (a == 1)
#define first_lX_sR(a) (a == 2)
#define first_lX_lR(a) (a == 3)
char sort_circ_ind(const vector< Circle >& circs){
    return ((char(circs[0][0] > circs[1][0]) << 1) || 
                (char(circs[0][2] > circs[1][2]) << 0));
}

#define TOPE 0
#define BOTE 1
#define TOPI 2
#define BOTI 3
void pvect(const vector<double>& vect,const string& head="Vect:"){
    cout << head <<  " "  << flush;
    for(size_t i = 0; i < vect.size(); ++i)
    {
        cout << " x_"<<i<<": " << vect[i];
    }
    cout << endl;
}
void find_com_tang_and_ep(vector< Circle >& circs,
    vector< vector<double> >& out_tuv, vector< vector<double> >& out_ep){
        char sort_in = sort_circ_ind(circs);
        pvect(circs[0]);
        pvect(circs[1]);
        //now circs[0] is always lower in x
        //bitflip to swap weather circs[0] was remebered as having larger R
        if(sort_in & 2){ sort_in = ~sort_in; circs[0].swap(circs[1]);}
        pvect(circs[0]);
        pvect(circs[1]);

        double dxdy[2];
        double len;
        vector<double> c2c_uv = unit_vect_2p(circs[0], circs[1], len, dxdy);
        cout << circs[0][0] << " HELLO" << endl;
        
        //if R of circ[1] is larger
        if(sort_in & 1){
            
            double sin_ext = sina_ext_ctang(len, circs[0], circs[1]);
            double sin_int = sina_int_ctang(len, circs[0], circs[1]);
            
            out_tuv[TOPE] = vrot_2D(c2c_uv, -1.0*sin_ext, SIN_ANG_FLAG);
            out_tuv[BOTE] = vrot_2D(c2c_uv, sin_ext, SIN_ANG_FLAG);
            out_tuv[TOPI] = vrot_2D(c2c_uv, sin_int, SIN_ANG_FLAG);
            out_tuv[BOTI] = vrot_2D(c2c_uv, -1.0*sin_int, SIN_ANG_FLAG);
            
            double sin_90p_ext = swap_cosa_sina(sin_ext);
            double sin_90na_int = swap_cosa_sina(sin_int);
            
            out_ep[TOPE] = vscale_2D(
                    vrot_2D(c2c_uv, -1.0*sin_90p_ext, SIN_ANG_FLAG), 
                    circs[0][2]);
            out_ep[BOTE] = vscale_2D(
                    vrot_2D(c2c_uv, sin_90p_ext, SIN_ANG_FLAG), circs[0][2]);
            out_ep[TOPI] = vscale_2D(
                    vrot_2D(c2c_uv, -1.0*sin_90na_int, SIN_ANG_FLAG), 
                    circs[0][2]);
            out_ep[BOTI] = vscale_2D(
                    vrot_2D(c2c_uv, sin_90na_int, SIN_ANG_FLAG), circs[0][2]);
            
        } else {
            //cos(angle_at_larger_circle) = sin(angle_at_smaller_circle)
            //want sin to remove ambiguity over rotations about 0
            double sin_90na_ext = sina_ext_ctang(len, circs[1], circs[0]);
            double sin_90na_int = sina_int_ctang(len, circs[1], circs[0]);
            
            out_tuv[TOPE] = vrot_2D(c2c_uv, sin_90na_ext, SIN_ANG_FLAG);
            out_tuv[BOTE] = vrot_2D(c2c_uv, -1.0*sin_90na_ext, SIN_ANG_FLAG);
            out_tuv[TOPI] = vrot_2D(c2c_uv, -1.0*sin_90na_int, SIN_ANG_FLAG);
            out_tuv[BOTI] = vrot_2D(c2c_uv, sin_90na_int, SIN_ANG_FLAG);
            
            double sin_ext = swap_cosa_sina(sin_ext);
            double sin_int = swap_cosa_sina(sin_int);
            out_ep[TOPE] = vscale_2D(
                    vrot_2D(c2c_uv, -1.0*sin_ext, SIN_ANG_FLAG), circs[0][2]);
            out_ep[BOTE] = vscale_2D(
                    vrot_2D(c2c_uv, sin_ext, SIN_ANG_FLAG), circs[0][2]);
            out_ep[TOPI] = vscale_2D(
                    vrot_2D(c2c_uv, -1.0*sin_int, SIN_ANG_FLAG), circs[0][2]);
            out_ep[BOTI] = vscale_2D(
                    vrot_2D(c2c_uv, sin_int, SIN_ANG_FLAG), circs[0][2]);
        }
        
        
    }