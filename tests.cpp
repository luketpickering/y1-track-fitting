#include <assert.h>
#include <iostream>
#include "read.hpp"
#include "fit.hpp"
#include <cmath>

using std::cout;
using std::endl;
using std::abs;
const double pi = 2.0*asin(1.0);

bool veq(const vector<double>& a, const vector<double>& b){    
    if( a.size() != b.size() ) {return false;}
    for(size_t i = 0; i < a.size(); ++i)
    {
        if (abs(a[i] - b[i]) > 0.000001) {
            return false;
        }
    }
    return true;
}

int main(){

    cout << "\n-----------Running tests-----------\n\n" << flush;

    assert( sizeof(unsigned short) == 2);
    assert( sizeof(unsigned char) == 1);
    
    //Corresponds to 52-3-2
    unsigned short test = 0xD1A;
    
    cout << "Read out data retrieval sanity check" << endl
         << "\tASSERT( " << (int)get_x(test) << " " << (int)get_y(test)
         << " " << (int)get_tdc(test) << " ) == ( 2 3 52 )" << endl;
    
    assert( (int)get_x(test) == 2 );
    assert( (int)get_y(test) == 3 );
    assert( (int)get_tdc(test) == 52);
    
    unsigned short test_out_expl[3];
    unsigned short test_out_fast[3];

    get_readout_fast(test,test_out_fast);
    get_readout_expl(test,test_out_expl);
    
    cout << "Fast Data parse test" << endl
         << "\tASSERT( " << (int)test_out_expl[0] << " "
         << (int)test_out_expl[1] << " " << (int)test_out_expl[2]
         << " ) == ( " << (int)test_out_fast[0] << " "
         << (int)test_out_fast[1] << " " << (int)test_out_fast[2]
         << " )" << endl;
    
    assert( test_out_expl[0] == test_out_fast[0]);
    assert( test_out_expl[1] == test_out_fast[1]);
    assert( test_out_expl[2] == test_out_fast[2]);

    vector<double> up(2,0.0);
    up[1] = 1.0;
    vector<double> left(2,-1.0);
    left[1] = 0.0;
    vector<double> down(2,0.0);
    down[1] = -1.0;

    vector<double> lfup = vrot_2D(up,pi*90.0/180.0,1);
    vector<double> dfup = vrot_2D(up,pi,1);
    vector<double> lfdwn = vrot_2D(down, pi*270.0/180.0,1);
    vector<double> lfdwn2 = vrot_2D(down, -90.0*pi/180.0,1);
    
    cout << "Vector Tests" << endl
         << "--Rotations by angles" << endl;
    cout << "Pi: " << pi << endl;
    pvect(left, "Left"); cout << "eq" << endl; pvect(lfup,"Rot(90)Up"); 
    cout << endl;
    pvect(down, "Down"); cout << "eq" << endl; pvect(dfup, "Rot(180)Down");
    cout << endl;
    pvect(left, "Left"); cout << "eq" << endl; pvect(lfdwn,"Rot(270)Down");
    cout << endl;
    pvect(left, "Left"); cout << "eq" << endl; pvect(lfdwn2,"Rot(-90)Down");
    cout << endl;
    
    assert(veq(left,lfup));
    assert(veq(down,dfup));
    assert(veq(left,lfdwn));
    assert(veq(left,lfdwn2));
    
    double c75 = cos(75.0*pi/180.0);
    double s75 = sin(75.0*pi/180.0);
    double c75fs = swap_cosa_sina(s75);
    double s75fc = swap_cosa_sina(c75);
    
    cout << "--Cos/Sine swap tests" << endl
         << "sin(75) = " << s75 << endl
         << "cos(75) = " << c75 << endl
         << "cos(75) from sin(75) = " << c75fs
         << endl << "sin(75) from cos(75) = "  
         << s75fc << endl;
    
    assert(abs(c75 - c75fs) < 0.0000001);
    assert(abs(s75 - s75fc) < 0.0000001);
    
    vector<double> right(2,1.0);
    right[1] = 0.0;
    vector<double> upright(2,sqrt(2.0)/2.0);
    double s45 = sin(45.0*pi/180.0);
    double c270 = cos(270.0*pi/180.0);
    
    vector<double> dfright_cos = vrot_2D(right, c270, 2);
    vector<double> urfup = vrot_2D(up,-1.0*s45,4);
    
    cout << "--Rotations by evaluated trig functions" << endl;
    pvect(upright, "UpRight"); cout << "eq" << endl; 
    pvect(urfup,"Rot(sin(45))Up"); 
    cout << endl;
    pvect(down, "Down"); cout << "eq" << endl; 
    pvect(dfright_cos, "Rot(cos(270))Down");
    cout << endl;
    
    assert(veq(upright,urfup));
    assert(veq(down,dfright_cos));
    
    double s50 = sin(50);
    double s90m50 = sin(90-50);
    double c50 = cos(50);
    
    
    
    
    
    
    cout << "\n-------------Succeded--------------\n" << endl;

}
