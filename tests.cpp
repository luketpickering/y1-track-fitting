#include <assert.h>
#include <iostream>
#include "read.hpp"

using std::cout;
using std::endl;

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
    
    cout << "\n-------------Succeded--------------\n" << endl;

}