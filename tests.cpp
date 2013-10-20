#include <assert.h>

void test_suite(){

    if (!quiet)
        cout << "-----------Running tests-----------\n\n" << flush;

    assert( sizeof(unsigned short) == 2);
    
    //Corresponds to 52-3-2
    unsigned short test = 0xD1A;
    if (!quiet){
    cout << "Read out data retrieval sanity check" << endl
         << "\tASSERT( " << (int)get_x(test) << " " << (int)get_y(test)
         << " " << (int)get_tdc(test) << " ) == ( 2 3 52 )\n" << endl;
    }
    
    assert( (int)get_x(test) == 2 );
    assert( (int)get_y(test) == 3 );
    assert( (int)get_tdc(test) == 52);
    
    if(!quiet)
        cout << "\n-------------Succeded--------------\n\n" << endl;

}