#include <iostream>
#include <sys/ioctl.h>

// ---------------------START USING STATEMENTS-----------------------
using std::endl;
using std::cout;
using std::flush;
// ----------------------END USING STATEMENTS------------------------

static unsigned char* snow;
static int sno_ctr = 0;
static int osc = 0;
static int cols = 80;
static int true_rows = 24;
static int rows;
static int row_len;
static int len;
static int ectr = 0;
static bool filled = false;

static void sno_shift_print(){
    cout << endl;
    for(int i = rows-1; i >= 0; --i)
        for(size_t j = 0; j < row_len; ++j)
            snow[(i+1)*row_len + j] = snow[i*row_len+j];

    for(size_t i = row_len; i < len; ++i)
    {
        if((i <= 2*row_len-1) || (i > (len-row_len-1))){
            for(size_t j = 0; j < 8; ++j)
                cout << "|";
            continue;
         }
        for(size_t j = 0; j < 8; ++j)
        {
            if(!(i%row_len) && j < 2){
                cout << "|";
                continue;
            }
            else if(!((i+1)%row_len) && j > 5){
                cout << "|";
                continue;
            }
            else if((snow[i] >> j)&1){
                cout << "*";
            }else{
                cout << " ";
            }
        }
    }
    cout << flush;   
}
void init_with_snow(){
    
    struct winsize w;
    ioctl(0,TIOCGWINSZ, &w);
    
    true_rows = w.ws_row;
    cols = w.ws_col;
    
    std::cerr << true_rows << " " << cols << endl;
    
    rows = true_rows+1;
    len = cols*rows/(sizeof(unsigned char)*8);
    row_len = cols/(sizeof(unsigned char)*8);
    snow = new unsigned char[len];
}
void read_with_snow(const unsigned char* dat, int arsize, int freq){
    
    ectr++;
    if(!filled || !(ectr%freq)){
        for(size_t i = osc; i < (arsize+osc) ; ++i){
            snow[sno_ctr] = dat[i%arsize];
            sno_ctr++;
            if(!filled){
                if (sno_ctr == len){
                    filled = true;
                    sno_ctr = 0;
                }else{
                    continue;
                }
            }else if(!(sno_ctr%row_len)){
                sno_shift_print();
                osc = i % arsize;
                sno_ctr = 0;
                break;
            }
        }
    }
}