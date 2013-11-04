#ifdef HIST   

const int max_cout = 1025;
const int bin_size = 5;
const int num_bins = max_cout/bin_size;

int hist[num_bins] = {0};
 
    while (get_event(event)) {
#ifdef VERBOSE
        if (! (num_events % 100000) ) {
            cout << "EVENT: " << num_events << endl;
            for (int i = 0; i < 8; i++) {
                print_ro( ( event+(3*i) ) );
            }
        }
#endif
        for (int i = 0; i < 8; i++) {
            unsigned short tdc = (event+(3*i))[2];
            int bin = floor(tdc/bin_size);
            if ( (bin <0) || (bin > (num_bins - 1)) ){
                cout << bin << endl;
                return -1;
            }
            hist[bin] += 1;
        }
        num_events++;
    }
    
    cout << "Found " << num_events << " events." << endl;
    
    int count = 0;
    for (int i = 0; i < num_bins; i++){
        cout << i << ", " << hist[i] << endl;
        count += hist[i];
        
    }
#endif     