// Define permutation
int perm[%(nexternal)s];
for(int i = 0; i < %(nexternal)s; i++) {
    perm[i]=i;
}

for(auto &me: mapFinalStates[selectedFinalState]) {
  
    double me_sum = 0;
    double me_mirror_sum = 0;
    
    for(int ihel = 0; ihel < %(ncomb)s; ihel++) {
        
        if(me.goodHel[ihel]) {
            
            double sum = 0.;
            calculate_wavefunctions(perm, helicities[ihel]);
            double meTemp = me.callback(*this); 
            sum += meTemp;
            me_sum += meTemp/me.denominator;

            if(me.hasMirrorProcess) {
                perm[0] = 1;
                perm[1] = 0;
                // Calculate wavefunctions
                calculate_wavefunctions(perm, helicities[ihel]);
                // Mirror back
                perm[0] = 0;
                perm[1] = 1;
                meTemp = me.callback(*this); 
                sum += meTemp;
                me_mirror_sum += meTemp/me.denominator;
            }

            if(!sum)
              me.goodHel[ihel] = false;
        }
    }

    for (auto const &initialState: me.initialStates) {
        result[initialState] = me_sum;
        if (me.hasMirrorProcess)
            result[std::make_pair(initialState.second, initialState.first)] = me_mirror_sum;
    }
}
