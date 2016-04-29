//==========================================================================
// Class member functions for calculating the matrix elements for
%(process_lines)s

//--------------------------------------------------------------------------
// Initialize process. 
  
%(process_class)s::%(process_class)s(Parameters_%(model_name)s &params):
    params(params)
    {
    %(initProc_lines)s

    %(finalstates_map)s
} 

//--------------------------------------------------------------------------
// Evaluate |M|^2, return a map of final states 

std::map< std::pair<int, int>, double > %(process_class)s::sigmaKin( const std::vector< std::vector<double> > &initialMomenta, const std::vector< std::pair<int, std::vector<double> > > &finalState ) { 

    // Set initial particle momenta
	momenta[0] = (double*) (&initialMomenta[0][0]);
    momenta[1] = (double*) (&initialMomenta[1][0]);

    // Suppose final particles are passed in the "correct" order
    std::vector<int> selectedFinalState(%(nexternal)s - 2);
    size_t index = 2;
    for (auto const &finalPart: finalState){
      selectedFinalState[index - 2] = finalPart.first;
      momenta[index] = (double*) (&finalPart.second[0]);
      index++;
    }

    // Set the event specific parameters 
    params.setDependentParameters();
    params.setDependentCouplings();

	// Initialise result object
    std::map<std::pair<int, int>, double> result;

    %(sigmaKin_lines)s

    return result;
}

//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------
// Evaluate |M|^2 for each subprocess

%(all_sigmaKin)s
