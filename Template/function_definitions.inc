// Initialize process. 
  
%(process_class)s::%(process_class)s(const ParameterSet& configuration) {
    
    std::string param_card = configuration.get<std::string>("card");
    params.reset(new Parameters_%(model_name)s(SLHA::Reader(param_card)));
  
    %(constructor_lines)s

    %(finalstates_map)s
} 

//--------------------------------------------------------------------------
// Evaluate |M|^2, return a map of final states 

std::map< std::pair<int, int>, double > %(process_class)s::compute( const std::pair< std::vector<double>, std::vector<double> > &initialMomenta, const std::vector< std::pair<int, std::vector<double> > > &finalState ) { 

    // Set initial particle momenta
	momenta[0] = (double*) (&initialMomenta.first[0]);
    momenta[1] = (double*) (&initialMomenta.second[0]);

    // Suppose final particles are passed in the "correct" order
    std::vector<int> selectedFinalState(%(nexternal)s - 2);
    for (size_t index = 0; index < (%(nexternal)s - 2); index++) {  
        selectedFinalState[index] = finalState[index].first;
        momenta[index + 2] = (double*) (&finalState[index].second[0]);
    }

    // Set the event specific parameters 
    params->updateParameters();
    params->updateCouplings();

	// Initialise result object
    std::map< std::pair<int, int>, double > result;

    %(matrix_averaging)s

    return result;
}

//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------
// Evaluate |M|^2 for each subprocess

%(matrix_evaluations)s
