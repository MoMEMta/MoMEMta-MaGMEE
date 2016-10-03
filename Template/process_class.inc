    class %(process_class)s: public momemta::MatrixElement 
    {
        public:
    
            // Constructor & destructor
            %(process_class)s(const std::string& param_card);
            virtual ~%(process_class)s() {};
    
            // Calculate flavour-independent parts of cross section.
            virtual momemta::MatrixElement::Result compute(
                const std::pair < std::vector<double>, std::vector<double> > &initialMomenta, 
                const std::vector < std::pair < int, std::vector<double> > > &finalState);
    
            virtual std::shared_ptr<momemta::MEParameters> getParameters() {
                return params;
            }

            // needed? const std::vector<double>& getMasses() const {return mME;}  
    
        private:
    
            // default constructor should be hidden
            %(process_class)s() = delete;
    
            // list of helicities combinations
            %(helicity_matrix)s
    
            // Private functions to calculate the matrix element for all subprocesses
            %(all_wavefunction_definitions)s
            %(all_matrix_definitions)s
    
            // map of final states 
            std::map<std::vector<int>, std::vector< SubProcess< %(process_class)s > > > mapFinalStates;
    
            // Reference to the model parameters instance passed in the constructor
            std::shared_ptr<Parameters_%(model_name)s> params; 
    
            // vector with external particle masses
            std::vector<std::reference_wrapper<double>> mME;
    
            // vector with momenta (to be changed each event)
            double* momenta[%(nexternal)s];
    };
