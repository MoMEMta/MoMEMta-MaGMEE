//==========================================================================
// A class for calculating the matrix elements for
%(process_lines)s
//--------------------------------------------------------------------------

//Forward declaration needed to template correctly __MatrixElement in the class
class %(process_class)s;

class %(process_class)s: public CPPProcess
{
  public:

    // Constructor & destructor
    %(process_class)s(Parameters_%(model_name)s &params);
    virtual ~%(process_class)s() {};

    // Calculate flavour-independent parts of cross section.
    virtual std::map< std::pair<int, int>, double> sigmaKin(const std::vector<std::vector<double>> &initialMomenta, const std::vector< std::pair<int, std::vector<double> > > &finalState);

    // Info on the subprocess.
    virtual std::string name()    const {return "%(process_definition)s";}

    const std::vector<double>& getMasses() const {return mME;}  

  private:

    // default constructor should be hidden
    %(process_class)s();

    // list of helicities combinations
    %(helicity_matrix)s

    // Private functions to calculate the matrix element for all subprocesses
    %(all_sigma_kin_definitions)s
    %(all_matrix_definitions)s

    // map of final states 
    std::map< std::vector<int>, std::vector< __MatrixElement<%(process_class)s> > > mapFinalStates;

    // Reference to the model parameters instance passed in the constructor
    Parameters_%(model_name)s& params;

    // vector with external particle masses
    std::vector<double> mME;

    // vector with momenta (to be changed each event)
    double* momenta[%(nexternal)s];
};

