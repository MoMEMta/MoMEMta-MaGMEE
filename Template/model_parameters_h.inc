//==========================================================================
// This file has been automatically generated for C++
%(info_lines)s
//==========================================================================

#ifndef Parameters_%(model_name)s_H
#define Parameters_%(model_name)s_H

#include <complex>
#include <memory>

#include "read_slha.h"

class Parameters_%(model_name)s
{

public:

// Constructor: from path to param card 
Parameters_%(model_name)s(const std::string param_card, const bool verbose = false);

// Define "zero"
double zero, ZERO;
%(independent_parameters)s
%(independent_couplings)s
%(dependent_parameters)s
%(dependent_couplings)s

// Set parameters that are unchanged during the run
void setIndependentParameters();
// Set couplings that are unchanged during the run
void setIndependentCouplings();
// Set parameters that are changed event by event
void setDependentParameters();
// Set couplings that are changed event by event
void setDependentCouplings();

// Print parameters that are unchanged during the run
void printIndependentParameters();
// Print couplings that are unchanged during the run
void printIndependentCouplings();
// Print parameters that are changed event by event
void printDependentParameters();
// Print couplings that are changed event by event
void printDependentCouplings();

SLHAReader slha;

};

#endif // Parameters_%(model_name)s_H
