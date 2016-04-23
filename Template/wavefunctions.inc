// Calculate all wavefunctions
static std::complex<double> w[%(nwavefuncs)s][18];

%(wavefunction_calls)s

// Calculate all amplitudes
%(amplitude_calls)s
