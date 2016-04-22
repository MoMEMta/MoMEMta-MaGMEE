double %(process_class)s::matrix_%(proc_name)s(%(matrix_args)s) { 

static std::complex<double> ztemp; 
static std::complex<double> jamp[%(ncolor)s];
// The color matrix
%(color_matrix_lines)s

// Calculate color flows
%(jamp_lines)s

// Sum and square the color flows to get the matrix element
double matrix = 0;
for(int i=0;i < %(ncolor)s; i++){
  ztemp = 0.;
  for(int j = 0; j < %(ncolor)s; j++)
    ztemp = ztemp + cf[i][j]*jamp[j];
  matrix = matrix+real(ztemp*conj(jamp[i]))/denom[i];
}

return matrix;
}
