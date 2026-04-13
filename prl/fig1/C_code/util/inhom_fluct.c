
/*
  piecewise fluctuation landscape w/linear interpolation, in 1 coordinate.
  xs lists the x coordinates of the reference points (in order from left to right). 
  fs lists their heights (in order).
  dfdxs lists their slopes (starting with (fs[1]-fs[0])/(xs[1]-xs[0]))
*/
typedef struct fparam_lin{
  int Nstep;         // number of modes
  double L;          // total system length
  double* xs;        // x values to interpolate between
  double* fs;        // f values to interpolate between
  double* dfdxs;     // (fs[i+1]-fs[i])/(xs[i+1]-xs[i])
} fparam_lin;

// suppose fs and xs are already provided; compute dfdxs
void compute_coefs_lin(fparam_lin* fp){
  double tmp;
  int i;
  fp[0].dfdxs = (double*) calloc(fp[0].Nstep, sizeof(double));
  for (i=0; i<fp[0].Nstep-1; i++){
    tmp = fp[0].xs[i+1] - fp[0].xs[i];
    if ((tmp<EPS) & (tmp>-EPS)) tmp = 1; // if there's no distance, don't divide by zero.
    fp[0].dfdxs[i]    = (fp[0].fs[i+1] - fp[0].fs[i]) / tmp;
  }
  tmp = fp[0].xs[0]+fp[0].L - fp[0].xs[i];
  if ((tmp<EPS) & (tmp>-EPS)) tmp = 1; // if there's no distance, don't divide by zero.
  fp[0].dfdxs[i]      = (fp[0].fs[0] - fp[0].fs[i]) / tmp;
}

double fofx_lin(fparam_lin fp, double x){
  int i;
  if (x<fp.xs[0]){
    return (x + fp.L - fp.xs[fp.Nstep-1]) * fp.dfdxs[fp.Nstep-1] + fp.fs[fp.Nstep-1];
  } else{
    for (i=1; i<fp.Nstep; i++){
      if (x<fp.xs[i]){
        return (x-fp.xs[i-1]) * fp.dfdxs[i-1] + fp.fs[i-1];
      }
    }
    return (x-fp.xs[fp.Nstep-1]) * fp.dfdxs[fp.Nstep-1] + fp.fs[fp.Nstep-1];
  }
}