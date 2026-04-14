

/*
  piecewise fluctuation landscape w/cubic interpolation, in 1 coordinate.
  xs lists the x coordinates of the reference points (in order from left to right). 
  If x is between xs[i] and xs[i+1], then returns:
      fs[i] + cs[i]*(x-xs[i])^2 + ds[i]*(x-xs[i])^3
*/
typedef struct fparam_cub{
  int Nstep;         // number of modes
  double L;          // total system length
  double* xs;        // x values to interpolate between
  double* fs;        // f values to interpolate between
  double* cs;        // -3*(fs[i] - fs[i+1])/(xs[i+1]-xs[i])^2
  double* ds;        // 2*(fs[i] - fs[i+1])/(xs[i+1]-xs[i])^3
} fparam_cub;

// suppose fs and xs are already provided; compute cs and ds
void compute_coefs_cub(fparam_cub* fp){
  double tmp;
  int i;
  fp[0].cs            = (double*) calloc(fp[0].Nstep, sizeof(double));
  fp[0].ds            = (double*) calloc(fp[0].Nstep, sizeof(double));
  for (i=0; i<fp[0].Nstep-1; i++){
    tmp                         = fp[0].xs[i+1] - fp[0].xs[i];
    if ((tmp<EPS) & (tmp>-EPS)) tmp = 1; // if there's no distance, don't divide by zero.
    fp[0].cs[i]       = -3 * (fp[0].fs[i] - fp[0].fs[i+1]) / pow(tmp,2);
    fp[0].ds[i]       =  2 * (fp[0].fs[i] - fp[0].fs[i+1]) / pow(tmp,3);
  }
  tmp                           = fp[0].xs[0] + fp[0].L - fp[0].xs[i];
  if ((tmp<EPS) & (tmp>-EPS)) tmp = 1; // if there's no distance, don't divide by zero.
  fp[0].ds[i]         =  2 * (fp[0].fs[i] - fp[0].fs[0]) / pow(tmp,3);
  fp[0].cs[i]         = -3 * (fp[0].fs[i] - fp[0].fs[0]) / pow(tmp,2);
}

double fofx_cub(fparam_cub fp, double x){
  int i,j;    // which interval it belongs to
  double dx;  // distance from left of interval
  if (x<fp.xs[0]){
    j = fp.Nstep-1;
    dx = x + fp.L - fp.xs[fp.Nstep-1];
  } else{
    j = fp.Nstep-1;
    dx = x-fp.xs[fp.Nstep-1];
    for (i=1; i<fp.Nstep; i++){
      if (x<fp.xs[i]){
        j = i-1;
        dx = x-fp.xs[i-1];
        break;
      }
    }
  }
  return fp.fs[j] + fp.cs[j]*pow(dx,2) + fp.ds[j]*pow(dx,3);
}

double dfdx_cub(fparam_cub fp, double x){
  int i,j;    // which interval it belongs to
  double dx;  // distance from left of interval
  if (x<fp.xs[0]){
    j = fp.Nstep-1;
    dx = x + fp.L - fp.xs[fp.Nstep-1];
  } else{
    j = fp.Nstep-1;
    dx = x-fp.xs[fp.Nstep-1];
    for (i=1; i<fp.Nstep; i++){
      if (x<fp.xs[i]){
        j = i-1;
        dx = x-fp.xs[i-1];
        break;
      }
    }
  }
  return 2*fp.cs[j]*dx + 3*fp.ds[j]*pow(dx,2);
}