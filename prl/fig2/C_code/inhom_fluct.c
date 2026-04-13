

/* 
  Trigonometric fluctuation landscape. Given int Nmode, double f0,
  and lists amps, coefs_x, phases_x, coefs_y, phases_y, returns:
    f = f0 + \sum_{i=1}^{Nmode} amps[i] * cos(coefs_x[i]*x + phases_x[i]) * cos(coefs_y[i]*y + phases_y[i]).

    Order of frequencies (x,y): (0,0), (1,0), (0,1), (2,0), (1,1), (0,2), (3,0), (2,1), (1,2), (0,3), etc.
*/
#ifdef FLUCT_TRIG
typedef struct fparam_trig{
  double f0;         // baseline value
  double Lx;         // length of system in x
  int Nmode;         // total number of modes
  double* amps;      // amplitudes of modes
  double* coefs_x;   // what to multiply x by inside cos(...)
  double* phases_x;  // phase shift to x
#ifdef TWOD
  double Ly;         // length of system in y
  double* coefs_y;   // what to multiply y by inside cos(...)
  double* phases_y;  // phase shift to y
#endif
} fparam_trig;

// suppose amps, phases_x, and phases_y are provided; compute coefs_x and coefs_y
// also scale phases by 2*pi
void compute_coefs_trig(fparam_trig* fp){
  fp[0].coefs_x       = (double*) calloc(fp[0].Nmode, sizeof(double));
  fp[0].coefs_y       = (double*) calloc(fp[0].Nmode, sizeof(double));
  for (i=0; i<fp[0].Nmode; i++){
    fp[0].phases_x[i] = fp[0].phases_x[i] * 2 * M_PI;
    fp[0].phases_y[i] = fp[0].phases_y[i] * 2 * M_PI;
  }
  int imax, ix, iy;
  int tot;
  imax=0;
  for (tot=1; imax<fp[0].Nmode; tot++){
    for (ix=tot; ix>=0; ix--){
      iy = tot-ix;
      fp[0].coefs_x[imax] = 2*M_PI*ix/fp[0].Lx;
      fp[0].coefs_y[imax] = 2*M_PI*iy/fp[0].Ly;
      imax += 1;
      if (imax>=fp[0].Nmode) break;
    }
  }
}


double fofx_trig(fparam_trig fp, double x
#ifdef TWOD
                 ,double y
#endif
  ){
  double f=fp.f0;    // the resulting fluctuation strength
  int i;             // mode index
  double tmp;

  for (i=0; i<fp.Nmode; i++){
    if ((fp.amps[i]<EPS) & (fp.amps[i]>-EPS)){
      continue;
    }
    tmp = fp.amps[i] * cos(fp.coefs_x[i]*x + fp.phases_x[i]);
#ifdef TWOD
    tmp = tmp * cos(fp.coefs_y[i]*y + fp.phases_y[i]);
#endif
    f +=  tmp;
  }
  return f;
}
#endif



/*
  piecewise fluctuation landscape w/linear interpolation, in 1 coordinate.
  xs lists the x coordinates of the reference points (in order from left to right). 
  fs lists their heights (in order).
  dfdxs lists their slopes (starting with (fs[1]-fs[0])/(xs[1]-xs[0]))
*/
#ifdef FLUCT_LIN
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
#endif


/*
  piecewise fluctuation landscape w/cubic interpolation, in 1 coordinate.
  xs lists the x coordinates of the reference points (in order from left to right). 
  If x is between xs[i] and xs[i+1], then returns:
      fs[i] + cs[i]*(x-xs[i])^2 + ds[i]*(x-xs[i])^3
*/
#ifdef FLUCT_CUB // piecewise fluctuation landscape w/cubic interpolation, in 1 coordinate
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

#endif



/*
  Triangular patch defined by corners (x1,y1),(x2,y2),(x3,y3).
  Inside patch, return f1. Outside, return f0.
*/
#ifdef FLUCT_POLY3
typedef struct fparam_poly3{
  double f0,f1;        // 2 levels of f
  double x1,y1;        // corner 1
  double x2,y2;        // corner 2
  double x3,y3;        // corner 3
} fparam_poly3;

double fofx_poly3(fparam_poly3 fp, double x, double y){
  if (inside_triangle(x,y,fp.x1,fp.y1,fp.x2,fp.y2,fp.x3,fp.y3)){
    return fp.f1;
  } else{
    return fp.f0;
  }
}
#endif


/*
  4-sided polygonal patch defined by corners (x1,y1),...,(x4,y4).
  Inside patch, return f1. Outside, return f0.
  Check if inside by checking if inside triangle (x1,y1),(x2,y2),(x3,y3) 
  or triangle (x3,y3),(x4,y4),(x1,y1).
*/
#ifdef FLUCT_POLY4
typedef struct fparam_poly4{
  double f0,f1;        // 2 levels of f
  double x1,y1;        // corner 1
  double x2,y2;        // corner 2
  double x3,y3;        // corner 3
  double x4,y4;        // corner 4
} fparam_poly4;

double fofx_poly4(fparam_poly4 fp, double x, double y){
  if (inside_triangle(x,y,fp.x1,fp.y1,fp.x2,fp.y2,fp.x3,fp.y3)){
    return fp.f1;
  } else if (inside_triangle(x,y,fp.x3,fp.y3,fp.x4,fp.y4,fp.x1,fp.y1)){
    return fp.f1;
  } else{
    return fp.f0;
  }
}
#endif



/*
  partial-annular patch defined by r1, r2, theta1, theta2.
*/
#ifdef FLUCT_ANNULUS
typedef struct fparam_annulus{
  double f0,f1;          // 2 levels of f
  double x0,y0;          // center around which radii are drawn
  double r1,r2;          // inner & outer radii
  double r12,r22;        // inner & outer radii squared
  double theta1,theta2;  // angular boundaries
  double dtheta;         // theta2-theta1, mod 2*pi)
} fparam_annulus;

double fofx_annulus(fparam_annulus fp, double x, double y){
  double dx,dy,r2,dtheta;
  dx = x-fp.x0;
  dy = y-fp.y0;
  r2 = dx*dx + dy*dy;
  if ((r2>fp.r12) & (r2<fp.r22)){
    dtheta = atan2(dy,dx) - fp.theta1;
    while (dtheta<0) dtheta += M_TWO_PI;
    while (dtheta>M_TWO_PI) dtheta -= M_TWO_PI;
    if (dtheta < fp.dtheta){
      return fp.f1;
    }
  }
  return fp.f0;
}
#endif