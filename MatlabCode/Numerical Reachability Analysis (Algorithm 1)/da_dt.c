/* This file was automatically generated by CasADi.
   The CasADi copyright holders make no ownership claim of its contents. */
#ifdef __cplusplus
extern "C" {
#endif

/* How to prefix internal symbols */
#ifdef CASADI_CODEGEN_PREFIX
  #define CASADI_NAMESPACE_CONCAT(NS, ID) _CASADI_NAMESPACE_CONCAT(NS, ID)
  #define _CASADI_NAMESPACE_CONCAT(NS, ID) NS ## ID
  #define CASADI_PREFIX(ID) CASADI_NAMESPACE_CONCAT(CODEGEN_PREFIX, ID)
#else
  #define CASADI_PREFIX(ID) da_dt_ ## ID
#endif

#include <math.h>
#include <stdio.h>
#include <string.h>
#ifdef MATLAB_MEX_FILE
#include <mex.h>
#endif

#ifndef casadi_real
#define casadi_real double
#endif

#ifndef casadi_int
#define casadi_int long long int
#endif

/* Add prefix to internal symbols */
#define casadi_f0 CASADI_PREFIX(f0)
#define casadi_fill CASADI_PREFIX(fill)
#define casadi_from_mex CASADI_PREFIX(from_mex)
#define casadi_s0 CASADI_PREFIX(s0)
#define casadi_s1 CASADI_PREFIX(s1)
#define casadi_s2 CASADI_PREFIX(s2)
#define casadi_sq CASADI_PREFIX(sq)
#define casadi_to_mex CASADI_PREFIX(to_mex)

/* Symbol visibility in DLLs */
#ifndef CASADI_SYMBOL_EXPORT
  #if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)
    #if defined(STATIC_LINKED)
      #define CASADI_SYMBOL_EXPORT
    #else
      #define CASADI_SYMBOL_EXPORT __declspec(dllexport)
    #endif
  #elif defined(__GNUC__) && defined(GCC_HASCLASSVISIBILITY)
    #define CASADI_SYMBOL_EXPORT __attribute__ ((visibility ("default")))
  #else
    #define CASADI_SYMBOL_EXPORT
  #endif
#endif

casadi_real casadi_sq(casadi_real x) { return x*x;}

void casadi_fill(casadi_real* x, casadi_int n, casadi_real alpha) {
  casadi_int i;
  if (x) {
    for (i=0; i<n; ++i) *x++ = alpha;
  }
}

#ifdef MATLAB_MEX_FILE
casadi_real* casadi_from_mex(const mxArray* p, casadi_real* y, const casadi_int* sp, casadi_real* w) {
  casadi_int nrow, ncol, is_sparse, c, k, p_nrow, p_ncol;
  const casadi_int *colind, *row;
  mwIndex *Jc, *Ir;
  const double* p_data;
  if (!mxIsDouble(p) || mxGetNumberOfDimensions(p)!=2)
    mexErrMsgIdAndTxt("Casadi:RuntimeError",
      "\"from_mex\" failed: Not a two-dimensional matrix of double precision.");
  nrow = *sp++;
  ncol = *sp++;
  colind = sp;
  row = sp+ncol+1;
  p_nrow = mxGetM(p);
  p_ncol = mxGetN(p);
  is_sparse = mxIsSparse(p);
  Jc = 0;
  Ir = 0;
  if (is_sparse) {
    Jc = mxGetJc(p);
    Ir = mxGetIr(p);
  }
  p_data = (const double*)mxGetData(p);
  if (p_nrow==1 && p_ncol==1) {
    casadi_int nnz;
    double v = is_sparse && Jc[1]==0 ? 0 : *p_data;
    nnz = sp[ncol];
    casadi_fill(y, nnz, v);
  } else {
    casadi_int tr = 0;
    if (nrow!=p_nrow || ncol!=p_ncol) {
      tr = nrow==p_ncol && ncol==p_nrow && (nrow==1 || ncol==1);
      if (!tr) mexErrMsgIdAndTxt("Casadi:RuntimeError",
                                 "\"from_mex\" failed: Dimension mismatch. "
                                 "Expected %d-by-%d, got %d-by-%d instead.",
                                 nrow, ncol, p_nrow, p_ncol);
    }
    if (is_sparse) {
      if (tr) {
        for (c=0; c<ncol; ++c)
          for (k=colind[c]; k<colind[c+1]; ++k) w[row[k]+c*nrow]=0;
        for (c=0; c<p_ncol; ++c)
          for (k=Jc[c]; k<(casadi_int) Jc[c+1]; ++k) w[c+Ir[k]*p_ncol] = p_data[k];
        for (c=0; c<ncol; ++c)
          for (k=colind[c]; k<colind[c+1]; ++k) y[k] = w[row[k]+c*nrow];
      } else {
        for (c=0; c<ncol; ++c) {
          for (k=colind[c]; k<colind[c+1]; ++k) w[row[k]]=0;
          for (k=Jc[c]; k<(casadi_int) Jc[c+1]; ++k) w[Ir[k]]=p_data[k];
          for (k=colind[c]; k<colind[c+1]; ++k) y[k]=w[row[k]];
        }
      }
    } else {
      for (c=0; c<ncol; ++c) {
        for (k=colind[c]; k<colind[c+1]; ++k) {
          y[k] = p_data[row[k]+c*nrow];
        }
      }
    }
  }
  return y;
}

#endif

#define casadi_to_double(x) ((double) x)

#ifdef MATLAB_MEX_FILE
mxArray* casadi_to_mex(const casadi_int* sp, const casadi_real* x) {
  casadi_int nrow, ncol, c, k;
#ifndef CASADI_MEX_NO_SPARSE
  casadi_int nnz;
#endif
  const casadi_int *colind, *row;
  mxArray *p;
  double *d;
#ifndef CASADI_MEX_NO_SPARSE
  casadi_int i;
  mwIndex *j;
#endif /* CASADI_MEX_NO_SPARSE */
  nrow = *sp++;
  ncol = *sp++;
  colind = sp;
  row = sp+ncol+1;
#ifndef CASADI_MEX_NO_SPARSE
  nnz = sp[ncol];
  if (nnz!=nrow*ncol) {
    p = mxCreateSparse(nrow, ncol, nnz, mxREAL);
    for (i=0, j=mxGetJc(p); i<=ncol; ++i) *j++ = *colind++;
    for (i=0, j=mxGetIr(p); i<nnz; ++i) *j++ = *row++;
    if (x) {
      d = (double*)mxGetData(p);
      for (i=0; i<nnz; ++i) *d++ = casadi_to_double(*x++);
    }
    return p;
  }
#endif /* CASADI_MEX_NO_SPARSE */
  p = mxCreateDoubleMatrix(nrow, ncol, mxREAL);
  if (x) {
    d = (double*)mxGetData(p);
    for (c=0; c<ncol; ++c) {
      for (k=colind[c]; k<colind[c+1]; ++k) {
        d[row[k]+c*nrow] = casadi_to_double(*x++);
      }
    }
  }
  return p;
}

#endif

#ifndef CASADI_PRINTF
#ifdef MATLAB_MEX_FILE
  #define CASADI_PRINTF mexPrintf
#else
  #define CASADI_PRINTF printf
#endif
#endif

static const casadi_int casadi_s0[7] = {3, 1, 0, 3, 0, 1, 2};
static const casadi_int casadi_s1[8] = {4, 1, 0, 4, 0, 1, 2, 3};
static const casadi_int casadi_s2[5] = {1, 1, 0, 1, 0};

/* da_dt:(r[3],v[3],q[4],t)->(adot) */
static int casadi_f0(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem) {
  casadi_real a0, a1, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a2, a20, a21, a22, a23, a24, a25, a26, a27, a28, a29, a3, a30, a31, a32, a33, a34, a35, a36, a37, a38, a4, a5, a6, a7, a8, a9;
  a0=2.;
  a1=arg[0]? arg[0][1] : 0;
  a2=arg[1]? arg[1][2] : 0;
  a3=(a1*a2);
  a4=arg[0]? arg[0][2] : 0;
  a5=arg[1]? arg[1][1] : 0;
  a6=(a4*a5);
  a3=(a3-a6);
  a3=casadi_sq(a3);
  a6=arg[1]? arg[1][0] : 0;
  a7=(a4*a6);
  a8=arg[0]? arg[0][0] : 0;
  a9=(a8*a2);
  a7=(a7-a9);
  a7=casadi_sq(a7);
  a3=(a3+a7);
  a7=(a8*a5);
  a9=(a1*a6);
  a7=(a7-a9);
  a7=casadi_sq(a7);
  a3=(a3+a7);
  a7=3.9860044150000002e+05;
  a3=(a3/a7);
  a9=1.;
  a10=2.5087779537745443e-06;
  a11=casadi_sq(a6);
  a12=casadi_sq(a5);
  a11=(a11+a12);
  a12=casadi_sq(a2);
  a11=(a11+a12);
  a12=casadi_sq(a8);
  a13=casadi_sq(a1);
  a12=(a12+a13);
  a13=casadi_sq(a4);
  a12=(a12+a13);
  a12=sqrt(a12);
  a7=(a7/a12);
  a11=(a11-a7);
  a7=(a11*a8);
  a12=(a8*a6);
  a13=(a1*a5);
  a12=(a12+a13);
  a13=(a4*a2);
  a12=(a12+a13);
  a13=(a12*a6);
  a7=(a7-a13);
  a7=(a10*a7);
  a13=casadi_sq(a7);
  a14=(a11*a1);
  a15=(a12*a5);
  a14=(a14-a15);
  a14=(a10*a14);
  a15=casadi_sq(a14);
  a13=(a13+a15);
  a11=(a11*a4);
  a12=(a12*a2);
  a11=(a11-a12);
  a10=(a10*a11);
  a11=casadi_sq(a10);
  a13=(a13+a11);
  a13=(a9-a13);
  a3=(a3/a13);
  a13=casadi_sq(a3);
  a13=(a0*a13);
  a11=(a1*a2);
  a12=(a4*a5);
  a11=(a11-a12);
  a11=casadi_sq(a11);
  a12=(a4*a6);
  a15=(a8*a2);
  a12=(a12-a15);
  a12=casadi_sq(a12);
  a11=(a11+a12);
  a12=(a8*a5);
  a15=(a1*a6);
  a12=(a12-a15);
  a12=casadi_sq(a12);
  a11=(a11+a12);
  a11=sqrt(a11);
  a13=(a13/a11);
  a11=casadi_sq(a7);
  a12=casadi_sq(a14);
  a11=(a11+a12);
  a12=casadi_sq(a10);
  a11=(a11+a12);
  a12=sqrt(a11);
  a15=(a9-a12);
  a16=(a9+a12);
  a15=(a15/a16);
  a15=sqrt(a15);
  a16=(a8*a6);
  a17=(a1*a5);
  a16=(a16+a17);
  a17=(a4*a2);
  a16=(a16+a17);
  a17=0.;
  a16=(a16<a17);
  a18=6.2831853071795862e+00;
  a19=(a7*a8);
  a20=(a14*a1);
  a19=(a19+a20);
  a20=(a10*a4);
  a19=(a19+a20);
  a7=casadi_sq(a7);
  a14=casadi_sq(a14);
  a7=(a7+a14);
  a10=casadi_sq(a10);
  a7=(a7+a10);
  a7=sqrt(a7);
  a10=casadi_sq(a8);
  a14=casadi_sq(a1);
  a10=(a10+a14);
  a14=casadi_sq(a4);
  a10=(a10+a14);
  a10=sqrt(a10);
  a7=(a7*a10);
  a19=(a19/a7);
  a19=acos(a19);
  a18=(a18-a19);
  a18=(a16?a18:0);
  a16=(!a16);
  a16=(a16?a19:0);
  a18=(a18+a16);
  a18=(a18/a0);
  a18=tan(a18);
  a15=(a15*a18);
  a15=atan(a15);
  a15=(a0*a15);
  a18=sin(a15);
  a18=(a12*a18);
  a15=(a15-a18);
  a15=sin(a15);
  a12=(a12*a15);
  a12=(a13*a12);
  a15=casadi_sq(a8);
  a18=casadi_sq(a1);
  a15=(a15+a18);
  a18=casadi_sq(a4);
  a15=(a15+a18);
  a15=sqrt(a15);
  a18=(a8/a15);
  a16=arg[2]? arg[2][1] : 0;
  a19=arg[2]? arg[2][0] : 0;
  a7=casadi_sq(a19);
  a10=casadi_sq(a16);
  a7=(a7+a10);
  a10=arg[2]? arg[2][2] : 0;
  a14=casadi_sq(a10);
  a7=(a7+a14);
  a14=arg[2]? arg[2][3] : 0;
  a20=casadi_sq(a14);
  a7=(a7+a20);
  a7=sqrt(a7);
  a16=(a16/a7);
  a14=(a14/a7);
  a20=(a16*a14);
  a19=(a19/a7);
  a10=(a10/a7);
  a7=(a19*a10);
  a20=(a20+a7);
  a20=(a0*a20);
  a7=1.7453292519943295e-02;
  a21=1.7487317400000001e+02;
  a22=-2.4109079999999999e-01;
  a23=2.4562962500000000e+06;
  a24=arg[3]? arg[3][0] : 0;
  a25=86400.;
  a24=(a24/a25);
  a23=(a23+a24);
  a24=2451545.;
  a23=(a23-a24);
  a24=36525.;
  a23=(a23/a24);
  a22=(a22*a23);
  a21=(a21+a22);
  a22=4.0670000000000002e-05;
  a24=casadi_sq(a23);
  a22=(a22*a24);
  a21=(a21+a22);
  a22=-1.3270000000000000e-06;
  a25=casadi_sq(a23);
  a25=(a25*a23);
  a22=(a22*a25);
  a21=(a21+a22);
  a21=(a7*a21);
  a22=(-a21);
  a22=cos(a22);
  a26=1.0293734800000000e+02;
  a27=3.2255570000000000e-01;
  a27=(a27*a23);
  a26=(a26+a27);
  a27=1.5025999999999999e-04;
  a27=(a27*a24);
  a26=(a26+a27);
  a27=4.7800000000000002e-07;
  a27=(a27*a25);
  a26=(a26+a27);
  a26=(a7*a26);
  a27=(a26-a21);
  a28=(-a27);
  a28=cos(a28);
  a29=(a22*a28);
  a30=(-a21);
  a30=sin(a30);
  a31=1.3054600000000000e-02;
  a31=(a31*a23);
  a32=-9.3100000000000006e-06;
  a32=(a32*a24);
  a31=(a31+a32);
  a32=-3.4000000000000000e-08;
  a32=(a32*a25);
  a31=(a31+a32);
  a31=(a7*a31);
  a32=(-a31);
  a32=cos(a32);
  a30=(a30*a32);
  a33=(-a27);
  a33=sin(a33);
  a34=(a30*a33);
  a29=(a29-a34);
  a34=1.4959802299063239e+08;
  a35=1.6708620000000000e-02;
  a36=-4.2036999999999997e-05;
  a36=(a36*a23);
  a35=(a35+a36);
  a36=-1.2360000000000000e-07;
  a36=(a36*a24);
  a35=(a35+a36);
  a36=3.9999999999999998e-11;
  a36=(a36*a25);
  a35=(a35+a36);
  a36=casadi_sq(a35);
  a36=(a9-a36);
  a34=(a34*a36);
  a36=1.0046644900000000e+02;
  a25=3.5999372851900000e+04;
  a25=(a25*a23);
  a36=(a36+a25);
  a25=-5.6799999999999998e-06;
  a25=(a25*a24);
  a36=(a36+a25);
  a7=(a7*a36);
  a7=(a7-a26);
  a26=(a0*a35);
  a36=casadi_sq(a35);
  a36=(a35*a36);
  a25=4.;
  a36=(a36/a25);
  a26=(a26-a36);
  a36=5.2083333333333336e-02;
  a24=casadi_sq(a35);
  a24=casadi_sq(a24);
  a24=(a35*a24);
  a36=(a36*a24);
  a26=(a26+a36);
  a36=sin(a7);
  a26=(a26*a36);
  a36=1.2500000000000000e+00;
  a24=casadi_sq(a35);
  a36=(a36*a24);
  a24=4.5833333333333331e-01;
  a23=casadi_sq(a35);
  a23=casadi_sq(a23);
  a24=(a24*a23);
  a36=(a36-a24);
  a24=(a0*a7);
  a24=sin(a24);
  a36=(a36*a24);
  a26=(a26+a36);
  a36=1.0833333333333333e+00;
  a24=casadi_sq(a35);
  a24=(a35*a24);
  a36=(a36*a24);
  a24=6.7187500000000000e-01;
  a23=casadi_sq(a35);
  a23=casadi_sq(a23);
  a23=(a35*a23);
  a24=(a24*a23);
  a36=(a36-a24);
  a24=3.;
  a24=(a24*a7);
  a24=sin(a24);
  a36=(a36*a24);
  a26=(a26+a36);
  a36=1.0729166666666667e+00;
  a24=casadi_sq(a35);
  a24=casadi_sq(a24);
  a36=(a36*a24);
  a25=(a25*a7);
  a25=sin(a25);
  a36=(a36*a25);
  a26=(a26+a36);
  a36=1.1427083333333334e+00;
  a25=casadi_sq(a35);
  a25=casadi_sq(a25);
  a25=(a35*a25);
  a36=(a36*a25);
  a25=5.;
  a25=(a25*a7);
  a25=sin(a25);
  a36=(a36*a25);
  a26=(a26+a36);
  a7=(a7+a26);
  a26=cos(a7);
  a35=(a35*a26);
  a35=(a9+a35);
  a34=(a34/a35);
  a35=cos(a7);
  a35=(a34*a35);
  a29=(a29*a35);
  a26=(-a27);
  a26=sin(a26);
  a22=(a22*a26);
  a27=(-a27);
  a27=cos(a27);
  a30=(a30*a27);
  a22=(a22+a30);
  a7=sin(a7);
  a34=(a34*a7);
  a22=(a22*a34);
  a29=(a29+a22);
  a29=(a29+a8);
  a22=casadi_sq(a29);
  a7=9.1748200035787253e-01;
  a30=(-a21);
  a30=cos(a30);
  a30=(a30*a32);
  a32=(a30*a27);
  a21=(-a21);
  a21=sin(a21);
  a26=(a21*a26);
  a32=(a32-a26);
  a32=(a32*a34);
  a21=(a21*a28);
  a30=(a30*a33);
  a21=(a21+a30);
  a21=(a21*a35);
  a32=(a32-a21);
  a21=(a7*a32);
  a30=-3.9777729827042280e-01;
  a31=(-a31);
  a31=sin(a31);
  a33=(a31*a33);
  a33=(a33*a35);
  a31=(a31*a27);
  a31=(a31*a34);
  a33=(a33-a31);
  a30=(a30*a33);
  a21=(a21+a30);
  a21=(a21+a1);
  a30=casadi_sq(a21);
  a22=(a22+a30);
  a30=3.9777729827042280e-01;
  a30=(a30*a32);
  a7=(a7*a33);
  a30=(a30+a7);
  a30=(a30+a4);
  a7=casadi_sq(a30);
  a22=(a22+a7);
  a22=sqrt(a22);
  a7=(a29/a22);
  a20=(a20*a7);
  a33=(a10*a14);
  a32=(a19*a16);
  a33=(a33-a32);
  a33=(a0*a33);
  a32=(a21/a22);
  a33=(a33*a32);
  a20=(a20+a33);
  a33=casadi_sq(a19);
  a31=casadi_sq(a16);
  a33=(a33-a31);
  a31=casadi_sq(a10);
  a33=(a33-a31);
  a31=casadi_sq(a14);
  a33=(a33+a31);
  a22=(a30/a22);
  a33=(a33*a22);
  a20=(a20+a33);
  a33=(-a20);
  a33=(a33<a17);
  a33=(!a33);
  a17=1.0000000000000000e-03;
  a31=1.6974242894395963e-04;
  a34=1.4959787069999999e+08;
  a29=(a29/a34);
  a29=casadi_sq(a29);
  a21=(a21/a34);
  a21=casadi_sq(a21);
  a29=(a29+a21);
  a30=(a30/a34);
  a30=casadi_sq(a30);
  a29=(a29+a30);
  a31=(a31/a29);
  a17=(a17*a31);
  a17=(a17*a20);
  a31=8.0000000000000004e-01;
  a29=casadi_sq(a19);
  a30=casadi_sq(a16);
  a29=(a29+a30);
  a30=casadi_sq(a10);
  a29=(a29-a30);
  a30=casadi_sq(a14);
  a29=(a29-a30);
  a29=(a29*a7);
  a30=(a16*a10);
  a34=(a19*a14);
  a30=(a30+a34);
  a30=(a0*a30);
  a30=(a30*a32);
  a29=(a29+a30);
  a30=(a16*a14);
  a34=(a19*a10);
  a30=(a30-a34);
  a30=(a0*a30);
  a30=(a30*a22);
  a29=(a29+a30);
  a29=(a31*a29);
  a29=(a17*a29);
  a30=casadi_sq(a19);
  a34=casadi_sq(a16);
  a30=(a30-a34);
  a34=casadi_sq(a10);
  a30=(a30+a34);
  a34=casadi_sq(a14);
  a30=(a30-a34);
  a34=casadi_sq(a19);
  a21=casadi_sq(a16);
  a34=(a34-a21);
  a21=casadi_sq(a10);
  a34=(a34-a21);
  a21=casadi_sq(a14);
  a34=(a34+a21);
  a21=(a30*a34);
  a27=(a10*a14);
  a35=(a19*a16);
  a27=(a27+a35);
  a27=(a0*a27);
  a35=(a10*a14);
  a28=(a19*a16);
  a35=(a35-a28);
  a35=(a0*a35);
  a28=(a27*a35);
  a21=(a21-a28);
  a28=casadi_sq(a19);
  a26=casadi_sq(a16);
  a28=(a28+a26);
  a26=casadi_sq(a10);
  a28=(a28-a26);
  a26=casadi_sq(a14);
  a28=(a28-a26);
  a26=(a30*a34);
  a36=(a27*a35);
  a26=(a26-a36);
  a26=(a28*a26);
  a36=(a16*a10);
  a25=(a19*a14);
  a36=(a36+a25);
  a36=(a0*a36);
  a25=(a16*a10);
  a24=(a19*a14);
  a25=(a25-a24);
  a25=(a0*a25);
  a24=(a25*a34);
  a23=(a16*a14);
  a37=(a19*a10);
  a23=(a23+a37);
  a23=(a0*a23);
  a37=(a27*a23);
  a24=(a24-a37);
  a24=(a36*a24);
  a26=(a26-a24);
  a24=(a16*a14);
  a37=(a19*a10);
  a24=(a24-a37);
  a24=(a0*a24);
  a37=(a25*a35);
  a38=(a30*a23);
  a37=(a37-a38);
  a37=(a24*a37);
  a26=(a26+a37);
  a21=(a21/a26);
  a21=(a29*a21);
  a21=(-a21);
  a21=(a33?a21:0);
  a37=(a16*a10);
  a38=(a19*a14);
  a37=(a37-a38);
  a37=(a0*a37);
  a37=(a37*a7);
  a7=casadi_sq(a19);
  a38=casadi_sq(a16);
  a7=(a7-a38);
  a38=casadi_sq(a10);
  a7=(a7+a38);
  a38=casadi_sq(a14);
  a7=(a7-a38);
  a7=(a7*a32);
  a37=(a37+a7);
  a10=(a10*a14);
  a19=(a19*a16);
  a10=(a10+a19);
  a0=(a0*a10);
  a0=(a0*a22);
  a37=(a37+a0);
  a37=(a31*a37);
  a37=(a17*a37);
  a0=(a36*a34);
  a22=(a24*a35);
  a0=(a0-a22);
  a0=(a0/a26);
  a0=(a37*a0);
  a0=(a33?a0:0);
  a21=(a21+a0);
  a0=5.0000000000000000e-01;
  a22=4.0000000000000002e-01;
  a22=(a22*a20);
  a0=(a0-a22);
  a31=(a31*a20);
  a0=(a0-a31);
  a17=(a17*a0);
  a0=(a36*a27);
  a31=(a24*a30);
  a0=(a0-a31);
  a0=(a0/a26);
  a0=(a17*a0);
  a0=(a33?a0:0);
  a21=(a21+a0);
  a0=(-a29);
  a0=(a33?a0:0);
  a21=(a21+a0);
  a0=(a18*a21);
  a31=(a1/a15);
  a20=(a25*a34);
  a22=(a27*a23);
  a20=(a20-a22);
  a20=(a20/a26);
  a20=(a29*a20);
  a20=(a33?a20:0);
  a34=(a28*a34);
  a22=(a24*a23);
  a34=(a34-a22);
  a34=(a34/a26);
  a34=(a37*a34);
  a34=(-a34);
  a34=(a33?a34:0);
  a20=(a20+a34);
  a27=(a28*a27);
  a24=(a24*a25);
  a27=(a27-a24);
  a27=(a27/a26);
  a27=(a17*a27);
  a27=(-a27);
  a27=(a33?a27:0);
  a20=(a20+a27);
  a27=(-a37);
  a27=(a33?a27:0);
  a20=(a20+a27);
  a27=(a31*a20);
  a0=(a0+a27);
  a15=(a4/a15);
  a27=(a25*a35);
  a24=(a30*a23);
  a27=(a27-a24);
  a27=(a27/a26);
  a29=(a29*a27);
  a29=(-a29);
  a29=(a33?a29:0);
  a35=(a28*a35);
  a23=(a36*a23);
  a35=(a35-a23);
  a35=(a35/a26);
  a37=(a37*a35);
  a37=(a33?a37:0);
  a29=(a29+a37);
  a28=(a28*a30);
  a36=(a36*a25);
  a28=(a28-a36);
  a28=(a28/a26);
  a28=(a17*a28);
  a28=(a33?a28:0);
  a29=(a29+a28);
  a33=(a33?a17:0);
  a29=(a29+a33);
  a33=(a15*a29);
  a0=(a0+a33);
  a12=(a12*a0);
  a9=(a9-a11);
  a3=(a3*a9);
  a9=casadi_sq(a8);
  a11=casadi_sq(a1);
  a9=(a9+a11);
  a11=casadi_sq(a4);
  a9=(a9+a11);
  a9=sqrt(a9);
  a3=(a3/a9);
  a13=(a13*a3);
  a3=(a4*a6);
  a9=(a8*a2);
  a3=(a3-a9);
  a2=(a1*a2);
  a4=(a4*a5);
  a2=(a2-a4);
  a4=casadi_sq(a2);
  a9=casadi_sq(a3);
  a4=(a4+a9);
  a8=(a8*a5);
  a1=(a1*a6);
  a8=(a8-a1);
  a1=casadi_sq(a8);
  a4=(a4+a1);
  a4=sqrt(a4);
  a3=(a3/a4);
  a1=(a3*a15);
  a8=(a8/a4);
  a6=(a8*a31);
  a1=(a1-a6);
  a6=casadi_sq(a1);
  a8=(a8*a18);
  a2=(a2/a4);
  a15=(a2*a15);
  a8=(a8-a15);
  a15=casadi_sq(a8);
  a6=(a6+a15);
  a2=(a2*a31);
  a3=(a3*a18);
  a2=(a2-a3);
  a3=casadi_sq(a2);
  a6=(a6+a3);
  a6=sqrt(a6);
  a1=(a1/a6);
  a1=(a1*a21);
  a8=(a8/a6);
  a8=(a8*a20);
  a1=(a1+a8);
  a2=(a2/a6);
  a2=(a2*a29);
  a1=(a1+a2);
  a13=(a13*a1);
  a12=(a12+a13);
  if (res[0]!=0) res[0][0]=a12;
  return 0;
}

CASADI_SYMBOL_EXPORT int da_dt(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem){
  return casadi_f0(arg, res, iw, w, mem);
}

CASADI_SYMBOL_EXPORT int da_dt_alloc_mem(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT int da_dt_init_mem(int mem) {
  return 0;
}

CASADI_SYMBOL_EXPORT void da_dt_free_mem(int mem) {
}

CASADI_SYMBOL_EXPORT int da_dt_checkout(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT void da_dt_release(int mem) {
}

CASADI_SYMBOL_EXPORT void da_dt_incref(void) {
}

CASADI_SYMBOL_EXPORT void da_dt_decref(void) {
}

CASADI_SYMBOL_EXPORT casadi_int da_dt_n_in(void) { return 4;}

CASADI_SYMBOL_EXPORT casadi_int da_dt_n_out(void) { return 1;}

CASADI_SYMBOL_EXPORT casadi_real da_dt_default_in(casadi_int i){
  switch (i) {
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* da_dt_name_in(casadi_int i){
  switch (i) {
    case 0: return "r";
    case 1: return "v";
    case 2: return "q";
    case 3: return "t";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* da_dt_name_out(casadi_int i){
  switch (i) {
    case 0: return "adot";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* da_dt_sparsity_in(casadi_int i) {
  switch (i) {
    case 0: return casadi_s0;
    case 1: return casadi_s0;
    case 2: return casadi_s1;
    case 3: return casadi_s2;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* da_dt_sparsity_out(casadi_int i) {
  switch (i) {
    case 0: return casadi_s2;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT int da_dt_work(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w) {
  if (sz_arg) *sz_arg = 4;
  if (sz_res) *sz_res = 1;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 0;
  return 0;
}

#ifdef MATLAB_MEX_FILE
void mex_da_dt(int resc, mxArray *resv[], int argc, const mxArray *argv[]) {
  casadi_int i;
  casadi_real w[51];
  casadi_int *iw = 0;
  const casadi_real* arg[4] = {0};
  casadi_real* res[1] = {0};
  if (argc>4) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"da_dt\" failed. Too many input arguments (%d, max 4)", argc);
  if (resc>1) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"da_dt\" failed. Too many output arguments (%d, max 1)", resc);
  if (--argc>=0) arg[0] = casadi_from_mex(argv[0], w, casadi_s0, w+12);
  if (--argc>=0) arg[1] = casadi_from_mex(argv[1], w+3, casadi_s0, w+12);
  if (--argc>=0) arg[2] = casadi_from_mex(argv[2], w+6, casadi_s1, w+12);
  if (--argc>=0) arg[3] = casadi_from_mex(argv[3], w+10, casadi_s2, w+12);
  --resc;
  res[0] = w+11;
  i = da_dt(arg, res, iw, w+12, 0);
  if (i) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"da_dt\" failed.");
  if (res[0]) resv[0] = casadi_to_mex(casadi_s2, res[0]);
}
#endif

casadi_int main_da_dt(casadi_int argc, char* argv[]) {
  casadi_int j;
  casadi_real* a;
  const casadi_real* r;
  casadi_int flag;
  casadi_int *iw = 0;
  casadi_real w[51];
  const casadi_real* arg[4];
  casadi_real* res[1];
  arg[0] = w+0;
  arg[1] = w+3;
  arg[2] = w+6;
  arg[3] = w+10;
  res[0] = w+11;
  a = w;
  for (j=0; j<11; ++j) if (scanf("%lg", a++)<=0) return 2;
  flag = da_dt(arg, res, iw, w+12, 0);
  if (flag) return flag;
  r = w+11;
  for (j=0; j<1; ++j) CASADI_PRINTF("%g ", *r++);
  CASADI_PRINTF("\n");
  return 0;
}


#ifdef MATLAB_MEX_FILE
void mexFunction(int resc, mxArray *resv[], int argc, const mxArray *argv[]) {
  char buf[6];
  int buf_ok = argc > 0 && !mxGetString(*argv, buf, sizeof(buf));
  if (!buf_ok) {
    mex_da_dt(resc, resv, argc, argv);
    return;
  } else if (strcmp(buf, "da_dt")==0) {
    mex_da_dt(resc, resv, argc-1, argv+1);
    return;
  }
  mexErrMsgTxt("First input should be a command string. Possible values: 'da_dt'");
}
#endif
int main(int argc, char* argv[]) {
  if (argc<2) {
    /* name error */
  } else if (strcmp(argv[1], "da_dt")==0) {
    return main_da_dt(argc-2, argv+2);
  }
  fprintf(stderr, "First input should be a command string. Possible values: 'da_dt'\nNote: you may use function.generate_input to create a command string.\n");
  return 1;
}
#ifdef __cplusplus
} /* extern "C" */
#endif