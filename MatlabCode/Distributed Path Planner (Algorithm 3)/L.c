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
  #define CASADI_PREFIX(ID) L_ ## ID
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

static const casadi_int casadi_s0[40] = {36, 1, 0, 36, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35};
static const casadi_int casadi_s1[5] = {1, 1, 0, 1, 0};
static const casadi_int casadi_s2[10] = {6, 1, 0, 6, 0, 1, 2, 3, 4, 5};

/* L:(X[36],dX[36],k,t,XChief[6])->(CurveLength) */
static int casadi_f0(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem) {
  casadi_real a0, a1, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a2, a20, a21, a22, a23, a24, a25, a26, a27, a28, a29, a3, a30, a31, a32, a33, a34, a35, a36, a37, a38, a39, a4, a40, a41, a42, a43, a44, a45, a5, a6, a7, a8, a9;
  a0=arg[1]? arg[1][0] : 0;
  a1=1.6840054676261389e+04;
  a2=arg[0]? arg[0][3] : 0;
  a3=(a1*a2);
  a4=(a0-a3);
  a5=arg[2]? arg[2][0] : 0;
  a4=(a4*a5);
  a0=(a0-a3);
  a4=(a4*a0);
  a0=arg[1]? arg[1][1] : 0;
  a3=arg[0]? arg[0][4] : 0;
  a6=(a1*a3);
  a7=(a0-a6);
  a7=(a7*a5);
  a0=(a0-a6);
  a7=(a7*a0);
  a4=(a4+a7);
  a7=arg[1]? arg[1][2] : 0;
  a0=arg[0]? arg[0][5] : 0;
  a6=(a1*a0);
  a8=(a7-a6);
  a8=(a8*a5);
  a7=(a7-a6);
  a8=(a8*a7);
  a4=(a4+a8);
  a8=arg[1]? arg[1][3] : 0;
  a7=-3.9860044150000002e+05;
  a6=arg[4]? arg[4][0] : 0;
  a9=casadi_sq(a6);
  a10=arg[4]? arg[4][1] : 0;
  a11=casadi_sq(a10);
  a9=(a9+a11);
  a11=arg[4]? arg[4][2] : 0;
  a12=casadi_sq(a11);
  a9=(a9+a12);
  a12=sqrt(a9);
  a13=arg[0]? arg[0][0] : 0;
  a14=(a12+a13);
  a15=casadi_sq(a14);
  a16=arg[0]? arg[0][1] : 0;
  a17=casadi_sq(a16);
  a15=(a15+a17);
  a17=arg[0]? arg[0][2] : 0;
  a18=casadi_sq(a17);
  a15=(a15+a18);
  a18=1.5000000000000000e+00;
  a15=pow(a15,a18);
  a15=(a7/a15);
  a14=(a15*a14);
  a19=3.9860044150000002e+05;
  a20=(a19/a9);
  a14=(a14+a20);
  a20=-2.;
  a21=arg[4]? arg[4][3] : 0;
  a22=(a11*a21);
  a23=arg[4]? arg[4][5] : 0;
  a24=(a6*a23);
  a22=(a22-a24);
  a24=(a10*a23);
  a25=arg[4]? arg[4][4] : 0;
  a26=(a11*a25);
  a24=(a24-a26);
  a26=casadi_sq(a24);
  a27=casadi_sq(a22);
  a26=(a26+a27);
  a27=(a6*a25);
  a28=(a10*a21);
  a27=(a27-a28);
  a28=casadi_sq(a27);
  a26=(a26+a28);
  a26=sqrt(a26);
  a22=(a22/a26);
  a28=casadi_sq(a6);
  a29=casadi_sq(a10);
  a28=(a28+a29);
  a29=casadi_sq(a11);
  a28=(a28+a29);
  a28=sqrt(a28);
  a29=(a11/a28);
  a30=(a22*a29);
  a27=(a27/a26);
  a31=(a10/a28);
  a32=(a27*a31);
  a30=(a30-a32);
  a32=casadi_sq(a30);
  a28=(a6/a28);
  a33=(a27*a28);
  a24=(a24/a26);
  a26=(a24*a29);
  a33=(a33-a26);
  a26=casadi_sq(a33);
  a32=(a32+a26);
  a26=(a24*a31);
  a34=(a22*a28);
  a26=(a26-a34);
  a34=casadi_sq(a26);
  a32=(a32+a34);
  a32=sqrt(a32);
  a30=(a30/a32);
  a34=(a10*a23);
  a35=(a11*a25);
  a34=(a34-a35);
  a35=(a34/a9);
  a36=(a30*a35);
  a33=(a33/a32);
  a37=(a11*a21);
  a38=(a6*a23);
  a37=(a37-a38);
  a38=(a37/a9);
  a39=(a33*a38);
  a36=(a36+a39);
  a26=(a26/a32);
  a32=(a6*a25);
  a39=(a10*a21);
  a32=(a32-a39);
  a39=(a32/a9);
  a40=(a26*a39);
  a36=(a36+a40);
  a40=(a20*a36);
  a41=(a40*a0);
  a42=(a24*a35);
  a43=(a22*a38);
  a42=(a42+a43);
  a43=(a27*a39);
  a42=(a42+a43);
  a43=(a20*a42);
  a44=(a43*a3);
  a41=(a41-a44);
  a35=(a28*a35);
  a38=(a31*a38);
  a35=(a35+a38);
  a39=(a29*a39);
  a35=(a35+a39);
  a39=(a35*a16);
  a38=(a36*a13);
  a39=(a39-a38);
  a38=(a36*a39);
  a44=(a42*a13);
  a45=(a35*a17);
  a44=(a44-a45);
  a45=(a42*a44);
  a38=(a38-a45);
  a41=(a41-a38);
  a6=(a6*a21);
  a10=(a10*a25);
  a6=(a6+a10);
  a11=(a11*a23);
  a6=(a6+a11);
  a6=(a20*a6);
  a34=(a6*a34);
  a11=casadi_sq(a9);
  a34=(a34/a11);
  a30=(a30*a34);
  a37=(a6*a37);
  a37=(a37/a11);
  a33=(a33*a37);
  a30=(a30+a33);
  a6=(a6*a32);
  a6=(a6/a11);
  a26=(a26*a6);
  a30=(a30+a26);
  a26=(a30*a17);
  a24=(a24*a34);
  a22=(a22*a37);
  a24=(a24+a22);
  a27=(a27*a6);
  a24=(a24+a27);
  a27=(a24*a16);
  a26=(a26-a27);
  a41=(a41-a26);
  a14=(a14+a41);
  a14=(a1*a14);
  a8=(a8-a14);
  a8=casadi_sq(a8);
  a4=(a4+a8);
  a8=arg[1]? arg[1][4] : 0;
  a14=(a15*a16);
  a41=(a20*a42);
  a26=(a41*a2);
  a27=(a20*a35);
  a0=(a27*a0);
  a26=(a26-a0);
  a0=(a36*a17);
  a22=(a42*a16);
  a0=(a0-a22);
  a22=(a42*a0);
  a39=(a35*a39);
  a22=(a22-a39);
  a26=(a26-a22);
  a22=(a24*a13);
  a28=(a28*a34);
  a31=(a31*a37);
  a28=(a28+a31);
  a29=(a29*a6);
  a28=(a28+a29);
  a29=(a28*a17);
  a22=(a22-a29);
  a26=(a26-a22);
  a14=(a14+a26);
  a14=(a1*a14);
  a8=(a8-a14);
  a8=casadi_sq(a8);
  a4=(a4+a8);
  a8=arg[1]? arg[1][5] : 0;
  a15=(a15*a17);
  a17=(a20*a35);
  a3=(a17*a3);
  a20=(a20*a36);
  a2=(a20*a2);
  a3=(a3-a2);
  a44=(a35*a44);
  a0=(a36*a0);
  a44=(a44-a0);
  a3=(a3-a44);
  a16=(a28*a16);
  a13=(a30*a13);
  a16=(a16-a13);
  a3=(a3-a16);
  a15=(a15+a3);
  a15=(a1*a15);
  a8=(a8-a15);
  a8=casadi_sq(a8);
  a4=(a4+a8);
  a8=arg[1]? arg[1][6] : 0;
  a15=arg[0]? arg[0][9] : 0;
  a3=(a1*a15);
  a16=(a8-a3);
  a16=(a16*a5);
  a8=(a8-a3);
  a16=(a16*a8);
  a4=(a4+a16);
  a16=arg[1]? arg[1][7] : 0;
  a8=arg[0]? arg[0][10] : 0;
  a3=(a1*a8);
  a13=(a16-a3);
  a13=(a13*a5);
  a16=(a16-a3);
  a13=(a13*a16);
  a4=(a4+a13);
  a13=arg[1]? arg[1][8] : 0;
  a16=arg[0]? arg[0][11] : 0;
  a3=(a1*a16);
  a44=(a13-a3);
  a44=(a44*a5);
  a13=(a13-a3);
  a44=(a44*a13);
  a4=(a4+a44);
  a44=arg[1]? arg[1][9] : 0;
  a13=arg[0]? arg[0][6] : 0;
  a3=(a12+a13);
  a0=casadi_sq(a3);
  a2=arg[0]? arg[0][7] : 0;
  a14=casadi_sq(a2);
  a0=(a0+a14);
  a14=arg[0]? arg[0][8] : 0;
  a26=casadi_sq(a14);
  a0=(a0+a26);
  a0=pow(a0,a18);
  a0=(a7/a0);
  a3=(a0*a3);
  a26=(a19/a9);
  a3=(a3+a26);
  a26=(a40*a16);
  a22=(a43*a8);
  a26=(a26-a22);
  a22=(a35*a2);
  a29=(a36*a13);
  a22=(a22-a29);
  a29=(a36*a22);
  a6=(a42*a13);
  a31=(a35*a14);
  a6=(a6-a31);
  a31=(a42*a6);
  a29=(a29-a31);
  a26=(a26-a29);
  a29=(a30*a14);
  a31=(a24*a2);
  a29=(a29-a31);
  a26=(a26-a29);
  a3=(a3+a26);
  a3=(a1*a3);
  a44=(a44-a3);
  a44=casadi_sq(a44);
  a4=(a4+a44);
  a44=arg[1]? arg[1][10] : 0;
  a3=(a0*a2);
  a26=(a41*a15);
  a16=(a27*a16);
  a26=(a26-a16);
  a16=(a36*a14);
  a29=(a42*a2);
  a16=(a16-a29);
  a29=(a42*a16);
  a22=(a35*a22);
  a29=(a29-a22);
  a26=(a26-a29);
  a29=(a24*a13);
  a22=(a28*a14);
  a29=(a29-a22);
  a26=(a26-a29);
  a3=(a3+a26);
  a3=(a1*a3);
  a44=(a44-a3);
  a44=casadi_sq(a44);
  a4=(a4+a44);
  a44=arg[1]? arg[1][11] : 0;
  a0=(a0*a14);
  a8=(a17*a8);
  a15=(a20*a15);
  a8=(a8-a15);
  a6=(a35*a6);
  a16=(a36*a16);
  a6=(a6-a16);
  a8=(a8-a6);
  a2=(a28*a2);
  a13=(a30*a13);
  a2=(a2-a13);
  a8=(a8-a2);
  a0=(a0+a8);
  a0=(a1*a0);
  a44=(a44-a0);
  a44=casadi_sq(a44);
  a4=(a4+a44);
  a44=arg[1]? arg[1][12] : 0;
  a0=arg[0]? arg[0][15] : 0;
  a8=(a1*a0);
  a2=(a44-a8);
  a2=(a2*a5);
  a44=(a44-a8);
  a2=(a2*a44);
  a4=(a4+a2);
  a2=arg[1]? arg[1][13] : 0;
  a44=arg[0]? arg[0][16] : 0;
  a8=(a1*a44);
  a13=(a2-a8);
  a13=(a13*a5);
  a2=(a2-a8);
  a13=(a13*a2);
  a4=(a4+a13);
  a13=arg[1]? arg[1][14] : 0;
  a2=arg[0]? arg[0][17] : 0;
  a8=(a1*a2);
  a6=(a13-a8);
  a6=(a6*a5);
  a13=(a13-a8);
  a6=(a6*a13);
  a4=(a4+a6);
  a6=arg[1]? arg[1][15] : 0;
  a13=arg[0]? arg[0][12] : 0;
  a8=(a12+a13);
  a16=casadi_sq(a8);
  a15=arg[0]? arg[0][13] : 0;
  a14=casadi_sq(a15);
  a16=(a16+a14);
  a14=arg[0]? arg[0][14] : 0;
  a3=casadi_sq(a14);
  a16=(a16+a3);
  a16=pow(a16,a18);
  a16=(a7/a16);
  a8=(a16*a8);
  a3=(a19/a9);
  a8=(a8+a3);
  a3=(a40*a2);
  a26=(a43*a44);
  a3=(a3-a26);
  a26=(a35*a15);
  a29=(a36*a13);
  a26=(a26-a29);
  a29=(a36*a26);
  a22=(a42*a13);
  a31=(a35*a14);
  a22=(a22-a31);
  a31=(a42*a22);
  a29=(a29-a31);
  a3=(a3-a29);
  a29=(a30*a14);
  a31=(a24*a15);
  a29=(a29-a31);
  a3=(a3-a29);
  a8=(a8+a3);
  a8=(a1*a8);
  a6=(a6-a8);
  a6=casadi_sq(a6);
  a4=(a4+a6);
  a6=arg[1]? arg[1][16] : 0;
  a8=(a16*a15);
  a3=(a41*a0);
  a2=(a27*a2);
  a3=(a3-a2);
  a2=(a36*a14);
  a29=(a42*a15);
  a2=(a2-a29);
  a29=(a42*a2);
  a26=(a35*a26);
  a29=(a29-a26);
  a3=(a3-a29);
  a29=(a24*a13);
  a26=(a28*a14);
  a29=(a29-a26);
  a3=(a3-a29);
  a8=(a8+a3);
  a8=(a1*a8);
  a6=(a6-a8);
  a6=casadi_sq(a6);
  a4=(a4+a6);
  a6=arg[1]? arg[1][17] : 0;
  a16=(a16*a14);
  a44=(a17*a44);
  a0=(a20*a0);
  a44=(a44-a0);
  a22=(a35*a22);
  a2=(a36*a2);
  a22=(a22-a2);
  a44=(a44-a22);
  a15=(a28*a15);
  a13=(a30*a13);
  a15=(a15-a13);
  a44=(a44-a15);
  a16=(a16+a44);
  a16=(a1*a16);
  a6=(a6-a16);
  a6=casadi_sq(a6);
  a4=(a4+a6);
  a6=arg[1]? arg[1][18] : 0;
  a16=arg[0]? arg[0][21] : 0;
  a44=(a1*a16);
  a15=(a6-a44);
  a15=(a15*a5);
  a6=(a6-a44);
  a15=(a15*a6);
  a4=(a4+a15);
  a15=arg[1]? arg[1][19] : 0;
  a6=arg[0]? arg[0][22] : 0;
  a44=(a1*a6);
  a13=(a15-a44);
  a13=(a13*a5);
  a15=(a15-a44);
  a13=(a13*a15);
  a4=(a4+a13);
  a13=arg[1]? arg[1][20] : 0;
  a15=arg[0]? arg[0][23] : 0;
  a44=(a1*a15);
  a22=(a13-a44);
  a22=(a22*a5);
  a13=(a13-a44);
  a22=(a22*a13);
  a4=(a4+a22);
  a22=arg[1]? arg[1][21] : 0;
  a13=arg[0]? arg[0][18] : 0;
  a44=(a12+a13);
  a2=casadi_sq(a44);
  a0=arg[0]? arg[0][19] : 0;
  a14=casadi_sq(a0);
  a2=(a2+a14);
  a14=arg[0]? arg[0][20] : 0;
  a8=casadi_sq(a14);
  a2=(a2+a8);
  a2=pow(a2,a18);
  a2=(a7/a2);
  a44=(a2*a44);
  a8=(a19/a9);
  a44=(a44+a8);
  a8=(a40*a15);
  a3=(a43*a6);
  a8=(a8-a3);
  a3=(a35*a0);
  a29=(a36*a13);
  a3=(a3-a29);
  a29=(a36*a3);
  a26=(a42*a13);
  a31=(a35*a14);
  a26=(a26-a31);
  a31=(a42*a26);
  a29=(a29-a31);
  a8=(a8-a29);
  a29=(a30*a14);
  a31=(a24*a0);
  a29=(a29-a31);
  a8=(a8-a29);
  a44=(a44+a8);
  a44=(a1*a44);
  a22=(a22-a44);
  a22=casadi_sq(a22);
  a4=(a4+a22);
  a22=arg[1]? arg[1][22] : 0;
  a44=(a2*a0);
  a8=(a41*a16);
  a15=(a27*a15);
  a8=(a8-a15);
  a15=(a36*a14);
  a29=(a42*a0);
  a15=(a15-a29);
  a29=(a42*a15);
  a3=(a35*a3);
  a29=(a29-a3);
  a8=(a8-a29);
  a29=(a24*a13);
  a3=(a28*a14);
  a29=(a29-a3);
  a8=(a8-a29);
  a44=(a44+a8);
  a44=(a1*a44);
  a22=(a22-a44);
  a22=casadi_sq(a22);
  a4=(a4+a22);
  a22=arg[1]? arg[1][23] : 0;
  a2=(a2*a14);
  a6=(a17*a6);
  a16=(a20*a16);
  a6=(a6-a16);
  a26=(a35*a26);
  a15=(a36*a15);
  a26=(a26-a15);
  a6=(a6-a26);
  a0=(a28*a0);
  a13=(a30*a13);
  a0=(a0-a13);
  a6=(a6-a0);
  a2=(a2+a6);
  a2=(a1*a2);
  a22=(a22-a2);
  a22=casadi_sq(a22);
  a4=(a4+a22);
  a22=arg[1]? arg[1][24] : 0;
  a2=arg[0]? arg[0][27] : 0;
  a6=(a1*a2);
  a0=(a22-a6);
  a0=(a0*a5);
  a22=(a22-a6);
  a0=(a0*a22);
  a4=(a4+a0);
  a0=arg[1]? arg[1][25] : 0;
  a22=arg[0]? arg[0][28] : 0;
  a6=(a1*a22);
  a13=(a0-a6);
  a13=(a13*a5);
  a0=(a0-a6);
  a13=(a13*a0);
  a4=(a4+a13);
  a13=arg[1]? arg[1][26] : 0;
  a0=arg[0]? arg[0][29] : 0;
  a6=(a1*a0);
  a26=(a13-a6);
  a26=(a26*a5);
  a13=(a13-a6);
  a26=(a26*a13);
  a4=(a4+a26);
  a26=arg[1]? arg[1][27] : 0;
  a13=arg[0]? arg[0][24] : 0;
  a6=(a12+a13);
  a15=casadi_sq(a6);
  a16=arg[0]? arg[0][25] : 0;
  a14=casadi_sq(a16);
  a15=(a15+a14);
  a14=arg[0]? arg[0][26] : 0;
  a44=casadi_sq(a14);
  a15=(a15+a44);
  a15=pow(a15,a18);
  a15=(a7/a15);
  a6=(a15*a6);
  a44=(a19/a9);
  a6=(a6+a44);
  a44=(a40*a0);
  a8=(a43*a22);
  a44=(a44-a8);
  a8=(a35*a16);
  a29=(a36*a13);
  a8=(a8-a29);
  a29=(a36*a8);
  a3=(a42*a13);
  a31=(a35*a14);
  a3=(a3-a31);
  a31=(a42*a3);
  a29=(a29-a31);
  a44=(a44-a29);
  a29=(a30*a14);
  a31=(a24*a16);
  a29=(a29-a31);
  a44=(a44-a29);
  a6=(a6+a44);
  a6=(a1*a6);
  a26=(a26-a6);
  a26=casadi_sq(a26);
  a4=(a4+a26);
  a26=arg[1]? arg[1][28] : 0;
  a6=(a15*a16);
  a44=(a41*a2);
  a0=(a27*a0);
  a44=(a44-a0);
  a0=(a36*a14);
  a29=(a42*a16);
  a0=(a0-a29);
  a29=(a42*a0);
  a8=(a35*a8);
  a29=(a29-a8);
  a44=(a44-a29);
  a29=(a24*a13);
  a8=(a28*a14);
  a29=(a29-a8);
  a44=(a44-a29);
  a6=(a6+a44);
  a6=(a1*a6);
  a26=(a26-a6);
  a26=casadi_sq(a26);
  a4=(a4+a26);
  a26=arg[1]? arg[1][29] : 0;
  a15=(a15*a14);
  a22=(a17*a22);
  a2=(a20*a2);
  a22=(a22-a2);
  a3=(a35*a3);
  a0=(a36*a0);
  a3=(a3-a0);
  a22=(a22-a3);
  a16=(a28*a16);
  a13=(a30*a13);
  a16=(a16-a13);
  a22=(a22-a16);
  a15=(a15+a22);
  a15=(a1*a15);
  a26=(a26-a15);
  a26=casadi_sq(a26);
  a4=(a4+a26);
  a26=arg[1]? arg[1][30] : 0;
  a15=arg[0]? arg[0][33] : 0;
  a22=(a1*a15);
  a16=(a26-a22);
  a16=(a16*a5);
  a26=(a26-a22);
  a16=(a16*a26);
  a4=(a4+a16);
  a16=arg[1]? arg[1][31] : 0;
  a26=arg[0]? arg[0][34] : 0;
  a22=(a1*a26);
  a13=(a16-a22);
  a13=(a13*a5);
  a16=(a16-a22);
  a13=(a13*a16);
  a4=(a4+a13);
  a13=arg[1]? arg[1][32] : 0;
  a16=arg[0]? arg[0][35] : 0;
  a22=(a1*a16);
  a3=(a13-a22);
  a3=(a3*a5);
  a13=(a13-a22);
  a3=(a3*a13);
  a4=(a4+a3);
  a3=arg[1]? arg[1][33] : 0;
  a13=arg[0]? arg[0][30] : 0;
  a12=(a12+a13);
  a22=casadi_sq(a12);
  a5=arg[0]? arg[0][31] : 0;
  a0=casadi_sq(a5);
  a22=(a22+a0);
  a0=arg[0]? arg[0][32] : 0;
  a2=casadi_sq(a0);
  a22=(a22+a2);
  a22=pow(a22,a18);
  a7=(a7/a22);
  a12=(a7*a12);
  a19=(a19/a9);
  a12=(a12+a19);
  a40=(a40*a16);
  a43=(a43*a26);
  a40=(a40-a43);
  a43=(a35*a5);
  a19=(a36*a13);
  a43=(a43-a19);
  a19=(a36*a43);
  a9=(a42*a13);
  a22=(a35*a0);
  a9=(a9-a22);
  a22=(a42*a9);
  a19=(a19-a22);
  a40=(a40-a19);
  a19=(a30*a0);
  a22=(a24*a5);
  a19=(a19-a22);
  a40=(a40-a19);
  a12=(a12+a40);
  a12=(a1*a12);
  a3=(a3-a12);
  a3=casadi_sq(a3);
  a4=(a4+a3);
  a3=arg[1]? arg[1][34] : 0;
  a12=(a7*a5);
  a41=(a41*a15);
  a27=(a27*a16);
  a41=(a41-a27);
  a27=(a36*a0);
  a16=(a42*a5);
  a27=(a27-a16);
  a42=(a42*a27);
  a43=(a35*a43);
  a42=(a42-a43);
  a41=(a41-a42);
  a24=(a24*a13);
  a42=(a28*a0);
  a24=(a24-a42);
  a41=(a41-a24);
  a12=(a12+a41);
  a12=(a1*a12);
  a3=(a3-a12);
  a3=casadi_sq(a3);
  a4=(a4+a3);
  a3=arg[1]? arg[1][35] : 0;
  a7=(a7*a0);
  a17=(a17*a26);
  a20=(a20*a15);
  a17=(a17-a20);
  a35=(a35*a9);
  a36=(a36*a27);
  a35=(a35-a36);
  a17=(a17-a35);
  a28=(a28*a5);
  a30=(a30*a13);
  a28=(a28-a30);
  a17=(a17-a28);
  a7=(a7+a17);
  a1=(a1*a7);
  a3=(a3-a1);
  a3=casadi_sq(a3);
  a4=(a4+a3);
  if (res[0]!=0) res[0][0]=a4;
  return 0;
}

CASADI_SYMBOL_EXPORT int L(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem){
  return casadi_f0(arg, res, iw, w, mem);
}

CASADI_SYMBOL_EXPORT int L_alloc_mem(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT int L_init_mem(int mem) {
  return 0;
}

CASADI_SYMBOL_EXPORT void L_free_mem(int mem) {
}

CASADI_SYMBOL_EXPORT int L_checkout(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT void L_release(int mem) {
}

CASADI_SYMBOL_EXPORT void L_incref(void) {
}

CASADI_SYMBOL_EXPORT void L_decref(void) {
}

CASADI_SYMBOL_EXPORT casadi_int L_n_in(void) { return 5;}

CASADI_SYMBOL_EXPORT casadi_int L_n_out(void) { return 1;}

CASADI_SYMBOL_EXPORT casadi_real L_default_in(casadi_int i){
  switch (i) {
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* L_name_in(casadi_int i){
  switch (i) {
    case 0: return "X";
    case 1: return "dX";
    case 2: return "k";
    case 3: return "t";
    case 4: return "XChief";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* L_name_out(casadi_int i){
  switch (i) {
    case 0: return "CurveLength";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* L_sparsity_in(casadi_int i) {
  switch (i) {
    case 0: return casadi_s0;
    case 1: return casadi_s0;
    case 2: return casadi_s1;
    case 3: return casadi_s1;
    case 4: return casadi_s2;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* L_sparsity_out(casadi_int i) {
  switch (i) {
    case 0: return casadi_s1;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT int L_work(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w) {
  if (sz_arg) *sz_arg = 5;
  if (sz_res) *sz_res = 1;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 0;
  return 0;
}

#ifdef MATLAB_MEX_FILE
void mex_L(int resc, mxArray *resv[], int argc, const mxArray *argv[]) {
  casadi_int i;
  casadi_real w[127];
  casadi_int *iw = 0;
  const casadi_real* arg[5] = {0};
  casadi_real* res[1] = {0};
  if (argc>5) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"L\" failed. Too many input arguments (%d, max 5)", argc);
  if (resc>1) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"L\" failed. Too many output arguments (%d, max 1)", resc);
  if (--argc>=0) arg[0] = casadi_from_mex(argv[0], w, casadi_s0, w+81);
  if (--argc>=0) arg[1] = casadi_from_mex(argv[1], w+36, casadi_s0, w+81);
  if (--argc>=0) arg[2] = casadi_from_mex(argv[2], w+72, casadi_s1, w+81);
  if (--argc>=0) arg[3] = casadi_from_mex(argv[3], w+73, casadi_s1, w+81);
  if (--argc>=0) arg[4] = casadi_from_mex(argv[4], w+74, casadi_s2, w+81);
  --resc;
  res[0] = w+80;
  i = L(arg, res, iw, w+81, 0);
  if (i) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"L\" failed.");
  if (res[0]) resv[0] = casadi_to_mex(casadi_s1, res[0]);
}
#endif

casadi_int main_L(casadi_int argc, char* argv[]) {
  casadi_int j;
  casadi_real* a;
  const casadi_real* r;
  casadi_int flag;
  casadi_int *iw = 0;
  casadi_real w[127];
  const casadi_real* arg[5];
  casadi_real* res[1];
  arg[0] = w+0;
  arg[1] = w+36;
  arg[2] = w+72;
  arg[3] = w+73;
  arg[4] = w+74;
  res[0] = w+80;
  a = w;
  for (j=0; j<80; ++j) if (scanf("%lg", a++)<=0) return 2;
  flag = L(arg, res, iw, w+81, 0);
  if (flag) return flag;
  r = w+80;
  for (j=0; j<1; ++j) CASADI_PRINTF("%g ", *r++);
  CASADI_PRINTF("\n");
  return 0;
}


#ifdef MATLAB_MEX_FILE
void mexFunction(int resc, mxArray *resv[], int argc, const mxArray *argv[]) {
  char buf[2];
  int buf_ok = argc > 0 && !mxGetString(*argv, buf, sizeof(buf));
  if (!buf_ok) {
    mex_L(resc, resv, argc, argv);
    return;
  } else if (strcmp(buf, "L")==0) {
    mex_L(resc, resv, argc-1, argv+1);
    return;
  }
  mexErrMsgTxt("First input should be a command string. Possible values: 'L'");
}
#endif
int main(int argc, char* argv[]) {
  if (argc<2) {
    /* name error */
  } else if (strcmp(argv[1], "L")==0) {
    return main_L(argc-2, argv+2);
  }
  fprintf(stderr, "First input should be a command string. Possible values: 'L'\nNote: you may use function.generate_input to create a command string.\n");
  return 1;
}
#ifdef __cplusplus
} /* extern "C" */
#endif
