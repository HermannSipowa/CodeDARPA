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

static const casadi_int casadi_s0[17] = {13, 1, 0, 13, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
static const casadi_int casadi_s1[5] = {1, 1, 0, 1, 0};
static const casadi_int casadi_s2[10] = {6, 1, 0, 6, 0, 1, 2, 3, 4, 5};

/* L:(X[13],dX[13],k1,k2,k3,t,XChief[6])->(CurveLength) */
static int casadi_f0(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem) {
  casadi_real a0, a1, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a2, a20, a21, a22, a23, a24, a25, a26, a27, a28, a29, a3, a30, a31, a32, a33, a34, a35, a36, a37, a38, a39, a4, a40, a41, a42, a43, a44, a45, a46, a47, a48, a49, a5, a50, a51, a52, a53, a54, a55, a56, a57, a58, a59, a6, a60, a61, a62, a63, a64, a65, a66, a67, a68, a69, a7, a70, a71, a72, a73, a74, a75, a76, a77, a78, a79, a8, a80, a81, a82, a83, a84, a85, a86, a87, a88, a89, a9, a90, a91, a92, a93, a94, a95;
  a0=arg[1]? arg[1][0] : 0;
  a1=1.1259418598742898e+05;
  a2=5.0000000000000000e-01;
  a3=arg[0]? arg[0][4] : 0;
  a4=(a2*a3);
  a5=arg[0]? arg[0][1] : 0;
  a6=arg[0]? arg[0][0] : 0;
  a7=casadi_sq(a6);
  a8=casadi_sq(a5);
  a7=(a7+a8);
  a8=arg[0]? arg[0][2] : 0;
  a9=casadi_sq(a8);
  a7=(a7+a9);
  a9=arg[0]? arg[0][3] : 0;
  a10=casadi_sq(a9);
  a7=(a7+a10);
  a7=sqrt(a7);
  a5=(a5/a7);
  a4=(a4*a5);
  a10=arg[0]? arg[0][5] : 0;
  a11=(a2*a10);
  a8=(a8/a7);
  a11=(a11*a8);
  a4=(a4+a11);
  a11=arg[0]? arg[0][6] : 0;
  a12=(a2*a11);
  a9=(a9/a7);
  a12=(a12*a9);
  a4=(a4+a12);
  a4=(a1*a4);
  a12=(a0+a4);
  a13=arg[2]? arg[2][0] : 0;
  a12=(a12*a13);
  a0=(a0+a4);
  a12=(a12*a0);
  a0=arg[1]? arg[1][1] : 0;
  a4=(a2*a3);
  a6=(a6/a7);
  a4=(a4*a6);
  a7=(a2*a11);
  a7=(a7*a8);
  a4=(a4+a7);
  a7=(a2*a10);
  a7=(a7*a9);
  a4=(a4-a7);
  a4=(a1*a4);
  a7=(a0-a4);
  a7=(a7*a13);
  a0=(a0-a4);
  a7=(a7*a0);
  a12=(a12+a7);
  a7=arg[1]? arg[1][2] : 0;
  a0=(a2*a10);
  a0=(a0*a6);
  a4=(a2*a11);
  a4=(a4*a5);
  a0=(a0-a4);
  a4=(a2*a3);
  a4=(a4*a9);
  a0=(a0+a4);
  a0=(a1*a0);
  a4=(a7-a0);
  a4=(a4*a13);
  a7=(a7-a0);
  a4=(a4*a7);
  a12=(a12+a4);
  a4=arg[1]? arg[1][3] : 0;
  a7=(a2*a11);
  a7=(a7*a6);
  a0=(a2*a10);
  a0=(a0*a5);
  a7=(a7+a0);
  a0=(a2*a3);
  a0=(a0*a8);
  a7=(a7-a0);
  a7=(a1*a7);
  a0=(a4-a7);
  a0=(a0*a13);
  a4=(a4-a7);
  a0=(a0*a4);
  a12=(a12+a0);
  a0=1.4692378328741967e-02;
  a4=arg[1]? arg[1][4] : 0;
  a7=1.2121212121212120e-01;
  a13=(a7*a10);
  a14=(a11*a13);
  a15=2.4242424242424240e-01;
  a16=(a15*a11);
  a17=(a10*a16);
  a14=(a14-a17);
  a14=(a14/a7);
  a14=(a1*a14);
  a17=(a4-a14);
  a17=(a0*a17);
  a4=(a4-a14);
  a17=(a17*a4);
  a12=(a12+a17);
  a17=arg[1]? arg[1][5] : 0;
  a16=(a3*a16);
  a4=(a7*a3);
  a11=(a11*a4);
  a16=(a16-a11);
  a16=(a16/a7);
  a16=(a1*a16);
  a7=(a17-a16);
  a0=(a0*a7);
  a17=(a17-a16);
  a0=(a0*a17);
  a12=(a12+a0);
  a0=5.8769513314967867e-02;
  a17=arg[1]? arg[1][6] : 0;
  a10=(a10*a4);
  a3=(a3*a13);
  a10=(a10-a3);
  a10=(a10/a15);
  a10=(a1*a10);
  a15=(a17-a10);
  a0=(a0*a15);
  a17=(a17-a10);
  a0=(a0*a17);
  a12=(a12+a0);
  a0=arg[1]? arg[1][7] : 0;
  a17=arg[0]? arg[0][10] : 0;
  a10=(a1*a17);
  a15=(a0-a10);
  a3=arg[3]? arg[3][0] : 0;
  a15=(a15*a3);
  a0=(a0-a10);
  a15=(a15*a0);
  a12=(a12+a15);
  a15=arg[1]? arg[1][8] : 0;
  a0=arg[0]? arg[0][11] : 0;
  a10=(a1*a0);
  a13=(a15-a10);
  a13=(a13*a3);
  a15=(a15-a10);
  a13=(a13*a15);
  a12=(a12+a13);
  a13=arg[1]? arg[1][9] : 0;
  a15=arg[0]? arg[0][12] : 0;
  a10=(a1*a15);
  a4=(a13-a10);
  a4=(a4*a3);
  a13=(a13-a10);
  a4=(a4*a13);
  a12=(a12+a4);
  a4=arg[1]? arg[1][10] : 0;
  a13=arg[6]? arg[6][0] : 0;
  a10=casadi_sq(a13);
  a3=arg[6]? arg[6][1] : 0;
  a16=casadi_sq(a3);
  a10=(a10+a16);
  a16=arg[6]? arg[6][2] : 0;
  a7=casadi_sq(a16);
  a10=(a10+a7);
  a10=sqrt(a10);
  a7=(a13/a10);
  a11=2.;
  a14=(a5*a8);
  a18=(a6*a9);
  a14=(a14+a18);
  a14=(a11*a14);
  a18=casadi_sq(a6);
  a19=casadi_sq(a5);
  a18=(a18-a19);
  a19=casadi_sq(a8);
  a18=(a18-a19);
  a19=casadi_sq(a9);
  a18=(a18+a19);
  a19=(a14*a18);
  a20=(a5*a9);
  a21=(a6*a8);
  a20=(a20-a21);
  a20=(a11*a20);
  a21=(a8*a9);
  a22=(a6*a5);
  a21=(a21-a22);
  a21=(a11*a21);
  a22=(a20*a21);
  a19=(a19-a22);
  a22=casadi_sq(a6);
  a23=casadi_sq(a5);
  a22=(a22+a23);
  a23=casadi_sq(a8);
  a22=(a22-a23);
  a23=casadi_sq(a9);
  a22=(a22-a23);
  a23=casadi_sq(a6);
  a24=casadi_sq(a5);
  a23=(a23-a24);
  a24=casadi_sq(a8);
  a23=(a23+a24);
  a24=casadi_sq(a9);
  a23=(a23-a24);
  a24=(a23*a18);
  a25=(a8*a9);
  a26=(a6*a5);
  a25=(a25+a26);
  a25=(a11*a25);
  a26=(a25*a21);
  a24=(a24-a26);
  a24=(a22*a24);
  a26=(a5*a8);
  a27=(a6*a9);
  a26=(a26-a27);
  a26=(a11*a26);
  a27=(a26*a18);
  a28=(a5*a9);
  a29=(a6*a8);
  a28=(a28+a29);
  a28=(a11*a28);
  a29=(a25*a28);
  a27=(a27-a29);
  a27=(a14*a27);
  a24=(a24-a27);
  a27=(a26*a21);
  a29=(a23*a28);
  a27=(a27-a29);
  a27=(a20*a27);
  a24=(a24+a27);
  a19=(a19/a24);
  a27=1.0000000000000000e-03;
  a29=5.4317577262067082e-03;
  a30=4.;
  a31=1.7453292519943295e-02;
  a32=1.7487317400000001e+02;
  a33=-2.4109079999999999e-01;
  a34=2.4562962500000000e+06;
  a35=arg[5]? arg[5][0] : 0;
  a35=(a1*a35);
  a36=86400.;
  a37=(a35/a36);
  a37=(a34+a37);
  a38=2451545.;
  a37=(a37-a38);
  a39=36525.;
  a37=(a37/a39);
  a40=(a33*a37);
  a40=(a32+a40);
  a41=4.0670000000000002e-05;
  a42=casadi_sq(a37);
  a43=(a41*a42);
  a40=(a40+a43);
  a43=-1.3270000000000000e-06;
  a44=casadi_sq(a37);
  a44=(a44*a37);
  a45=(a43*a44);
  a40=(a40+a45);
  a40=(a31*a40);
  a45=(-a40);
  a45=cos(a45);
  a46=1.0293734800000000e+02;
  a47=3.2255570000000000e-01;
  a48=(a47*a37);
  a48=(a46+a48);
  a49=1.5025999999999999e-04;
  a50=(a49*a42);
  a48=(a48+a50);
  a50=4.7800000000000002e-07;
  a51=(a50*a44);
  a48=(a48+a51);
  a48=(a31*a48);
  a51=(a48-a40);
  a52=(-a51);
  a52=cos(a52);
  a53=(a45*a52);
  a54=(-a40);
  a54=sin(a54);
  a55=1.3054600000000000e-02;
  a56=(a55*a37);
  a57=-9.3100000000000006e-06;
  a58=(a57*a42);
  a56=(a56+a58);
  a58=-3.4000000000000000e-08;
  a59=(a58*a44);
  a56=(a56+a59);
  a56=(a31*a56);
  a59=(-a56);
  a59=cos(a59);
  a54=(a54*a59);
  a60=(-a51);
  a60=sin(a60);
  a61=(a54*a60);
  a53=(a53-a61);
  a61=1.4959802299063239e+08;
  a62=1.;
  a63=1.6708620000000000e-02;
  a64=-4.2036999999999997e-05;
  a65=(a64*a37);
  a65=(a63+a65);
  a66=-1.2360000000000000e-07;
  a67=(a66*a42);
  a65=(a65+a67);
  a67=3.9999999999999998e-11;
  a44=(a67*a44);
  a65=(a65+a44);
  a44=casadi_sq(a65);
  a44=(a62-a44);
  a44=(a61*a44);
  a68=1.0046644900000000e+02;
  a69=3.5999372851900000e+04;
  a37=(a69*a37);
  a37=(a68+a37);
  a70=-5.6799999999999998e-06;
  a42=(a70*a42);
  a37=(a37+a42);
  a37=(a31*a37);
  a37=(a37-a48);
  a48=(a11*a65);
  a42=casadi_sq(a65);
  a42=(a65*a42);
  a42=(a42/a30);
  a48=(a48-a42);
  a42=5.2083333333333336e-02;
  a71=casadi_sq(a65);
  a71=casadi_sq(a71);
  a71=(a65*a71);
  a71=(a42*a71);
  a48=(a48+a71);
  a71=sin(a37);
  a48=(a48*a71);
  a71=1.2500000000000000e+00;
  a72=casadi_sq(a65);
  a72=(a71*a72);
  a73=4.5833333333333331e-01;
  a74=casadi_sq(a65);
  a74=casadi_sq(a74);
  a74=(a73*a74);
  a72=(a72-a74);
  a74=(a11*a37);
  a74=sin(a74);
  a72=(a72*a74);
  a48=(a48+a72);
  a72=1.0833333333333333e+00;
  a74=casadi_sq(a65);
  a74=(a65*a74);
  a74=(a72*a74);
  a75=6.7187500000000000e-01;
  a76=casadi_sq(a65);
  a76=casadi_sq(a76);
  a76=(a65*a76);
  a76=(a75*a76);
  a74=(a74-a76);
  a76=3.;
  a77=(a76*a37);
  a77=sin(a77);
  a74=(a74*a77);
  a48=(a48+a74);
  a74=1.0729166666666667e+00;
  a77=casadi_sq(a65);
  a77=casadi_sq(a77);
  a77=(a74*a77);
  a78=(a30*a37);
  a78=sin(a78);
  a77=(a77*a78);
  a48=(a48+a77);
  a77=1.1427083333333334e+00;
  a78=casadi_sq(a65);
  a78=casadi_sq(a78);
  a78=(a65*a78);
  a78=(a77*a78);
  a79=5.;
  a80=(a79*a37);
  a80=sin(a80);
  a78=(a78*a80);
  a48=(a48+a78);
  a37=(a37+a48);
  a48=cos(a37);
  a65=(a65*a48);
  a65=(a62+a65);
  a44=(a44/a65);
  a65=cos(a37);
  a65=(a44*a65);
  a53=(a53*a65);
  a48=(-a51);
  a48=sin(a48);
  a45=(a45*a48);
  a51=(-a51);
  a51=cos(a51);
  a54=(a54*a51);
  a45=(a45+a54);
  a37=sin(a37);
  a44=(a44*a37);
  a45=(a45*a44);
  a53=(a53+a45);
  a45=arg[6]? arg[6][4] : 0;
  a37=(a13*a45);
  a54=arg[6]? arg[6][3] : 0;
  a78=(a3*a54);
  a37=(a37-a78);
  a78=arg[6]? arg[6][5] : 0;
  a80=(a3*a78);
  a81=(a16*a45);
  a80=(a80-a81);
  a81=casadi_sq(a80);
  a82=(a16*a54);
  a83=(a13*a78);
  a82=(a82-a83);
  a83=casadi_sq(a82);
  a81=(a81+a83);
  a83=casadi_sq(a37);
  a81=(a81+a83);
  a81=sqrt(a81);
  a37=(a37/a81);
  a83=(a37*a7);
  a80=(a80/a81);
  a84=(a16/a10);
  a85=(a80*a84);
  a83=(a83-a85);
  a82=(a82/a81);
  a81=(a82*a84);
  a10=(a3/a10);
  a85=(a37*a10);
  a81=(a81-a85);
  a85=casadi_sq(a81);
  a86=casadi_sq(a83);
  a85=(a85+a86);
  a86=(a80*a10);
  a87=(a82*a7);
  a86=(a86-a87);
  a87=casadi_sq(a86);
  a85=(a85+a87);
  a85=sqrt(a85);
  a83=(a83/a85);
  a87=(a83*a37);
  a86=(a86/a85);
  a88=(a86*a82);
  a87=(a87-a88);
  a88=(a83*a37);
  a89=(a86*a82);
  a88=(a88-a89);
  a88=(a7*a88);
  a81=(a81/a85);
  a85=(a81*a37);
  a89=(a86*a80);
  a85=(a85-a89);
  a85=(a10*a85);
  a88=(a88-a85);
  a85=(a81*a82);
  a89=(a83*a80);
  a85=(a85-a89);
  a85=(a84*a85);
  a88=(a88+a85);
  a87=(a87/a88);
  a85=casadi_sq(a13);
  a89=casadi_sq(a3);
  a85=(a85+a89);
  a89=casadi_sq(a16);
  a85=(a85+a89);
  a89=sqrt(a85);
  a90=arg[0]? arg[0][7] : 0;
  a89=(a89+a90);
  a87=(a87*a89);
  a91=(a10*a37);
  a92=(a84*a82);
  a91=(a91-a92);
  a91=(a91/a88);
  a92=arg[0]? arg[0][8] : 0;
  a91=(a91*a92);
  a87=(a87-a91);
  a91=(a10*a86);
  a93=(a84*a83);
  a91=(a91-a93);
  a91=(a91/a88);
  a93=arg[0]? arg[0][9] : 0;
  a91=(a91*a93);
  a87=(a87+a91);
  a53=(a53+a87);
  a87=1.4959787069999999e+08;
  a91=(a53/a87);
  a91=casadi_sq(a91);
  a94=9.1748200035787253e-01;
  a95=(-a40);
  a95=cos(a95);
  a95=(a95*a59);
  a59=(a95*a51);
  a40=(-a40);
  a40=sin(a40);
  a48=(a40*a48);
  a59=(a59-a48);
  a59=(a59*a44);
  a40=(a40*a52);
  a95=(a95*a60);
  a40=(a40+a95);
  a40=(a40*a65);
  a59=(a59-a40);
  a40=(a94*a59);
  a95=-3.9777729827042280e-01;
  a56=(-a56);
  a56=sin(a56);
  a60=(a56*a60);
  a60=(a60*a65);
  a56=(a56*a51);
  a56=(a56*a44);
  a60=(a60-a56);
  a56=(a95*a60);
  a40=(a40+a56);
  a56=(a7*a37);
  a44=(a84*a80);
  a56=(a56-a44);
  a56=(a56/a88);
  a56=(a56*a92);
  a44=(a81*a37);
  a51=(a86*a80);
  a44=(a44-a51);
  a44=(a44/a88);
  a44=(a44*a89);
  a56=(a56-a44);
  a44=(a7*a86);
  a51=(a84*a81);
  a44=(a44-a51);
  a44=(a44/a88);
  a44=(a44*a93);
  a56=(a56-a44);
  a40=(a40+a56);
  a56=(a40/a87);
  a56=casadi_sq(a56);
  a91=(a91+a56);
  a56=3.9777729827042280e-01;
  a59=(a56*a59);
  a60=(a94*a60);
  a59=(a59+a60);
  a60=(a81*a82);
  a44=(a83*a80);
  a60=(a60-a44);
  a60=(a60/a88);
  a60=(a60*a89);
  a44=(a7*a82);
  a51=(a10*a80);
  a44=(a44-a51);
  a44=(a44/a88);
  a44=(a44*a92);
  a60=(a60-a44);
  a44=(a7*a83);
  a51=(a10*a81);
  a44=(a44-a51);
  a44=(a44/a88);
  a44=(a44*a93);
  a60=(a60+a44);
  a59=(a59+a60);
  a87=(a59/a87);
  a87=casadi_sq(a87);
  a91=(a91+a87);
  a91=(a30*a91);
  a29=(a29/a91);
  a27=(a27*a29);
  a29=(a5*a9);
  a91=(a6*a8);
  a29=(a29+a91);
  a29=(a11*a29);
  a91=casadi_sq(a53);
  a87=casadi_sq(a40);
  a91=(a91+a87);
  a87=casadi_sq(a59);
  a91=(a91+a87);
  a91=sqrt(a91);
  a53=(a53/a91);
  a29=(a29*a53);
  a87=(a8*a9);
  a60=(a6*a5);
  a87=(a87-a60);
  a87=(a11*a87);
  a40=(a40/a91);
  a87=(a87*a40);
  a29=(a29+a87);
  a87=casadi_sq(a6);
  a60=casadi_sq(a5);
  a87=(a87-a60);
  a60=casadi_sq(a8);
  a87=(a87-a60);
  a60=casadi_sq(a9);
  a87=(a87+a60);
  a59=(a59/a91);
  a87=(a87*a59);
  a29=(a29+a87);
  a27=(a27*a29);
  a87=8.0000000000000004e-01;
  a91=(a5*a8);
  a60=(a6*a9);
  a91=(a91-a60);
  a91=(a11*a91);
  a91=(a91*a53);
  a60=casadi_sq(a6);
  a44=casadi_sq(a5);
  a60=(a60-a44);
  a44=casadi_sq(a8);
  a60=(a60+a44);
  a44=casadi_sq(a9);
  a60=(a60-a44);
  a60=(a60*a40);
  a91=(a91+a60);
  a60=(a8*a9);
  a44=(a6*a5);
  a60=(a60+a44);
  a60=(a11*a60);
  a60=(a60*a59);
  a91=(a91+a60);
  a91=(a87*a91);
  a91=(a27*a91);
  a19=(a19*a91);
  a60=(a23*a18);
  a44=(a25*a21);
  a60=(a60-a44);
  a60=(a60/a24);
  a44=casadi_sq(a6);
  a88=casadi_sq(a5);
  a44=(a44+a88);
  a88=casadi_sq(a8);
  a44=(a44-a88);
  a88=casadi_sq(a9);
  a44=(a44-a88);
  a44=(a44*a53);
  a53=(a5*a8);
  a88=(a6*a9);
  a53=(a53+a88);
  a53=(a11*a53);
  a53=(a53*a40);
  a44=(a44+a53);
  a5=(a5*a9);
  a6=(a6*a8);
  a5=(a5-a6);
  a5=(a11*a5);
  a5=(a5*a59);
  a44=(a44+a5);
  a44=(a87*a44);
  a44=(a27*a44);
  a60=(a60*a44);
  a19=(a19-a60);
  a60=(a14*a25);
  a5=(a20*a23);
  a60=(a60-a5);
  a60=(a60/a24);
  a5=4.0000000000000002e-01;
  a5=(a5*a29);
  a2=(a2-a5);
  a87=(a87*a29);
  a2=(a2-a87);
  a27=(a27*a2);
  a60=(a60*a27);
  a19=(a19+a60);
  a60=-1.9094623163594840e+09;
  a35=(a35/a36);
  a34=(a34+a35);
  a34=(a34-a38);
  a34=(a34/a39);
  a33=(a33*a34);
  a32=(a32+a33);
  a33=casadi_sq(a34);
  a41=(a41*a33);
  a32=(a32+a41);
  a41=casadi_sq(a34);
  a41=(a41*a34);
  a43=(a43*a41);
  a32=(a32+a43);
  a32=(a31*a32);
  a43=(-a32);
  a43=cos(a43);
  a47=(a47*a34);
  a46=(a46+a47);
  a49=(a49*a33);
  a46=(a46+a49);
  a50=(a50*a41);
  a46=(a46+a50);
  a46=(a31*a46);
  a50=(a46-a32);
  a49=(-a50);
  a49=cos(a49);
  a47=(a43*a49);
  a39=(-a32);
  a39=sin(a39);
  a55=(a55*a34);
  a57=(a57*a33);
  a55=(a55+a57);
  a58=(a58*a41);
  a55=(a55+a58);
  a55=(a31*a55);
  a58=(-a55);
  a58=cos(a58);
  a39=(a39*a58);
  a57=(-a50);
  a57=sin(a57);
  a38=(a39*a57);
  a47=(a47-a38);
  a64=(a64*a34);
  a63=(a63+a64);
  a66=(a66*a33);
  a63=(a63+a66);
  a67=(a67*a41);
  a63=(a63+a67);
  a67=casadi_sq(a63);
  a67=(a62-a67);
  a67=(a61*a67);
  a69=(a69*a34);
  a68=(a68+a69);
  a70=(a70*a33);
  a68=(a68+a70);
  a31=(a31*a68);
  a31=(a31-a46);
  a46=(a11*a63);
  a68=casadi_sq(a63);
  a68=(a63*a68);
  a68=(a68/a30);
  a46=(a46-a68);
  a68=casadi_sq(a63);
  a68=casadi_sq(a68);
  a68=(a63*a68);
  a42=(a42*a68);
  a46=(a46+a42);
  a42=sin(a31);
  a46=(a46*a42);
  a42=casadi_sq(a63);
  a71=(a71*a42);
  a42=casadi_sq(a63);
  a42=casadi_sq(a42);
  a73=(a73*a42);
  a71=(a71-a73);
  a73=(a11*a31);
  a73=sin(a73);
  a71=(a71*a73);
  a46=(a46+a71);
  a71=casadi_sq(a63);
  a71=(a63*a71);
  a72=(a72*a71);
  a71=casadi_sq(a63);
  a71=casadi_sq(a71);
  a71=(a63*a71);
  a75=(a75*a71);
  a72=(a72-a75);
  a75=(a76*a31);
  a75=sin(a75);
  a72=(a72*a75);
  a46=(a46+a72);
  a72=casadi_sq(a63);
  a72=casadi_sq(a72);
  a74=(a74*a72);
  a30=(a30*a31);
  a30=sin(a30);
  a74=(a74*a30);
  a46=(a46+a74);
  a74=casadi_sq(a63);
  a74=casadi_sq(a74);
  a74=(a63*a74);
  a77=(a77*a74);
  a79=(a79*a31);
  a79=sin(a79);
  a77=(a77*a79);
  a46=(a46+a77);
  a31=(a31+a46);
  a46=cos(a31);
  a46=(a63*a46);
  a46=(a62+a46);
  a67=(a67/a46);
  a46=cos(a31);
  a46=(a67*a46);
  a77=(a47*a46);
  a79=(-a50);
  a79=sin(a79);
  a43=(a43*a79);
  a50=(-a50);
  a50=cos(a50);
  a39=(a39*a50);
  a43=(a43+a39);
  a39=sin(a31);
  a67=(a67*a39);
  a39=(a43*a67);
  a77=(a77+a39);
  a77=(a77+a13);
  a39=(a60*a77);
  a74=casadi_sq(a77);
  a30=(-a32);
  a30=cos(a30);
  a30=(a30*a58);
  a58=(a30*a50);
  a32=(-a32);
  a32=sin(a32);
  a79=(a32*a79);
  a58=(a58-a79);
  a79=(a58*a67);
  a32=(a32*a49);
  a30=(a30*a57);
  a32=(a32+a30);
  a30=(a32*a46);
  a79=(a79-a30);
  a30=(a94*a79);
  a55=(-a55);
  a55=sin(a55);
  a57=(a55*a57);
  a46=(a57*a46);
  a55=(a55*a50);
  a67=(a55*a67);
  a46=(a46-a67);
  a67=(a95*a46);
  a30=(a30+a67);
  a30=(a30+a3);
  a67=casadi_sq(a30);
  a74=(a74+a67);
  a79=(a56*a79);
  a46=(a94*a46);
  a79=(a79+a46);
  a79=(a79+a16);
  a46=casadi_sq(a79);
  a74=(a74+a46);
  a46=sqrt(a74);
  a46=(a46*a74);
  a39=(a39/a46);
  a19=(a19+a39);
  a74=(a7*a19);
  a67=(a26*a18);
  a50=(a25*a28);
  a67=(a67-a50);
  a67=(a67/a24);
  a67=(a67*a44);
  a18=(a22*a18);
  a50=(a20*a28);
  a18=(a18-a50);
  a18=(a18/a24);
  a18=(a18*a91);
  a67=(a67-a18);
  a25=(a22*a25);
  a20=(a20*a26);
  a25=(a25-a20);
  a25=(a25/a24);
  a25=(a25*a27);
  a67=(a67-a25);
  a25=(a60*a30);
  a25=(a25/a46);
  a67=(a67+a25);
  a20=(a10*a67);
  a74=(a74+a20);
  a20=(a22*a21);
  a18=(a14*a28);
  a20=(a20-a18);
  a20=(a20/a24);
  a20=(a20*a91);
  a21=(a26*a21);
  a28=(a23*a28);
  a21=(a21-a28);
  a21=(a21/a24);
  a21=(a21*a44);
  a20=(a20-a21);
  a22=(a22*a23);
  a14=(a14*a26);
  a22=(a22-a14);
  a22=(a22/a24);
  a22=(a22*a27);
  a20=(a20+a22);
  a60=(a60*a79);
  a60=(a60/a46);
  a20=(a20+a60);
  a46=(a84*a20);
  a74=(a74+a46);
  a46=-3.9860044150000002e+05;
  a22=casadi_sq(a89);
  a27=casadi_sq(a92);
  a22=(a22+a27);
  a27=casadi_sq(a93);
  a22=(a22+a27);
  a27=sqrt(a22);
  a27=(a27*a22);
  a46=(a46/a27);
  a89=(a46*a89);
  a27=3.9860044150000002e+05;
  a27=(a27/a85);
  a89=(a89+a27);
  a74=(a74+a89);
  a89=-2.;
  a27=(a3*a78);
  a22=(a16*a45);
  a27=(a27-a22);
  a22=(a27/a85);
  a24=casadi_sq(a27);
  a14=(a16*a54);
  a26=(a13*a78);
  a14=(a14-a26);
  a26=casadi_sq(a14);
  a24=(a24+a26);
  a26=(a13*a45);
  a23=(a3*a54);
  a26=(a26-a23);
  a23=casadi_sq(a26);
  a24=(a24+a23);
  a23=sqrt(a24);
  a21=(a27/a23);
  a44=(a39*a21);
  a28=(a14/a23);
  a91=(a25*a28);
  a44=(a44+a91);
  a91=(a26/a23);
  a18=(a60*a91);
  a44=(a44+a18);
  a18=(a44*a13);
  a18=(a18/a23);
  a22=(a22-a18);
  a18=(a81*a22);
  a50=(a14/a85);
  a49=(a44*a3);
  a49=(a49/a23);
  a50=(a50-a49);
  a49=(a83*a50);
  a18=(a18+a49);
  a49=(a26/a85);
  a44=(a44*a16);
  a44=(a44/a23);
  a49=(a49-a44);
  a44=(a86*a49);
  a18=(a18+a44);
  a44=(a89*a18);
  a44=(a44*a15);
  a72=(a80*a22);
  a75=(a82*a50);
  a72=(a72+a75);
  a75=(a37*a49);
  a72=(a72+a75);
  a75=(a89*a72);
  a75=(a75*a0);
  a44=(a44-a75);
  a22=(a7*a22);
  a50=(a10*a50);
  a22=(a22+a50);
  a49=(a84*a49);
  a22=(a22+a49);
  a49=(a22*a92);
  a50=(a18*a90);
  a49=(a49-a50);
  a50=(a18*a49);
  a75=(a72*a90);
  a71=(a22*a93);
  a75=(a75-a71);
  a71=(a72*a75);
  a50=(a50-a71);
  a44=(a44-a50);
  a50=(a16*a25);
  a71=(a3*a60);
  a50=(a50-a71);
  a71=(a50/a85);
  a73=(a13*a54);
  a42=(a3*a45);
  a73=(a73+a42);
  a42=(a16*a78);
  a73=(a73+a42);
  a11=(a11*a73);
  a73=(a11*a27);
  a42=casadi_sq(a85);
  a73=(a73/a42);
  a71=(a71-a73);
  a73=(a27*a50);
  a68=(a13*a60);
  a70=(a16*a39);
  a68=(a68-a70);
  a70=(a14*a68);
  a73=(a73+a70);
  a70=(a3*a39);
  a33=(a13*a25);
  a70=(a70-a33);
  a33=(a26*a70);
  a73=(a73+a33);
  a33=(a39*a21);
  a69=(a25*a28);
  a33=(a33+a69);
  a69=(a60*a91);
  a33=(a33+a69);
  a73=(a73*a33);
  a33=(a73*a13);
  a69=(a23*a24);
  a33=(a33/a69);
  a71=(a71+a33);
  a33=1.9094623163594840e+09;
  a34=1.3271244001798700e+11;
  a41=casadi_sq(a63);
  a62=(a62-a41);
  a61=(a61*a62);
  a34=(a34/a61);
  a34=sqrt(a34);
  a61=cos(a31);
  a63=(a63+a61);
  a63=(a34*a63);
  a43=(a43*a63);
  a31=sin(a31);
  a34=(a34*a31);
  a47=(a47*a34);
  a43=(a43-a47);
  a43=(a43+a54);
  a47=(a77*a43);
  a32=(a32*a34);
  a58=(a58*a63);
  a32=(a32+a58);
  a58=(a94*a32);
  a57=(a57*a34);
  a55=(a55*a63);
  a57=(a57+a55);
  a95=(a95*a57);
  a58=(a58-a95);
  a58=(a58+a45);
  a95=(a30*a58);
  a47=(a47+a95);
  a56=(a56*a32);
  a94=(a94*a57);
  a56=(a56-a94);
  a56=(a56+a78);
  a94=(a79*a56);
  a47=(a47+a94);
  a76=(a76*a47);
  a47=(a76*a77);
  a77=casadi_sq(a77);
  a94=casadi_sq(a30);
  a77=(a77+a94);
  a94=casadi_sq(a79);
  a77=(a77+a94);
  a94=sqrt(a77);
  a57=casadi_sq(a77);
  a57=(a94*a57);
  a47=(a47/a57);
  a94=(a94*a77);
  a43=(a43/a94);
  a47=(a47-a43);
  a47=(a33*a47);
  a47=(a47*a21);
  a30=(a76*a30);
  a30=(a30/a57);
  a58=(a58/a94);
  a30=(a30-a58);
  a30=(a33*a30);
  a30=(a30*a28);
  a47=(a47+a30);
  a76=(a76*a79);
  a76=(a76/a57);
  a56=(a56/a94);
  a76=(a76-a56);
  a33=(a33*a76);
  a33=(a33*a91);
  a47=(a47+a33);
  a33=(a50/a23);
  a50=(a27*a50);
  a76=(a14*a68);
  a50=(a50+a76);
  a76=(a26*a70);
  a50=(a50+a76);
  a27=(a50*a27);
  a24=(a23*a24);
  a27=(a27/a24);
  a33=(a33-a27);
  a33=(a39*a33);
  a27=(a68/a23);
  a76=(a50*a14);
  a76=(a76/a24);
  a27=(a27-a76);
  a27=(a25*a27);
  a33=(a33+a27);
  a27=(a70/a23);
  a50=(a50*a26);
  a50=(a50/a24);
  a27=(a27-a50);
  a27=(a60*a27);
  a33=(a33+a27);
  a47=(a47-a33);
  a13=(a47*a13);
  a13=(a13/a23);
  a71=(a71+a13);
  a39=(a39*a21);
  a25=(a25*a28);
  a39=(a39+a25);
  a60=(a60*a91);
  a39=(a39+a60);
  a54=(a39*a54);
  a54=(a54/a23);
  a71=(a71-a54);
  a54=(a81*a71);
  a68=(a68/a85);
  a14=(a11*a14);
  a14=(a14/a42);
  a68=(a68-a14);
  a14=(a73*a3);
  a14=(a14/a69);
  a68=(a68+a14);
  a3=(a47*a3);
  a3=(a3/a23);
  a68=(a68+a3);
  a45=(a39*a45);
  a45=(a45/a23);
  a68=(a68-a45);
  a45=(a83*a68);
  a54=(a54+a45);
  a70=(a70/a85);
  a11=(a11*a26);
  a11=(a11/a42);
  a70=(a70-a11);
  a73=(a73*a16);
  a73=(a73/a69);
  a70=(a70+a73);
  a47=(a47*a16);
  a47=(a47/a23);
  a70=(a70+a47);
  a39=(a39*a78);
  a39=(a39/a23);
  a70=(a70-a39);
  a39=(a86*a70);
  a54=(a54+a39);
  a39=(a54*a93);
  a23=(a80*a71);
  a78=(a82*a68);
  a23=(a23+a78);
  a78=(a37*a70);
  a23=(a23+a78);
  a78=(a23*a92);
  a39=(a39-a78);
  a44=(a44-a39);
  a74=(a74+a44);
  a74=(a1*a74);
  a44=(a4-a74);
  a39=arg[4]? arg[4][0] : 0;
  a44=(a44*a39);
  a4=(a4-a74);
  a44=(a44*a4);
  a12=(a12+a44);
  a44=arg[1]? arg[1][11] : 0;
  a81=(a81*a19);
  a83=(a83*a67);
  a81=(a81+a83);
  a86=(a86*a20);
  a81=(a81+a86);
  a86=(a46*a92);
  a81=(a81+a86);
  a86=(a89*a72);
  a86=(a86*a17);
  a83=(a89*a22);
  a83=(a83*a15);
  a86=(a86-a83);
  a83=(a18*a93);
  a15=(a72*a92);
  a83=(a83-a15);
  a72=(a72*a83);
  a49=(a22*a49);
  a72=(a72-a49);
  a86=(a86-a72);
  a23=(a23*a90);
  a7=(a7*a71);
  a10=(a10*a68);
  a7=(a7+a10);
  a84=(a84*a70);
  a7=(a7+a84);
  a84=(a7*a93);
  a23=(a23-a84);
  a86=(a86-a23);
  a81=(a81+a86);
  a81=(a1*a81);
  a86=(a44-a81);
  a86=(a86*a39);
  a44=(a44-a81);
  a86=(a86*a44);
  a12=(a12+a86);
  a86=arg[1]? arg[1][12] : 0;
  a80=(a80*a19);
  a82=(a82*a67);
  a80=(a80+a82);
  a37=(a37*a20);
  a80=(a80+a37);
  a46=(a46*a93);
  a80=(a80+a46);
  a46=(a89*a22);
  a46=(a46*a0);
  a89=(a89*a18);
  a89=(a89*a17);
  a46=(a46-a89);
  a22=(a22*a75);
  a18=(a18*a83);
  a22=(a22-a18);
  a46=(a46-a22);
  a7=(a7*a92);
  a54=(a54*a90);
  a7=(a7-a54);
  a46=(a46-a7);
  a80=(a80+a46);
  a1=(a1*a80);
  a80=(a86-a1);
  a80=(a80*a39);
  a86=(a86-a1);
  a80=(a80*a86);
  a12=(a12+a80);
  if (res[0]!=0) res[0][0]=a12;
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

CASADI_SYMBOL_EXPORT casadi_int L_n_in(void) { return 7;}

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
    case 2: return "k1";
    case 3: return "k2";
    case 4: return "k3";
    case 5: return "t";
    case 6: return "XChief";
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
    case 4: return casadi_s1;
    case 5: return casadi_s1;
    case 6: return casadi_s2;
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
  if (sz_arg) *sz_arg = 7;
  if (sz_res) *sz_res = 1;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 0;
  return 0;
}

#ifdef MATLAB_MEX_FILE
void mex_L(int resc, mxArray *resv[], int argc, const mxArray *argv[]) {
  casadi_int i;
  casadi_real w[133];
  casadi_int *iw = 0;
  const casadi_real* arg[7] = {0};
  casadi_real* res[1] = {0};
  if (argc>7) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"L\" failed. Too many input arguments (%d, max 7)", argc);
  if (resc>1) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"L\" failed. Too many output arguments (%d, max 1)", resc);
  if (--argc>=0) arg[0] = casadi_from_mex(argv[0], w, casadi_s0, w+37);
  if (--argc>=0) arg[1] = casadi_from_mex(argv[1], w+13, casadi_s0, w+37);
  if (--argc>=0) arg[2] = casadi_from_mex(argv[2], w+26, casadi_s1, w+37);
  if (--argc>=0) arg[3] = casadi_from_mex(argv[3], w+27, casadi_s1, w+37);
  if (--argc>=0) arg[4] = casadi_from_mex(argv[4], w+28, casadi_s1, w+37);
  if (--argc>=0) arg[5] = casadi_from_mex(argv[5], w+29, casadi_s1, w+37);
  if (--argc>=0) arg[6] = casadi_from_mex(argv[6], w+30, casadi_s2, w+37);
  --resc;
  res[0] = w+36;
  i = L(arg, res, iw, w+37, 0);
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
  casadi_real w[133];
  const casadi_real* arg[7];
  casadi_real* res[1];
  arg[0] = w+0;
  arg[1] = w+13;
  arg[2] = w+26;
  arg[3] = w+27;
  arg[4] = w+28;
  arg[5] = w+29;
  arg[6] = w+30;
  res[0] = w+36;
  a = w;
  for (j=0; j<36; ++j) if (scanf("%lg", a++)<=0) return 2;
  flag = L(arg, res, iw, w+37, 0);
  if (flag) return flag;
  r = w+36;
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
