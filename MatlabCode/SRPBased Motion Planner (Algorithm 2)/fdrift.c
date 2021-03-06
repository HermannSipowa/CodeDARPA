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
  #define CASADI_PREFIX(ID) fdrift_ ## ID
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

static const casadi_int casadi_s0[5] = {1, 1, 0, 1, 0};
static const casadi_int casadi_s1[17] = {13, 1, 0, 13, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
static const casadi_int casadi_s2[10] = {6, 1, 0, 6, 0, 1, 2, 3, 4, 5};

/* fdrift:(i0,i1[13],i2[6])->(o0[13]) */
static int casadi_f0(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem) {
  casadi_real a0, a1, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a2, a20, a21, a22, a23, a24, a25, a26, a27, a28, a29, a3, a30, a31, a32, a33, a34, a35, a36, a37, a38, a39, a4, a40, a41, a42, a43, a44, a45, a46, a47, a48, a49, a5, a50, a51, a52, a53, a54, a55, a56, a57, a58, a59, a6, a60, a61, a62, a63, a64, a65, a66, a67, a68, a69, a7, a70, a71, a72, a73, a74, a75, a76, a77, a78, a79, a8, a80, a81, a82, a83, a84, a85, a86, a87, a88, a89, a9, a90, a91, a92, a93;
  a0=1.1259418598742898e+05;
  a1=5.0000000000000000e-01;
  a2=arg[1]? arg[1][4] : 0;
  a3=(a1*a2);
  a4=arg[1]? arg[1][1] : 0;
  a5=arg[1]? arg[1][0] : 0;
  a6=casadi_sq(a5);
  a7=casadi_sq(a4);
  a6=(a6+a7);
  a7=arg[1]? arg[1][2] : 0;
  a8=casadi_sq(a7);
  a6=(a6+a8);
  a8=arg[1]? arg[1][3] : 0;
  a9=casadi_sq(a8);
  a6=(a6+a9);
  a6=sqrt(a6);
  a4=(a4/a6);
  a3=(a3*a4);
  a9=arg[1]? arg[1][5] : 0;
  a10=(a1*a9);
  a7=(a7/a6);
  a10=(a10*a7);
  a3=(a3+a10);
  a10=arg[1]? arg[1][6] : 0;
  a11=(a1*a10);
  a8=(a8/a6);
  a11=(a11*a8);
  a3=(a3+a11);
  a3=(a0*a3);
  a3=(-a3);
  if (res[0]!=0) res[0][0]=a3;
  a3=(a1*a2);
  a5=(a5/a6);
  a3=(a3*a5);
  a6=(a1*a10);
  a6=(a6*a7);
  a3=(a3+a6);
  a6=(a1*a9);
  a6=(a6*a8);
  a3=(a3-a6);
  a3=(a0*a3);
  if (res[0]!=0) res[0][1]=a3;
  a3=(a1*a9);
  a3=(a3*a5);
  a6=(a1*a10);
  a6=(a6*a4);
  a3=(a3-a6);
  a6=(a1*a2);
  a6=(a6*a8);
  a3=(a3+a6);
  a3=(a0*a3);
  if (res[0]!=0) res[0][2]=a3;
  a3=(a1*a10);
  a3=(a3*a5);
  a6=(a1*a9);
  a6=(a6*a4);
  a3=(a3+a6);
  a6=(a1*a2);
  a6=(a6*a7);
  a3=(a3-a6);
  a3=(a0*a3);
  if (res[0]!=0) res[0][3]=a3;
  a3=1.2121212121212120e-01;
  a6=(a3*a9);
  a11=(a10*a6);
  a12=2.4242424242424240e-01;
  a13=(a12*a10);
  a14=(a9*a13);
  a11=(a11-a14);
  a11=(a11/a3);
  a11=(a0*a11);
  if (res[0]!=0) res[0][4]=a11;
  a13=(a2*a13);
  a11=(a3*a2);
  a10=(a10*a11);
  a13=(a13-a10);
  a13=(a13/a3);
  a13=(a0*a13);
  if (res[0]!=0) res[0][5]=a13;
  a9=(a9*a11);
  a2=(a2*a6);
  a9=(a9-a2);
  a9=(a9/a12);
  a9=(a0*a9);
  if (res[0]!=0) res[0][6]=a9;
  a9=arg[1]? arg[1][10] : 0;
  a12=(a0*a9);
  if (res[0]!=0) res[0][7]=a12;
  a12=arg[1]? arg[1][11] : 0;
  a2=(a0*a12);
  if (res[0]!=0) res[0][8]=a2;
  a2=arg[1]? arg[1][12] : 0;
  a6=(a0*a2);
  if (res[0]!=0) res[0][9]=a6;
  a6=arg[2]? arg[2][0] : 0;
  a11=casadi_sq(a6);
  a13=arg[2]? arg[2][1] : 0;
  a3=casadi_sq(a13);
  a11=(a11+a3);
  a3=arg[2]? arg[2][2] : 0;
  a10=casadi_sq(a3);
  a11=(a11+a10);
  a11=sqrt(a11);
  a10=(a6/a11);
  a14=2.;
  a15=(a4*a7);
  a16=(a5*a8);
  a15=(a15+a16);
  a15=(a14*a15);
  a16=casadi_sq(a5);
  a17=casadi_sq(a4);
  a16=(a16-a17);
  a17=casadi_sq(a7);
  a16=(a16-a17);
  a17=casadi_sq(a8);
  a16=(a16+a17);
  a17=(a15*a16);
  a18=(a4*a8);
  a19=(a5*a7);
  a18=(a18-a19);
  a18=(a14*a18);
  a19=(a7*a8);
  a20=(a5*a4);
  a19=(a19-a20);
  a19=(a14*a19);
  a20=(a18*a19);
  a17=(a17-a20);
  a20=casadi_sq(a5);
  a21=casadi_sq(a4);
  a20=(a20+a21);
  a21=casadi_sq(a7);
  a20=(a20-a21);
  a21=casadi_sq(a8);
  a20=(a20-a21);
  a21=casadi_sq(a5);
  a22=casadi_sq(a4);
  a21=(a21-a22);
  a22=casadi_sq(a7);
  a21=(a21+a22);
  a22=casadi_sq(a8);
  a21=(a21-a22);
  a22=(a21*a16);
  a23=(a7*a8);
  a24=(a5*a4);
  a23=(a23+a24);
  a23=(a14*a23);
  a24=(a23*a19);
  a22=(a22-a24);
  a22=(a20*a22);
  a24=(a4*a7);
  a25=(a5*a8);
  a24=(a24-a25);
  a24=(a14*a24);
  a25=(a24*a16);
  a26=(a4*a8);
  a27=(a5*a7);
  a26=(a26+a27);
  a26=(a14*a26);
  a27=(a23*a26);
  a25=(a25-a27);
  a25=(a15*a25);
  a22=(a22-a25);
  a25=(a24*a19);
  a27=(a21*a26);
  a25=(a25-a27);
  a25=(a18*a25);
  a22=(a22+a25);
  a17=(a17/a22);
  a25=1.0000000000000000e-03;
  a27=5.4317577262067082e-03;
  a28=4.;
  a29=1.7453292519943295e-02;
  a30=1.7487317400000001e+02;
  a31=-2.4109079999999999e-01;
  a32=2.4562962500000000e+06;
  a33=arg[0]? arg[0][0] : 0;
  a33=(a0*a33);
  a34=86400.;
  a35=(a33/a34);
  a35=(a32+a35);
  a36=2451545.;
  a35=(a35-a36);
  a37=36525.;
  a35=(a35/a37);
  a38=(a31*a35);
  a38=(a30+a38);
  a39=4.0670000000000002e-05;
  a40=casadi_sq(a35);
  a41=(a39*a40);
  a38=(a38+a41);
  a41=-1.3270000000000000e-06;
  a42=casadi_sq(a35);
  a42=(a42*a35);
  a43=(a41*a42);
  a38=(a38+a43);
  a38=(a29*a38);
  a43=(-a38);
  a43=cos(a43);
  a44=1.0293734800000000e+02;
  a45=3.2255570000000000e-01;
  a46=(a45*a35);
  a46=(a44+a46);
  a47=1.5025999999999999e-04;
  a48=(a47*a40);
  a46=(a46+a48);
  a48=4.7800000000000002e-07;
  a49=(a48*a42);
  a46=(a46+a49);
  a46=(a29*a46);
  a49=(a46-a38);
  a50=(-a49);
  a50=cos(a50);
  a51=(a43*a50);
  a52=(-a38);
  a52=sin(a52);
  a53=1.3054600000000000e-02;
  a54=(a53*a35);
  a55=-9.3100000000000006e-06;
  a56=(a55*a40);
  a54=(a54+a56);
  a56=-3.4000000000000000e-08;
  a57=(a56*a42);
  a54=(a54+a57);
  a54=(a29*a54);
  a57=(-a54);
  a57=cos(a57);
  a52=(a52*a57);
  a58=(-a49);
  a58=sin(a58);
  a59=(a52*a58);
  a51=(a51-a59);
  a59=1.4959802299063239e+08;
  a60=1.;
  a61=1.6708620000000000e-02;
  a62=-4.2036999999999997e-05;
  a63=(a62*a35);
  a63=(a61+a63);
  a64=-1.2360000000000000e-07;
  a65=(a64*a40);
  a63=(a63+a65);
  a65=3.9999999999999998e-11;
  a42=(a65*a42);
  a63=(a63+a42);
  a42=casadi_sq(a63);
  a42=(a60-a42);
  a42=(a59*a42);
  a66=1.0046644900000000e+02;
  a67=3.5999372851900000e+04;
  a35=(a67*a35);
  a35=(a66+a35);
  a68=-5.6799999999999998e-06;
  a40=(a68*a40);
  a35=(a35+a40);
  a35=(a29*a35);
  a35=(a35-a46);
  a46=(a14*a63);
  a40=casadi_sq(a63);
  a40=(a63*a40);
  a40=(a40/a28);
  a46=(a46-a40);
  a40=5.2083333333333336e-02;
  a69=casadi_sq(a63);
  a69=casadi_sq(a69);
  a69=(a63*a69);
  a69=(a40*a69);
  a46=(a46+a69);
  a69=sin(a35);
  a46=(a46*a69);
  a69=1.2500000000000000e+00;
  a70=casadi_sq(a63);
  a70=(a69*a70);
  a71=4.5833333333333331e-01;
  a72=casadi_sq(a63);
  a72=casadi_sq(a72);
  a72=(a71*a72);
  a70=(a70-a72);
  a72=(a14*a35);
  a72=sin(a72);
  a70=(a70*a72);
  a46=(a46+a70);
  a70=1.0833333333333333e+00;
  a72=casadi_sq(a63);
  a72=(a63*a72);
  a72=(a70*a72);
  a73=6.7187500000000000e-01;
  a74=casadi_sq(a63);
  a74=casadi_sq(a74);
  a74=(a63*a74);
  a74=(a73*a74);
  a72=(a72-a74);
  a74=3.;
  a75=(a74*a35);
  a75=sin(a75);
  a72=(a72*a75);
  a46=(a46+a72);
  a72=1.0729166666666667e+00;
  a75=casadi_sq(a63);
  a75=casadi_sq(a75);
  a75=(a72*a75);
  a76=(a28*a35);
  a76=sin(a76);
  a75=(a75*a76);
  a46=(a46+a75);
  a75=1.1427083333333334e+00;
  a76=casadi_sq(a63);
  a76=casadi_sq(a76);
  a76=(a63*a76);
  a76=(a75*a76);
  a77=5.;
  a78=(a77*a35);
  a78=sin(a78);
  a76=(a76*a78);
  a46=(a46+a76);
  a35=(a35+a46);
  a46=cos(a35);
  a63=(a63*a46);
  a63=(a60+a63);
  a42=(a42/a63);
  a63=cos(a35);
  a63=(a42*a63);
  a51=(a51*a63);
  a46=(-a49);
  a46=sin(a46);
  a43=(a43*a46);
  a49=(-a49);
  a49=cos(a49);
  a52=(a52*a49);
  a43=(a43+a52);
  a35=sin(a35);
  a42=(a42*a35);
  a43=(a43*a42);
  a51=(a51+a43);
  a43=arg[2]? arg[2][4] : 0;
  a35=(a6*a43);
  a52=arg[2]? arg[2][3] : 0;
  a76=(a13*a52);
  a35=(a35-a76);
  a76=arg[2]? arg[2][5] : 0;
  a78=(a13*a76);
  a79=(a3*a43);
  a78=(a78-a79);
  a79=casadi_sq(a78);
  a80=(a3*a52);
  a81=(a6*a76);
  a80=(a80-a81);
  a81=casadi_sq(a80);
  a79=(a79+a81);
  a81=casadi_sq(a35);
  a79=(a79+a81);
  a79=sqrt(a79);
  a35=(a35/a79);
  a81=(a35*a10);
  a78=(a78/a79);
  a82=(a3/a11);
  a83=(a78*a82);
  a81=(a81-a83);
  a80=(a80/a79);
  a79=(a80*a82);
  a11=(a13/a11);
  a83=(a35*a11);
  a79=(a79-a83);
  a83=casadi_sq(a79);
  a84=casadi_sq(a81);
  a83=(a83+a84);
  a84=(a78*a11);
  a85=(a80*a10);
  a84=(a84-a85);
  a85=casadi_sq(a84);
  a83=(a83+a85);
  a83=sqrt(a83);
  a81=(a81/a83);
  a85=(a81*a35);
  a84=(a84/a83);
  a86=(a84*a80);
  a85=(a85-a86);
  a86=(a81*a35);
  a87=(a84*a80);
  a86=(a86-a87);
  a86=(a10*a86);
  a79=(a79/a83);
  a83=(a79*a35);
  a87=(a84*a78);
  a83=(a83-a87);
  a83=(a11*a83);
  a86=(a86-a83);
  a83=(a79*a80);
  a87=(a81*a78);
  a83=(a83-a87);
  a83=(a82*a83);
  a86=(a86+a83);
  a85=(a85/a86);
  a83=casadi_sq(a6);
  a87=casadi_sq(a13);
  a83=(a83+a87);
  a87=casadi_sq(a3);
  a83=(a83+a87);
  a87=sqrt(a83);
  a88=arg[1]? arg[1][7] : 0;
  a87=(a87+a88);
  a85=(a85*a87);
  a89=(a11*a35);
  a90=(a82*a80);
  a89=(a89-a90);
  a89=(a89/a86);
  a90=arg[1]? arg[1][8] : 0;
  a89=(a89*a90);
  a85=(a85-a89);
  a89=(a11*a84);
  a91=(a82*a81);
  a89=(a89-a91);
  a89=(a89/a86);
  a91=arg[1]? arg[1][9] : 0;
  a89=(a89*a91);
  a85=(a85+a89);
  a51=(a51+a85);
  a85=1.4959787069999999e+08;
  a89=(a51/a85);
  a89=casadi_sq(a89);
  a92=9.1748200035787253e-01;
  a93=(-a38);
  a93=cos(a93);
  a93=(a93*a57);
  a57=(a93*a49);
  a38=(-a38);
  a38=sin(a38);
  a46=(a38*a46);
  a57=(a57-a46);
  a57=(a57*a42);
  a38=(a38*a50);
  a93=(a93*a58);
  a38=(a38+a93);
  a38=(a38*a63);
  a57=(a57-a38);
  a38=(a92*a57);
  a93=-3.9777729827042280e-01;
  a54=(-a54);
  a54=sin(a54);
  a58=(a54*a58);
  a58=(a58*a63);
  a54=(a54*a49);
  a54=(a54*a42);
  a58=(a58-a54);
  a54=(a93*a58);
  a38=(a38+a54);
  a54=(a10*a35);
  a42=(a82*a78);
  a54=(a54-a42);
  a54=(a54/a86);
  a54=(a54*a90);
  a42=(a79*a35);
  a49=(a84*a78);
  a42=(a42-a49);
  a42=(a42/a86);
  a42=(a42*a87);
  a54=(a54-a42);
  a42=(a10*a84);
  a49=(a82*a79);
  a42=(a42-a49);
  a42=(a42/a86);
  a42=(a42*a91);
  a54=(a54-a42);
  a38=(a38+a54);
  a54=(a38/a85);
  a54=casadi_sq(a54);
  a89=(a89+a54);
  a54=3.9777729827042280e-01;
  a57=(a54*a57);
  a58=(a92*a58);
  a57=(a57+a58);
  a58=(a79*a80);
  a42=(a81*a78);
  a58=(a58-a42);
  a58=(a58/a86);
  a58=(a58*a87);
  a42=(a10*a80);
  a49=(a11*a78);
  a42=(a42-a49);
  a42=(a42/a86);
  a42=(a42*a90);
  a58=(a58-a42);
  a42=(a10*a81);
  a49=(a11*a79);
  a42=(a42-a49);
  a42=(a42/a86);
  a42=(a42*a91);
  a58=(a58+a42);
  a57=(a57+a58);
  a85=(a57/a85);
  a85=casadi_sq(a85);
  a89=(a89+a85);
  a89=(a28*a89);
  a27=(a27/a89);
  a25=(a25*a27);
  a27=(a4*a8);
  a89=(a5*a7);
  a27=(a27+a89);
  a27=(a14*a27);
  a89=casadi_sq(a51);
  a85=casadi_sq(a38);
  a89=(a89+a85);
  a85=casadi_sq(a57);
  a89=(a89+a85);
  a89=sqrt(a89);
  a51=(a51/a89);
  a27=(a27*a51);
  a85=(a7*a8);
  a58=(a5*a4);
  a85=(a85-a58);
  a85=(a14*a85);
  a38=(a38/a89);
  a85=(a85*a38);
  a27=(a27+a85);
  a85=casadi_sq(a5);
  a58=casadi_sq(a4);
  a85=(a85-a58);
  a58=casadi_sq(a7);
  a85=(a85-a58);
  a58=casadi_sq(a8);
  a85=(a85+a58);
  a57=(a57/a89);
  a85=(a85*a57);
  a27=(a27+a85);
  a25=(a25*a27);
  a85=8.0000000000000004e-01;
  a89=(a4*a7);
  a58=(a5*a8);
  a89=(a89-a58);
  a89=(a14*a89);
  a89=(a89*a51);
  a58=casadi_sq(a5);
  a42=casadi_sq(a4);
  a58=(a58-a42);
  a42=casadi_sq(a7);
  a58=(a58+a42);
  a42=casadi_sq(a8);
  a58=(a58-a42);
  a58=(a58*a38);
  a89=(a89+a58);
  a58=(a7*a8);
  a42=(a5*a4);
  a58=(a58+a42);
  a58=(a14*a58);
  a58=(a58*a57);
  a89=(a89+a58);
  a89=(a85*a89);
  a89=(a25*a89);
  a17=(a17*a89);
  a58=(a21*a16);
  a42=(a23*a19);
  a58=(a58-a42);
  a58=(a58/a22);
  a42=casadi_sq(a5);
  a86=casadi_sq(a4);
  a42=(a42+a86);
  a86=casadi_sq(a7);
  a42=(a42-a86);
  a86=casadi_sq(a8);
  a42=(a42-a86);
  a42=(a42*a51);
  a51=(a4*a7);
  a86=(a5*a8);
  a51=(a51+a86);
  a51=(a14*a51);
  a51=(a51*a38);
  a42=(a42+a51);
  a4=(a4*a8);
  a5=(a5*a7);
  a4=(a4-a5);
  a4=(a14*a4);
  a4=(a4*a57);
  a42=(a42+a4);
  a42=(a85*a42);
  a42=(a25*a42);
  a58=(a58*a42);
  a17=(a17-a58);
  a58=(a15*a23);
  a4=(a18*a21);
  a58=(a58-a4);
  a58=(a58/a22);
  a4=4.0000000000000002e-01;
  a4=(a4*a27);
  a1=(a1-a4);
  a85=(a85*a27);
  a1=(a1-a85);
  a25=(a25*a1);
  a58=(a58*a25);
  a17=(a17+a58);
  a58=-1.9094623163594840e+09;
  a33=(a33/a34);
  a32=(a32+a33);
  a32=(a32-a36);
  a32=(a32/a37);
  a31=(a31*a32);
  a30=(a30+a31);
  a31=casadi_sq(a32);
  a39=(a39*a31);
  a30=(a30+a39);
  a39=casadi_sq(a32);
  a39=(a39*a32);
  a41=(a41*a39);
  a30=(a30+a41);
  a30=(a29*a30);
  a41=(-a30);
  a41=cos(a41);
  a45=(a45*a32);
  a44=(a44+a45);
  a47=(a47*a31);
  a44=(a44+a47);
  a48=(a48*a39);
  a44=(a44+a48);
  a44=(a29*a44);
  a48=(a44-a30);
  a47=(-a48);
  a47=cos(a47);
  a45=(a41*a47);
  a37=(-a30);
  a37=sin(a37);
  a53=(a53*a32);
  a55=(a55*a31);
  a53=(a53+a55);
  a56=(a56*a39);
  a53=(a53+a56);
  a53=(a29*a53);
  a56=(-a53);
  a56=cos(a56);
  a37=(a37*a56);
  a55=(-a48);
  a55=sin(a55);
  a36=(a37*a55);
  a45=(a45-a36);
  a62=(a62*a32);
  a61=(a61+a62);
  a64=(a64*a31);
  a61=(a61+a64);
  a65=(a65*a39);
  a61=(a61+a65);
  a65=casadi_sq(a61);
  a65=(a60-a65);
  a65=(a59*a65);
  a67=(a67*a32);
  a66=(a66+a67);
  a68=(a68*a31);
  a66=(a66+a68);
  a29=(a29*a66);
  a29=(a29-a44);
  a44=(a14*a61);
  a66=casadi_sq(a61);
  a66=(a61*a66);
  a66=(a66/a28);
  a44=(a44-a66);
  a66=casadi_sq(a61);
  a66=casadi_sq(a66);
  a66=(a61*a66);
  a40=(a40*a66);
  a44=(a44+a40);
  a40=sin(a29);
  a44=(a44*a40);
  a40=casadi_sq(a61);
  a69=(a69*a40);
  a40=casadi_sq(a61);
  a40=casadi_sq(a40);
  a71=(a71*a40);
  a69=(a69-a71);
  a71=(a14*a29);
  a71=sin(a71);
  a69=(a69*a71);
  a44=(a44+a69);
  a69=casadi_sq(a61);
  a69=(a61*a69);
  a70=(a70*a69);
  a69=casadi_sq(a61);
  a69=casadi_sq(a69);
  a69=(a61*a69);
  a73=(a73*a69);
  a70=(a70-a73);
  a73=(a74*a29);
  a73=sin(a73);
  a70=(a70*a73);
  a44=(a44+a70);
  a70=casadi_sq(a61);
  a70=casadi_sq(a70);
  a72=(a72*a70);
  a28=(a28*a29);
  a28=sin(a28);
  a72=(a72*a28);
  a44=(a44+a72);
  a72=casadi_sq(a61);
  a72=casadi_sq(a72);
  a72=(a61*a72);
  a75=(a75*a72);
  a77=(a77*a29);
  a77=sin(a77);
  a75=(a75*a77);
  a44=(a44+a75);
  a29=(a29+a44);
  a44=cos(a29);
  a44=(a61*a44);
  a44=(a60+a44);
  a65=(a65/a44);
  a44=cos(a29);
  a44=(a65*a44);
  a75=(a45*a44);
  a77=(-a48);
  a77=sin(a77);
  a41=(a41*a77);
  a48=(-a48);
  a48=cos(a48);
  a37=(a37*a48);
  a41=(a41+a37);
  a37=sin(a29);
  a65=(a65*a37);
  a37=(a41*a65);
  a75=(a75+a37);
  a75=(a75+a6);
  a37=(a58*a75);
  a72=casadi_sq(a75);
  a28=(-a30);
  a28=cos(a28);
  a28=(a28*a56);
  a56=(a28*a48);
  a30=(-a30);
  a30=sin(a30);
  a77=(a30*a77);
  a56=(a56-a77);
  a77=(a56*a65);
  a30=(a30*a47);
  a28=(a28*a55);
  a30=(a30+a28);
  a28=(a30*a44);
  a77=(a77-a28);
  a28=(a92*a77);
  a53=(-a53);
  a53=sin(a53);
  a55=(a53*a55);
  a44=(a55*a44);
  a53=(a53*a48);
  a65=(a53*a65);
  a44=(a44-a65);
  a65=(a93*a44);
  a28=(a28+a65);
  a28=(a28+a13);
  a65=casadi_sq(a28);
  a72=(a72+a65);
  a77=(a54*a77);
  a44=(a92*a44);
  a77=(a77+a44);
  a77=(a77+a3);
  a44=casadi_sq(a77);
  a72=(a72+a44);
  a44=sqrt(a72);
  a44=(a44*a72);
  a37=(a37/a44);
  a17=(a17+a37);
  a72=(a10*a17);
  a65=(a24*a16);
  a48=(a23*a26);
  a65=(a65-a48);
  a65=(a65/a22);
  a65=(a65*a42);
  a16=(a20*a16);
  a48=(a18*a26);
  a16=(a16-a48);
  a16=(a16/a22);
  a16=(a16*a89);
  a65=(a65-a16);
  a23=(a20*a23);
  a18=(a18*a24);
  a23=(a23-a18);
  a23=(a23/a22);
  a23=(a23*a25);
  a65=(a65-a23);
  a23=(a58*a28);
  a23=(a23/a44);
  a65=(a65+a23);
  a18=(a11*a65);
  a72=(a72+a18);
  a18=(a20*a19);
  a16=(a15*a26);
  a18=(a18-a16);
  a18=(a18/a22);
  a18=(a18*a89);
  a19=(a24*a19);
  a26=(a21*a26);
  a19=(a19-a26);
  a19=(a19/a22);
  a19=(a19*a42);
  a18=(a18-a19);
  a20=(a20*a21);
  a15=(a15*a24);
  a20=(a20-a15);
  a20=(a20/a22);
  a20=(a20*a25);
  a18=(a18+a20);
  a58=(a58*a77);
  a58=(a58/a44);
  a18=(a18+a58);
  a44=(a82*a18);
  a72=(a72+a44);
  a44=-3.9860044150000002e+05;
  a20=casadi_sq(a87);
  a25=casadi_sq(a90);
  a20=(a20+a25);
  a25=casadi_sq(a91);
  a20=(a20+a25);
  a25=sqrt(a20);
  a25=(a25*a20);
  a44=(a44/a25);
  a87=(a44*a87);
  a25=3.9860044150000002e+05;
  a25=(a25/a83);
  a87=(a87+a25);
  a72=(a72+a87);
  a87=-2.;
  a25=(a13*a76);
  a20=(a3*a43);
  a25=(a25-a20);
  a20=(a25/a83);
  a22=casadi_sq(a25);
  a15=(a3*a52);
  a24=(a6*a76);
  a15=(a15-a24);
  a24=casadi_sq(a15);
  a22=(a22+a24);
  a24=(a6*a43);
  a21=(a13*a52);
  a24=(a24-a21);
  a21=casadi_sq(a24);
  a22=(a22+a21);
  a21=sqrt(a22);
  a19=(a25/a21);
  a42=(a37*a19);
  a26=(a15/a21);
  a89=(a23*a26);
  a42=(a42+a89);
  a89=(a24/a21);
  a16=(a58*a89);
  a42=(a42+a16);
  a16=(a42*a6);
  a16=(a16/a21);
  a20=(a20-a16);
  a16=(a79*a20);
  a48=(a15/a83);
  a47=(a42*a13);
  a47=(a47/a21);
  a48=(a48-a47);
  a47=(a81*a48);
  a16=(a16+a47);
  a47=(a24/a83);
  a42=(a42*a3);
  a42=(a42/a21);
  a47=(a47-a42);
  a42=(a84*a47);
  a16=(a16+a42);
  a42=(a87*a16);
  a42=(a42*a2);
  a70=(a78*a20);
  a73=(a80*a48);
  a70=(a70+a73);
  a73=(a35*a47);
  a70=(a70+a73);
  a73=(a87*a70);
  a73=(a73*a12);
  a42=(a42-a73);
  a20=(a10*a20);
  a48=(a11*a48);
  a20=(a20+a48);
  a47=(a82*a47);
  a20=(a20+a47);
  a47=(a20*a90);
  a48=(a16*a88);
  a47=(a47-a48);
  a48=(a16*a47);
  a73=(a70*a88);
  a69=(a20*a91);
  a73=(a73-a69);
  a69=(a70*a73);
  a48=(a48-a69);
  a42=(a42-a48);
  a48=(a3*a23);
  a69=(a13*a58);
  a48=(a48-a69);
  a69=(a48/a83);
  a71=(a6*a52);
  a40=(a13*a43);
  a71=(a71+a40);
  a40=(a3*a76);
  a71=(a71+a40);
  a14=(a14*a71);
  a71=(a14*a25);
  a40=casadi_sq(a83);
  a71=(a71/a40);
  a69=(a69-a71);
  a71=(a25*a48);
  a66=(a6*a58);
  a68=(a3*a37);
  a66=(a66-a68);
  a68=(a15*a66);
  a71=(a71+a68);
  a68=(a13*a37);
  a31=(a6*a23);
  a68=(a68-a31);
  a31=(a24*a68);
  a71=(a71+a31);
  a31=(a37*a19);
  a67=(a23*a26);
  a31=(a31+a67);
  a67=(a58*a89);
  a31=(a31+a67);
  a71=(a71*a31);
  a31=(a71*a6);
  a67=(a21*a22);
  a31=(a31/a67);
  a69=(a69+a31);
  a31=1.9094623163594840e+09;
  a32=1.3271244001798700e+11;
  a39=casadi_sq(a61);
  a60=(a60-a39);
  a59=(a59*a60);
  a32=(a32/a59);
  a32=sqrt(a32);
  a59=cos(a29);
  a61=(a61+a59);
  a61=(a32*a61);
  a41=(a41*a61);
  a29=sin(a29);
  a32=(a32*a29);
  a45=(a45*a32);
  a41=(a41-a45);
  a41=(a41+a52);
  a45=(a75*a41);
  a30=(a30*a32);
  a56=(a56*a61);
  a30=(a30+a56);
  a56=(a92*a30);
  a55=(a55*a32);
  a53=(a53*a61);
  a55=(a55+a53);
  a93=(a93*a55);
  a56=(a56-a93);
  a56=(a56+a43);
  a93=(a28*a56);
  a45=(a45+a93);
  a54=(a54*a30);
  a92=(a92*a55);
  a54=(a54-a92);
  a54=(a54+a76);
  a92=(a77*a54);
  a45=(a45+a92);
  a74=(a74*a45);
  a45=(a74*a75);
  a75=casadi_sq(a75);
  a92=casadi_sq(a28);
  a75=(a75+a92);
  a92=casadi_sq(a77);
  a75=(a75+a92);
  a92=sqrt(a75);
  a55=casadi_sq(a75);
  a55=(a92*a55);
  a45=(a45/a55);
  a92=(a92*a75);
  a41=(a41/a92);
  a45=(a45-a41);
  a45=(a31*a45);
  a45=(a45*a19);
  a28=(a74*a28);
  a28=(a28/a55);
  a56=(a56/a92);
  a28=(a28-a56);
  a28=(a31*a28);
  a28=(a28*a26);
  a45=(a45+a28);
  a74=(a74*a77);
  a74=(a74/a55);
  a54=(a54/a92);
  a74=(a74-a54);
  a31=(a31*a74);
  a31=(a31*a89);
  a45=(a45+a31);
  a31=(a48/a21);
  a48=(a25*a48);
  a74=(a15*a66);
  a48=(a48+a74);
  a74=(a24*a68);
  a48=(a48+a74);
  a25=(a48*a25);
  a22=(a21*a22);
  a25=(a25/a22);
  a31=(a31-a25);
  a31=(a37*a31);
  a25=(a66/a21);
  a74=(a48*a15);
  a74=(a74/a22);
  a25=(a25-a74);
  a25=(a23*a25);
  a31=(a31+a25);
  a25=(a68/a21);
  a48=(a48*a24);
  a48=(a48/a22);
  a25=(a25-a48);
  a25=(a58*a25);
  a31=(a31+a25);
  a45=(a45-a31);
  a6=(a45*a6);
  a6=(a6/a21);
  a69=(a69+a6);
  a37=(a37*a19);
  a23=(a23*a26);
  a37=(a37+a23);
  a58=(a58*a89);
  a37=(a37+a58);
  a52=(a37*a52);
  a52=(a52/a21);
  a69=(a69-a52);
  a52=(a79*a69);
  a66=(a66/a83);
  a15=(a14*a15);
  a15=(a15/a40);
  a66=(a66-a15);
  a15=(a71*a13);
  a15=(a15/a67);
  a66=(a66+a15);
  a13=(a45*a13);
  a13=(a13/a21);
  a66=(a66+a13);
  a43=(a37*a43);
  a43=(a43/a21);
  a66=(a66-a43);
  a43=(a81*a66);
  a52=(a52+a43);
  a68=(a68/a83);
  a14=(a14*a24);
  a14=(a14/a40);
  a68=(a68-a14);
  a71=(a71*a3);
  a71=(a71/a67);
  a68=(a68+a71);
  a45=(a45*a3);
  a45=(a45/a21);
  a68=(a68+a45);
  a37=(a37*a76);
  a37=(a37/a21);
  a68=(a68-a37);
  a37=(a84*a68);
  a52=(a52+a37);
  a37=(a52*a91);
  a21=(a78*a69);
  a76=(a80*a66);
  a21=(a21+a76);
  a76=(a35*a68);
  a21=(a21+a76);
  a76=(a21*a90);
  a37=(a37-a76);
  a42=(a42-a37);
  a72=(a72+a42);
  a72=(a0*a72);
  if (res[0]!=0) res[0][10]=a72;
  a79=(a79*a17);
  a81=(a81*a65);
  a79=(a79+a81);
  a84=(a84*a18);
  a79=(a79+a84);
  a84=(a44*a90);
  a79=(a79+a84);
  a84=(a87*a70);
  a84=(a84*a9);
  a81=(a87*a20);
  a81=(a81*a2);
  a84=(a84-a81);
  a81=(a16*a91);
  a2=(a70*a90);
  a81=(a81-a2);
  a70=(a70*a81);
  a47=(a20*a47);
  a70=(a70-a47);
  a84=(a84-a70);
  a21=(a21*a88);
  a10=(a10*a69);
  a11=(a11*a66);
  a10=(a10+a11);
  a82=(a82*a68);
  a10=(a10+a82);
  a82=(a10*a91);
  a21=(a21-a82);
  a84=(a84-a21);
  a79=(a79+a84);
  a79=(a0*a79);
  if (res[0]!=0) res[0][11]=a79;
  a78=(a78*a17);
  a80=(a80*a65);
  a78=(a78+a80);
  a35=(a35*a18);
  a78=(a78+a35);
  a44=(a44*a91);
  a78=(a78+a44);
  a44=(a87*a20);
  a44=(a44*a12);
  a87=(a87*a16);
  a87=(a87*a9);
  a44=(a44-a87);
  a20=(a20*a73);
  a16=(a16*a81);
  a20=(a20-a16);
  a44=(a44-a20);
  a10=(a10*a90);
  a52=(a52*a88);
  a10=(a10-a52);
  a44=(a44-a10);
  a78=(a78+a44);
  a0=(a0*a78);
  if (res[0]!=0) res[0][12]=a0;
  return 0;
}

CASADI_SYMBOL_EXPORT int fdrift(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem){
  return casadi_f0(arg, res, iw, w, mem);
}

CASADI_SYMBOL_EXPORT int fdrift_alloc_mem(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT int fdrift_init_mem(int mem) {
  return 0;
}

CASADI_SYMBOL_EXPORT void fdrift_free_mem(int mem) {
}

CASADI_SYMBOL_EXPORT int fdrift_checkout(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT void fdrift_release(int mem) {
}

CASADI_SYMBOL_EXPORT void fdrift_incref(void) {
}

CASADI_SYMBOL_EXPORT void fdrift_decref(void) {
}

CASADI_SYMBOL_EXPORT casadi_int fdrift_n_in(void) { return 3;}

CASADI_SYMBOL_EXPORT casadi_int fdrift_n_out(void) { return 1;}

CASADI_SYMBOL_EXPORT casadi_real fdrift_default_in(casadi_int i){
  switch (i) {
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* fdrift_name_in(casadi_int i){
  switch (i) {
    case 0: return "i0";
    case 1: return "i1";
    case 2: return "i2";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* fdrift_name_out(casadi_int i){
  switch (i) {
    case 0: return "o0";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* fdrift_sparsity_in(casadi_int i) {
  switch (i) {
    case 0: return casadi_s0;
    case 1: return casadi_s1;
    case 2: return casadi_s2;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* fdrift_sparsity_out(casadi_int i) {
  switch (i) {
    case 0: return casadi_s1;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT int fdrift_work(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w) {
  if (sz_arg) *sz_arg = 3;
  if (sz_res) *sz_res = 1;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 0;
  return 0;
}

#ifdef MATLAB_MEX_FILE
void mex_fdrift(int resc, mxArray *resv[], int argc, const mxArray *argv[]) {
  casadi_int i;
  casadi_real w[127];
  casadi_int *iw = 0;
  const casadi_real* arg[3] = {0};
  casadi_real* res[1] = {0};
  if (argc>3) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"fdrift\" failed. Too many input arguments (%d, max 3)", argc);
  if (resc>1) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"fdrift\" failed. Too many output arguments (%d, max 1)", resc);
  if (--argc>=0) arg[0] = casadi_from_mex(argv[0], w, casadi_s0, w+33);
  if (--argc>=0) arg[1] = casadi_from_mex(argv[1], w+1, casadi_s1, w+33);
  if (--argc>=0) arg[2] = casadi_from_mex(argv[2], w+14, casadi_s2, w+33);
  --resc;
  res[0] = w+20;
  i = fdrift(arg, res, iw, w+33, 0);
  if (i) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"fdrift\" failed.");
  if (res[0]) resv[0] = casadi_to_mex(casadi_s1, res[0]);
}
#endif

casadi_int main_fdrift(casadi_int argc, char* argv[]) {
  casadi_int j;
  casadi_real* a;
  const casadi_real* r;
  casadi_int flag;
  casadi_int *iw = 0;
  casadi_real w[127];
  const casadi_real* arg[3];
  casadi_real* res[1];
  arg[0] = w+0;
  arg[1] = w+1;
  arg[2] = w+14;
  res[0] = w+20;
  a = w;
  for (j=0; j<20; ++j) if (scanf("%lg", a++)<=0) return 2;
  flag = fdrift(arg, res, iw, w+33, 0);
  if (flag) return flag;
  r = w+20;
  for (j=0; j<13; ++j) CASADI_PRINTF("%g ", *r++);
  CASADI_PRINTF("\n");
  return 0;
}


#ifdef MATLAB_MEX_FILE
void mexFunction(int resc, mxArray *resv[], int argc, const mxArray *argv[]) {
  char buf[7];
  int buf_ok = argc > 0 && !mxGetString(*argv, buf, sizeof(buf));
  if (!buf_ok) {
    mex_fdrift(resc, resv, argc, argv);
    return;
  } else if (strcmp(buf, "fdrift")==0) {
    mex_fdrift(resc, resv, argc-1, argv+1);
    return;
  }
  mexErrMsgTxt("First input should be a command string. Possible values: 'fdrift'");
}
#endif
int main(int argc, char* argv[]) {
  if (argc<2) {
    /* name error */
  } else if (strcmp(argv[1], "fdrift")==0) {
    return main_fdrift(argc-2, argv+2);
  }
  fprintf(stderr, "First input should be a command string. Possible values: 'fdrift'\nNote: you may use function.generate_input to create a command string.\n");
  return 1;
}
#ifdef __cplusplus
} /* extern "C" */
#endif
