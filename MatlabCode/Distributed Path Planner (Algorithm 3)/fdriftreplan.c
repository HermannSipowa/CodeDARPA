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
  #define CASADI_PREFIX(ID) fdriftreplan_ ## ID
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
static const casadi_int casadi_s1[10] = {6, 1, 0, 6, 0, 1, 2, 3, 4, 5};

/* fdriftreplan:(i0,i1[6],i2[6])->(o0[6]) */
static int casadi_f0(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem) {
  casadi_real a0, a1, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a2, a20, a21, a22, a23, a24, a25, a26, a27, a28, a29, a3, a30, a31, a32, a33, a34, a35, a36, a4, a5, a6, a7, a8, a9;
  a0=1.6840054676261389e+04;
  a1=arg[1]? arg[1][3] : 0;
  a2=(a0*a1);
  if (res[0]!=0) res[0][0]=a2;
  a2=arg[1]? arg[1][4] : 0;
  a3=(a0*a2);
  if (res[0]!=0) res[0][1]=a3;
  a3=arg[1]? arg[1][5] : 0;
  a4=(a0*a3);
  if (res[0]!=0) res[0][2]=a4;
  a4=-3.9860044150000002e+05;
  a5=arg[2]? arg[2][0] : 0;
  a6=casadi_sq(a5);
  a7=arg[2]? arg[2][1] : 0;
  a8=casadi_sq(a7);
  a6=(a6+a8);
  a8=arg[2]? arg[2][2] : 0;
  a9=casadi_sq(a8);
  a6=(a6+a9);
  a9=sqrt(a6);
  a10=arg[1]? arg[1][0] : 0;
  a9=(a9+a10);
  a11=casadi_sq(a9);
  a12=arg[1]? arg[1][1] : 0;
  a13=casadi_sq(a12);
  a11=(a11+a13);
  a13=arg[1]? arg[1][2] : 0;
  a14=casadi_sq(a13);
  a11=(a11+a14);
  a14=1.5000000000000000e+00;
  a11=pow(a11,a14);
  a4=(a4/a11);
  a9=(a4*a9);
  a11=3.9860044150000002e+05;
  a11=(a11/a6);
  a9=(a9+a11);
  a11=-2.;
  a14=arg[2]? arg[2][3] : 0;
  a15=(a8*a14);
  a16=arg[2]? arg[2][5] : 0;
  a17=(a5*a16);
  a15=(a15-a17);
  a17=(a7*a16);
  a18=arg[2]? arg[2][4] : 0;
  a19=(a8*a18);
  a17=(a17-a19);
  a19=casadi_sq(a17);
  a20=casadi_sq(a15);
  a19=(a19+a20);
  a20=(a5*a18);
  a21=(a7*a14);
  a20=(a20-a21);
  a21=casadi_sq(a20);
  a19=(a19+a21);
  a19=sqrt(a19);
  a15=(a15/a19);
  a21=casadi_sq(a5);
  a22=casadi_sq(a7);
  a21=(a21+a22);
  a22=casadi_sq(a8);
  a21=(a21+a22);
  a21=sqrt(a21);
  a22=(a8/a21);
  a23=(a15*a22);
  a20=(a20/a19);
  a24=(a7/a21);
  a25=(a20*a24);
  a23=(a23-a25);
  a25=casadi_sq(a23);
  a21=(a5/a21);
  a26=(a20*a21);
  a17=(a17/a19);
  a19=(a17*a22);
  a26=(a26-a19);
  a19=casadi_sq(a26);
  a25=(a25+a19);
  a19=(a17*a24);
  a27=(a15*a21);
  a19=(a19-a27);
  a27=casadi_sq(a19);
  a25=(a25+a27);
  a25=sqrt(a25);
  a23=(a23/a25);
  a27=(a7*a16);
  a28=(a8*a18);
  a27=(a27-a28);
  a28=(a27/a6);
  a29=(a23*a28);
  a26=(a26/a25);
  a30=(a8*a14);
  a31=(a5*a16);
  a30=(a30-a31);
  a31=(a30/a6);
  a32=(a26*a31);
  a29=(a29+a32);
  a19=(a19/a25);
  a25=(a5*a18);
  a32=(a7*a14);
  a25=(a25-a32);
  a32=(a25/a6);
  a33=(a19*a32);
  a29=(a29+a33);
  a33=(a11*a29);
  a33=(a33*a3);
  a34=(a17*a28);
  a35=(a15*a31);
  a34=(a34+a35);
  a35=(a20*a32);
  a34=(a34+a35);
  a35=(a11*a34);
  a35=(a35*a2);
  a33=(a33-a35);
  a28=(a21*a28);
  a31=(a24*a31);
  a28=(a28+a31);
  a32=(a22*a32);
  a28=(a28+a32);
  a32=(a28*a12);
  a31=(a29*a10);
  a32=(a32-a31);
  a31=(a29*a32);
  a35=(a34*a10);
  a36=(a28*a13);
  a35=(a35-a36);
  a36=(a34*a35);
  a31=(a31-a36);
  a33=(a33-a31);
  a5=(a5*a14);
  a7=(a7*a18);
  a5=(a5+a7);
  a8=(a8*a16);
  a5=(a5+a8);
  a5=(a11*a5);
  a27=(a5*a27);
  a6=casadi_sq(a6);
  a27=(a27/a6);
  a23=(a23*a27);
  a30=(a5*a30);
  a30=(a30/a6);
  a26=(a26*a30);
  a23=(a23+a26);
  a5=(a5*a25);
  a5=(a5/a6);
  a19=(a19*a5);
  a23=(a23+a19);
  a19=(a23*a13);
  a17=(a17*a27);
  a15=(a15*a30);
  a17=(a17+a15);
  a20=(a20*a5);
  a17=(a17+a20);
  a20=(a17*a12);
  a19=(a19-a20);
  a33=(a33-a19);
  a9=(a9+a33);
  a9=(a0*a9);
  if (res[0]!=0) res[0][3]=a9;
  a9=(a4*a12);
  a33=(a11*a34);
  a33=(a33*a1);
  a19=(a11*a28);
  a19=(a19*a3);
  a33=(a33-a19);
  a19=(a29*a13);
  a3=(a34*a12);
  a19=(a19-a3);
  a34=(a34*a19);
  a32=(a28*a32);
  a34=(a34-a32);
  a33=(a33-a34);
  a17=(a17*a10);
  a21=(a21*a27);
  a24=(a24*a30);
  a21=(a21+a24);
  a22=(a22*a5);
  a21=(a21+a22);
  a22=(a21*a13);
  a17=(a17-a22);
  a33=(a33-a17);
  a9=(a9+a33);
  a9=(a0*a9);
  if (res[0]!=0) res[0][4]=a9;
  a4=(a4*a13);
  a13=(a11*a28);
  a13=(a13*a2);
  a11=(a11*a29);
  a11=(a11*a1);
  a13=(a13-a11);
  a28=(a28*a35);
  a29=(a29*a19);
  a28=(a28-a29);
  a13=(a13-a28);
  a21=(a21*a12);
  a23=(a23*a10);
  a21=(a21-a23);
  a13=(a13-a21);
  a4=(a4+a13);
  a0=(a0*a4);
  if (res[0]!=0) res[0][5]=a0;
  return 0;
}

CASADI_SYMBOL_EXPORT int fdriftreplan(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem){
  return casadi_f0(arg, res, iw, w, mem);
}

CASADI_SYMBOL_EXPORT int fdriftreplan_alloc_mem(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT int fdriftreplan_init_mem(int mem) {
  return 0;
}

CASADI_SYMBOL_EXPORT void fdriftreplan_free_mem(int mem) {
}

CASADI_SYMBOL_EXPORT int fdriftreplan_checkout(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT void fdriftreplan_release(int mem) {
}

CASADI_SYMBOL_EXPORT void fdriftreplan_incref(void) {
}

CASADI_SYMBOL_EXPORT void fdriftreplan_decref(void) {
}

CASADI_SYMBOL_EXPORT casadi_int fdriftreplan_n_in(void) { return 3;}

CASADI_SYMBOL_EXPORT casadi_int fdriftreplan_n_out(void) { return 1;}

CASADI_SYMBOL_EXPORT casadi_real fdriftreplan_default_in(casadi_int i){
  switch (i) {
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* fdriftreplan_name_in(casadi_int i){
  switch (i) {
    case 0: return "i0";
    case 1: return "i1";
    case 2: return "i2";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* fdriftreplan_name_out(casadi_int i){
  switch (i) {
    case 0: return "o0";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* fdriftreplan_sparsity_in(casadi_int i) {
  switch (i) {
    case 0: return casadi_s0;
    case 1: return casadi_s1;
    case 2: return casadi_s1;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* fdriftreplan_sparsity_out(casadi_int i) {
  switch (i) {
    case 0: return casadi_s1;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT int fdriftreplan_work(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w) {
  if (sz_arg) *sz_arg = 3;
  if (sz_res) *sz_res = 1;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 0;
  return 0;
}

#ifdef MATLAB_MEX_FILE
void mex_fdriftreplan(int resc, mxArray *resv[], int argc, const mxArray *argv[]) {
  casadi_int i;
  casadi_real w[56];
  casadi_int *iw = 0;
  const casadi_real* arg[3] = {0};
  casadi_real* res[1] = {0};
  if (argc>3) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"fdriftreplan\" failed. Too many input arguments (%d, max 3)", argc);
  if (resc>1) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"fdriftreplan\" failed. Too many output arguments (%d, max 1)", resc);
  if (--argc>=0) arg[0] = casadi_from_mex(argv[0], w, casadi_s0, w+19);
  if (--argc>=0) arg[1] = casadi_from_mex(argv[1], w+1, casadi_s1, w+19);
  if (--argc>=0) arg[2] = casadi_from_mex(argv[2], w+7, casadi_s1, w+19);
  --resc;
  res[0] = w+13;
  i = fdriftreplan(arg, res, iw, w+19, 0);
  if (i) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"fdriftreplan\" failed.");
  if (res[0]) resv[0] = casadi_to_mex(casadi_s1, res[0]);
}
#endif

casadi_int main_fdriftreplan(casadi_int argc, char* argv[]) {
  casadi_int j;
  casadi_real* a;
  const casadi_real* r;
  casadi_int flag;
  casadi_int *iw = 0;
  casadi_real w[56];
  const casadi_real* arg[3];
  casadi_real* res[1];
  arg[0] = w+0;
  arg[1] = w+1;
  arg[2] = w+7;
  res[0] = w+13;
  a = w;
  for (j=0; j<13; ++j) if (scanf("%lg", a++)<=0) return 2;
  flag = fdriftreplan(arg, res, iw, w+19, 0);
  if (flag) return flag;
  r = w+13;
  for (j=0; j<6; ++j) CASADI_PRINTF("%g ", *r++);
  CASADI_PRINTF("\n");
  return 0;
}


#ifdef MATLAB_MEX_FILE
void mexFunction(int resc, mxArray *resv[], int argc, const mxArray *argv[]) {
  char buf[13];
  int buf_ok = argc > 0 && !mxGetString(*argv, buf, sizeof(buf));
  if (!buf_ok) {
    mex_fdriftreplan(resc, resv, argc, argv);
    return;
  } else if (strcmp(buf, "fdriftreplan")==0) {
    mex_fdriftreplan(resc, resv, argc-1, argv+1);
    return;
  }
  mexErrMsgTxt("First input should be a command string. Possible values: 'fdriftreplan'");
}
#endif
int main(int argc, char* argv[]) {
  if (argc<2) {
    /* name error */
  } else if (strcmp(argv[1], "fdriftreplan")==0) {
    return main_fdriftreplan(argc-2, argv+2);
  }
  fprintf(stderr, "First input should be a command string. Possible values: 'fdriftreplan'\nNote: you may use function.generate_input to create a command string.\n");
  return 1;
}
#ifdef __cplusplus
} /* extern "C" */
#endif
