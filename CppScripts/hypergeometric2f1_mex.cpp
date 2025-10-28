#include "mex.h"
#include <complex>
#include <flint/acb.h>
#include <flint/acb_hypgeom.h>

typedef std::complex<double> dcomplex;

dcomplex hypergeometric2F1_flint(dcomplex a, dcomplex b, dcomplex c, dcomplex z, slong prec = 64)
{
    acb_t A, B, C, Z, R;
    acb_init(A); acb_init(B); acb_init(C); acb_init(Z); acb_init(R);

    acb_set_d_d(A, a.real(), a.imag());
    acb_set_d_d(B, b.real(), b.imag());
    acb_set_d_d(C, c.real(), c.imag());
    acb_set_d_d(Z, z.real(), z.imag());

    // Compute 2F1
    acb_hypgeom_2f1(R, A, B, C, Z, 0, prec);

    // Correct way: use arb_midref
    double real_res = arf_get_d(arb_midref(acb_realref(R)), ARF_RND_NEAR);
    double imag_res = arf_get_d(arb_midref(acb_imagref(R)), ARF_RND_NEAR);

    acb_clear(A); acb_clear(B); acb_clear(C); acb_clear(Z); acb_clear(R);

    return dcomplex(real_res, imag_res);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs < 4)
        mexErrMsgTxt("Usage: y = hypergeometric2f1_mex(a,b,c,z)");

    dcomplex a(mxGetScalar(prhs[0]),0);
    dcomplex b(mxGetScalar(prhs[1]),0);
    dcomplex c(mxGetScalar(prhs[2]),0);
    dcomplex z(mxGetScalar(prhs[3]),0);

    dcomplex res = hypergeometric2F1_flint(a,b,c,z);

    plhs[0] = mxCreateDoubleMatrix(1,1,mxCOMPLEX);
    *mxGetPr(plhs[0]) = res.real();
    *mxGetPi(plhs[0]) = res.imag();
}
