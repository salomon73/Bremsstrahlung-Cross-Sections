#include "mex.h"
#include <complex>
#include <vector>
#include <flint/acb.h>
#include <flint/acb_hypgeom.h>

typedef std::complex<double> dcomplex;

// Vectorized version
void hypergeometric2F1_flint_vec(const std::vector<dcomplex>& a,
                                 const std::vector<dcomplex>& b,
                                 const std::vector<dcomplex>& c,
                                 const std::vector<dcomplex>& z,
                                 std::vector<dcomplex>& out)
{
    size_t n = a.size();
    out.resize(n);

    for (size_t i = 0; i < n; ++i)
    {
        out[i] = hypergeometric2F1_mex(a[i], b[i], c[i], z[i]);
    }
}

// MEX entry point
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 4)
        mexErrMsgTxt("Usage: y = hypergeometric2f1_mex_vec(a,b,c,z)");

    size_t n = mxGetNumberOfElements(prhs[0]);
    std::vector<dcomplex> a(n), b(n), c(n), z(n);

    double* a_pr = mxGetPr(prhs[0]); double* a_pi = mxGetPi(prhs[0]);
    double* b_pr = mxGetPr(prhs[1]); double* b_pi = mxGetPi(prhs[1]);
    double* c_pr = mxGetPr(prhs[2]); double* c_pi = mxGetPi(prhs[2]);
    double* z_pr = mxGetPr(prhs[3]); double* z_pi = mxGetPi(prhs[3]);

    for (size_t i = 0; i < n; ++i)
    {
        a[i] = dcomplex(a_pr[i], a_pi ? a_pi[i] : 0.0);
        b[i] = dcomplex(b_pr[i], b_pi ? b_pi[i] : 0.0);
        c[i] = dcomplex(c_pr[i], c_pi ? c_pi[i] : 0.0);
        z[i] = dcomplex(z_pr[i], z_pi ? z_pi[i] : 0.0);
    }

    std::vector<dcomplex> out;
    hypergeometric2F1_flint_vec(a, b, c, z, out);

    plhs[0] = mxCreateDoubleMatrix(n, 1, mxCOMPLEX);
    double* pr = mxGetPr(plhs[0]);
    double* pi = mxGetPi(plhs[0]);
    for (size_t i = 0; i < n; ++i)
    {
        pr[i] = out[i].real();
        pi[i] = out[i].imag();
    }
}
