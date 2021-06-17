/* WARNING: THIS FUNCTION PERFORMS NO ERROR CHECKING WHATSOEVER! Don't even
 * think of calling genvar_mex itself! It is designed to be called via the
 * genvma.m Matlab wrapper function (utils subfolder); all error checking and
 * documentation happens there. To compile, see Makefile in this folder */

#include "mex.h"
#include <string.h> // for memcpy

#define UNUSED __attribute__ ((unused))

void mexFunction(int UNUSED nlhs, mxArray *plhs[], int UNUSED nrhs, const mxArray *prhs[])
{
	const mxArray* const pc = prhs[0];
	const mxArray* const px = prhs[1];

	const double* const c = mxGetPr(pc);
	const double* const x = mxGetPr(px);

	const size_t n  = mxGetM(px);
	const size_t m  = mxGetN(px);
	const mxLogical cvec = *mxGetLogicals(prhs[2]);
	const mwSize p  = cvec ? mxGetNumberOfElements(pc) : mxGetNumberOfDimensions(pc) < 3 ? 1 : mxGetDimensions(pc)[2];

	double* const y = mxGetPr(plhs[0] = mxCreateDoubleMatrix(n,m,mxREAL)); /* allocate output array y  */

	const size_t mmp = m < p ? m : p;
	const size_t nsq = n*n;

	memcpy(y,x,n*m*sizeof(double)); /* copy input x to output y */

	if (cvec) { /* apply c to each row in turn */
		if (n == 1) { /* y is c vector: special case */
			for (size_t t=0; t<mmp; ++t) {
				double yvar = 0.0;
				for (size_t k=0; k<t; ++k) {
					yvar += c[k]*x[t-(k+1)]; /* lag is k+1 !!! */
				}
				y[t] += yvar;
			}
			for (size_t t=mmp; t<m; ++t) {
				double yvar = 0.0;
				for (size_t k=0; k<p; ++k) {
					yvar += c[k]*x[t-(k+1)]; /* lag is k+1 !!! */
				}
				y[t] += yvar;
			}
		}
		else {
			for (size_t i=0; i<n; ++i) {
				const double* const xi = x+i;
				double* const yi = y+i;
				for (size_t t=0; t<mmp; ++t) {
					double yvar = 0.0;
					for (size_t k=0; k<t; ++k) {
						yvar += c[k]*xi[n*(t-(k+1))]; /* lag is k+1 !!! */
					}
					yi[n*t] += yvar;
				}
				for (size_t t=mmp; t<m; ++t) {
					double yvar = 0.0;
					for (size_t k=0; k<p; ++k) {
						yvar += c[k]*xi[n*(t-(k+1))]; /* lag is k+1 !!! */
					}
					yi[n*t] += yvar;
				}
			}
		}
	}
	else {
		for (size_t t=0; t<mmp; ++t) {
			const double* const xt = x+n*t;
			double* const yt = y+n*t;
			for (size_t k=0; k<t; ++k) {
				const double* const ck = c+nsq*k;
				const double* const xlag = xt-n*(k+1); /* lag is k+1 !!! */
				for (size_t i=0; i<n; ++i) {
					const double* const cki = ck+i;
					double yvar = 0.0;
					for (size_t j=0; j<n; ++j) {
						yvar += cki[n*j]*xlag[j];
					}
					yt[i] += yvar;
				}
			}
		}
		for (size_t t=mmp; t<m; ++t) {
			const double* const xt = x+n*t;
			double* const yt = y+n*t;
			for (size_t k=0; k<p; ++k) {
				const double* const ck = c+nsq*k;
				const double* const xlag = xt-n*(k+1); /* lag is k+1 !!! */
				for (size_t i=0; i<n; ++i) {
					const double* const cki = ck+i;
					double yvar = 0.0;
					for (size_t j=0; j<n; ++j) {
						yvar += cki[n*j]*xlag[j];
					}
					yt[i] += yvar;
				}
			}
		}
	}
}
