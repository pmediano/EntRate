/* WARNING: THIS FUNCTION PERFORMS NO ERROR CHECKING WHATSOEVER! Don't even
 * think of calling mvfilter_mex itself! It is designed to be called via the
 * mvfilter.m Matlab wrapper function (utils subfolder); all error checking and
 * documentation happens there. To compile, see Makefile in this folder */

#include "mex.h"
#include <string.h> // for memcpy

#define UNUSED __attribute__ ((unused))

#define MIN(a,b) ((a) < (b) ? (a) : (b))

void mexFunction(int UNUSED nlhs, mxArray *plhs[], int UNUSED nrhs, const mxArray *prhs[])
{
	const double* const b = mxGetPr(prhs[0]);
	const bool bvec = (bool)(*mxGetLogicals(prhs[1]));
	const size_t q = (size_t)mxGetScalar(prhs[2]);

	const double* const a = mxGetPr(prhs[3]);
	const bool avec = (bool)(*mxGetLogicals(prhs[4]));
	const size_t p = (size_t)mxGetScalar(prhs[5]);

	const double* const x = mxGetPr(prhs[6]);

	const size_t n = mxGetM(prhs[6]);
	const size_t m = mxGetN(prhs[6]);

	double* const y = mxGetPr(plhs[0] = mxCreateDoubleMatrix(n,m,mxREAL)); /* allocate output array y  */

	const size_t nsq = n*n;

	memcpy(y,x,n*m*sizeof(double)); // copy input x to output y

	if (q > 0) { // MA component
		if (p > 0) { // AR component
			if (bvec) {
				if (avec) { // bvec, avec
					for (size_t t=0; t<m; ++t) {
						const double* const xt = x+n*t;
						double*       const yt = y+n*t;
						for (size_t k=0; k<MIN(q,t); ++k) {
							const double* const xtlagk = xt-n*(k+1); // lag is k+1 !!!
							for (size_t i=0; i<n; ++i) yt[i] += b[k]*xtlagk[i];
						}
						for (size_t k=0; k<MIN(p,t); ++k) {
							const double* const ytlagk = yt-n*(k+1); // lag is k+1 !!!
							for (size_t i=0; i<n; ++i) yt[i] += a[k]*ytlagk[i];
						}
					}
				}
				else { // bvec, ~avec
					for (size_t t=0; t<m; ++t) {
						const double* const xt = x+n*t;
						double*       const yt = y+n*t;
						for (size_t k=0; k<MIN(q,t); ++k) {
							const double* const xtlagk = xt-n*(k+1); // lag is k+1 !!!
							for (size_t i=0; i<n; ++i) yt[i] += b[k]*xtlagk[i];
						}
						for (size_t k=0; k<MIN(p,t); ++k) {
							const double* const ak = a+nsq*k;
							const double* const ytlagk = yt-n*(k+1); // lag is k+1 !!!
							for (size_t i=0; i<n; ++i) for (size_t j=0; j<n; ++j) yt[i] += ak[i+n*j]*ytlagk[j];
						}
					}
				}
			}
			else { // ~bvec
				if (avec) { // ~bvec, avec
					for (size_t t=0; t<m; ++t) {
						const double* const xt = x+n*t;
						double*       const yt = y+n*t;
						for (size_t k=0; k<MIN(q,t); ++k) {
							const double* const bk = b+nsq*k;
							const double* const xtlagk = xt-n*(k+1); // lag is k+1 !!!
							for (size_t i=0; i<n; ++i) for (size_t j=0; j<n; ++j) yt[i] += bk[i+n*j]*xtlagk[j];
						}
						for (size_t k=0; k<MIN(p,t); ++k) {
							const double* const ytlagk = yt-n*(k+1); // lag is k+1 !!!
							for (size_t i=0; i<n; ++i) yt[i] += a[k]*ytlagk[i];
						}
					}
				}
				else { // ~bvec, ~avec
					for (size_t t=0; t<m; ++t) {
						const double* const xt = x+n*t;
						double*       const yt = y+n*t;
						for (size_t k=0; k<MIN(q,t); ++k) {
							const double* const bk = b+nsq*k;
							const double* const xtlagk = xt-n*(k+1); // lag is k+1 !!!
							for (size_t i=0; i<n; ++i) for (size_t j=0; j<n; ++j) yt[i] += bk[i+n*j]*xtlagk[j];
						}
						for (size_t k=0; k<MIN(p,t); ++k) {
							const double* const ak = a+nsq*k;
							const double* const ytlagk = yt-n*(k+1); // lag is k+1 !!!
							for (size_t i=0; i<n; ++i) for (size_t j=0; j<n; ++j) yt[i] += ak[i+n*j]*ytlagk[j];
						}
					}
				}
			}
		}
		else { // p == 0 : pure MA
			const size_t mmq = m < q ? m : q;
			if (bvec) {
				for (size_t t=0; t<mmq; ++t) {
					const double* const xt = x+n*t;
					double*       const yt = y+n*t;
					for (size_t k=0; k<t; ++k) {
						const double* const xtlagk = xt-n*(k+1); // lag is k+1 !!!
						for (size_t i=0; i<n; ++i) yt[i] += b[k]*xtlagk[i];
					}
				}
				for (size_t t=mmq; t<m; ++t) {
					const double* const xt = x+n*t;
					double*       const yt = y+n*t;
					for (size_t k=0; k<q; ++k) {
						const double* const xtlagk = xt-n*(k+1); // lag is k+1 !!!
						for (size_t i=0; i<n; ++i) yt[i] += b[k]*xtlagk[i];
					}
				}
			}
			else { // ~bvec
				for (size_t t=0; t<mmq; ++t) {
					const double* const xt = x+n*t;
					double*       const yt = y+n*t;
					for (size_t k=0; k<t; ++k) {
						const double* const bk = b+nsq*k;
						const double* const xtlagk = xt-n*(k+1); // lag is k+1 !!!
						for (size_t i=0; i<n; ++i) for (size_t j=0; j<n; ++j) yt[i] += bk[i+n*j]*xtlagk[j];
					}
				}
				for (size_t t=mmq; t<m; ++t) {
					const double* const xt = x+n*t;
					double*       const yt = y+n*t;
					for (size_t k=0; k<q; ++k) {
						const double* const bk = b+nsq*k;
						const double* const xtlagk = xt-n*(k+1); // lag is k+1 !!!
						for (size_t i=0; i<n; ++i) for (size_t j=0; j<n; ++j) yt[i] += bk[i+n*j]*xtlagk[j];
					}
				}
			}
		}
	}
	else { // q == 0 : pure AR
		const size_t mmp = m < p ? m : p;
		if (avec) {
			for (size_t t=0; t<mmp; ++t) {
				double* const yt = y+n*t;
				for (size_t k=0; k<t; ++k) {
					const double* const ytlagk = yt-n*(k+1); // lag is k+1 !!!
					for (size_t i=0; i<n; ++i) yt[i] += a[k]*ytlagk[i];
				}
			}
			for (size_t t=mmp; t<m; ++t) {
				double* const yt = y+n*t;
				for (size_t k=0; k<p; ++k) {
					const double* const ytlagk = yt-n*(k+1); // lag is k+1 !!!
					for (size_t i=0; i<n; ++i) yt[i] += a[k]*ytlagk[i];
				}
			}
		}
		else { // ~avec
			for (size_t t=0; t<mmp; ++t) {
				double* const yt = y+n*t;
				for (size_t k=0; k<t; ++k) {
					const double* const ak = a+nsq*k;
					const double* const ytlagk = yt-n*(k+1); // lag is k+1 !!!
					for (size_t i=0; i<n; ++i) for (size_t j=0; j<n; ++j) yt[i] += ak[i+n*j]*ytlagk[j];
				}
			}
			for (size_t t=mmp; t<m; ++t) {
				double* const yt = y+n*t;
				for (size_t k=0; k<p; ++k) {
					const double* const ak = a+nsq*k;
					const double* const ytlagk = yt-n*(k+1); // lag is k+1 !!!
					for (size_t i=0; i<n; ++i) for (size_t j=0; j<n; ++j) yt[i] += ak[i+n*j]*ytlagk[j];
				}
			}
		}
	}
}

#undef MIN
