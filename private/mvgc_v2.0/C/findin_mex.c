/* WARNING: THIS FUNCTION PERFORMS NO ERROR CHECKING WHATSOEVER! Don't even
 * think of calling findin_mex itself! It is designed to be called via the
 * findin.m Matlab wrapper function (utils subfolder); all error checking and
 * documentation happens there. To compile, see Makefile in this folder */

#include "mex.h"

#define UNUSED __attribute__ ((unused))

void mexFunction(int UNUSED nlhs, mxArray *plhs[], int UNUSED nrhs, const mxArray *prhs[])
{
	const double* const a  = mxGetPr(prhs[0]);
	const double* const b  = mxGetPr(prhs[1]);

	const mwSize na = mxGetNumberOfElements(prhs[0]);
	const mwSize nb = mxGetNumberOfElements(prhs[1]);

	double* const idx = mxGetPr(plhs[0] = mxCreateDoubleMatrix(1,na,mxREAL));

	for (size_t i=0;i<na;++i) {
		const double ai = a[i];
		for (size_t j=0;j<nb;++j) {
			if (ai == b[j]) {
				idx[i] = (double)(j+1);
			}
		}
	}
}
