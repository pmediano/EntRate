#include "mex.h"

/* WARNING: THIS FUNCTION PERFORMS VIRTUALLY NO ERROR CHECKING!
   To compile, see mvgc_makemex.m in the utils subfolder */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double* kt        = mxGetPr(prhs[0]);
	double* rr        = mxGetPr(prhs[1]);
	double* yy        = mxGetPr(prhs[2]);
	const size_t nr   = mxGetN(prhs[1]);
	const size_t ny   = mxGetN(prhs[2]);
	const size_t pny  = mxGetM(prhs[0]);
	const size_t nr1  = mxGetN(prhs[0]);

	const size_t p    = pny/ny;
	const size_t n    = nr+ny;
	const size_t pn   = p*n;

	if (nr1  != nr ) mexErrMsgTxt("size mismatch (nr)"  );
	if (p*ny != pny) mexErrMsgTxt("size mismatch (p,ny)");

	double* kr = mxGetPr(plhs[0] = mxCreateDoubleMatrix(pn,nr,mxREAL));

	size_t* r = mxCalloc(nr,sizeof(size_t));
	size_t* y = mxCalloc(ny,sizeof(size_t));

	size_t i,j,k,l;

	for (i=0; i<nr; ++i) r[i] = (size_t)rr[i]-1;
	for (j=0; j<ny; ++j) y[j] = (size_t)yy[j]-1;

	for (i=0; i<nr; ++i) {
		kr[r[i]+pn*i] = 1.0;
	}

	l = 0;
	for (k=0; k<pny; k += ny) {
		for (i=0; i<nr; ++i) {
			for (j=0; j<ny; ++j) {
				kr[l+y[j]+pn*i] = kt[k+j+pny*i];
			}
		}
		l += n;
	}

	mxFree(r);
	mxFree(y);
}
