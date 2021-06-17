/* WARNING: THIS FUNCTION PERFORMS NO ERROR CHECKING WHATSOEVER! Don't even
 * think of calling sdde_mex itself! It is designed to be called via the
 * sdde_mex.m Matlab wrapper function (utils subfolder); all error checking and
 * documentation happens there. To compile (Linux):

mex -O -R2018a CFLAGS="\$CFLAGS -std=c99 -march=native -O3 -flto -Wall -Wextra -Wconversion -Winline -pedantic-errors -D_POSIX_C_SOURCE=199309L -D_BSD_SOURCE -DLINUX -static" -outdir ../mex sdde_mex.c

*/

#include "mex.h"
#include "matrix.h"
#include <string.h> // for memcpy
#include <inttypes.h>

#define ct_assert(e) extern char (*ct_assert(void)) [sizeof(char[1 - 2*!(e)])]

#define UNUSED __attribute__ ((unused))

#define FUINT_T uint_fast32_t
#define MXGETFUINTS mxGetUint64s

void mexFunction(int UNUSED nlhs, mxArray *plhs[], int UNUSED nrhs, const mxArray *prhs[])
{
	ct_assert(sizeof(FUINT_T) == sizeof(uint64_t)); // must match int type in sdde.m !!!

	const double*  const a = mxGetPr(prhs[0]);
	const FUINT_T* const i = MXGETFUINTS(prhs[1]);
	const FUINT_T* const j = MXGETFUINTS(prhs[2]);
	const double*  const s = mxGetPr(prhs[3]);
	const FUINT_T* const d = MXGETFUINTS(prhs[4]);
	const double*  const x = mxGetPr(prhs[5]);

	const FUINT_T c = mxGetNumberOfElements(prhs[1]);
	const FUINT_T n = mxGetM(prhs[5]);
	const FUINT_T m = mxGetN(prhs[5]);

	// mexPrintf("\nfastest 32 bit unsigned type (bytes) : %u\n",sizeof(uint_fast32_t));
	// mexPrintf("fastest 64 bit unsigned type (bytes) : %u\n",sizeof(uint_fast64_t));
	// mexPrintf("size_t type (bytes)                  : %u\n",sizeof(size_t));

	double* const y = mxGetPr(plhs[0] = mxCreateDoubleMatrix(n,m,mxREAL)); // create output array

	memcpy(y,x,n*m*sizeof(double)); // copy input x to output y

	for (FUINT_T t = 1; t < m; ++t) {
		double* const yt = y+n*t;
		for (FUINT_T r=0; r<n; ++r) *(yt+r) += a[r]*(*(yt-n+r));
		for (FUINT_T k=0; k<c; ++k) if (t >= d[k]) *(yt+i[k]) += s[k]*(*(yt-n*d[k]+j[k]));
	}
}
