/**
 * LZ76 complexity algorithm for Matlab's MEX. Compile using
 *
 * mex LZ76.cpp
 *
 * from the Matlab prompt. The function expects a 1D array of logicals,
 * such as
 *
 * LZ76(rand([100,1]) < 0.5)
 *
 * Authors: Pedro Mediano
 *          Fernando Rosas
 *          Joel Dyer
 */
#include "mex.h"

double LZc(bool fSS[], int n) {

    /* Computes the LZ76 complexity index of binary sequence fSS using the
       algorithm described in F. Kaspar, H. G. Schuster, "Easily-calculable 
       measure for the complexity of spatiotemporal patterns", Physical Review
       A, Volume 36, Number 2 (1987).

       fSS - flattened binary significant source matrix. */

    // Shortcut for the case of length 1 strings
    if (n == 1) {
        return 1;
    }

    // i - starting position of search for extending substring in string we 
    //     have so far.
    // k - number of characters into new substring which is extending the
    //     string we have so far.
    // l - starting position of new substring which is extending the string we
    //     have so far.
    int i = 0, k = 1, l = 1;
    // c - LZ76 complexity index.
    // k_max - keeps track of longest possible extending substring which can
    //         be reproduced with existing substrings
    int c = 1, k_max = 1;
    while (true) {

        // Iterate along string to find existing character which matches
        // upcoming character
        if (fSS[i + k - 1] == fSS[l + k - 1]) {

            // If there's a match, check that the next character in the
            // existing substring matches the next upcoming character. This is
            // achieved by incrementing k
            k = k + 1;
 
            // If we've reached the end of fSS while searching for the
            // extending substring in what we have of the total string so far,
            // then terminate search.
            if (l + k >= n) {
                c = c + 1;
                break;
            }

        }

        else {

            if (k > k_max) {
               k_max = k;
            }

            // Start searching for extending substring using the next element
            // in the string we have so far as starting point
            i = i + 1;

            // If we've reached the beginning of the extending substring while
            // searching for the extending substring in what we have of the
            // total string so far, then the extending substring doesn't
            // already exist
            if (i == l) {

                c = c + 1;
                l = l + k_max;

                // We've reached the end of the string
                if (l >= n) {
                    break;
                }

                else {
                    // Start at the beginning of the string we have so far
                    i = 0;
                    k = 1;
                    k_max = 1;
                }
            }

            // If we haven't reached extending substring, we just want to reset
            // k = 1 so we start at beginning of extending substring again
            else {
                k = 1;
            }
        }
    }

    return c;

}


void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  bool *fSS;
  double *result;
  size_t mrows,ncols;
  
  /* Check for proper number of arguments. */
  if(nrhs!=1) {
    mexErrMsgIdAndTxt( "MATLAB:lz76:invalidNumInputs",
            "One input required.");
  } else if(nlhs>1) {
    mexErrMsgIdAndTxt( "MATLAB:lz76:maxlhs",
            "Too many output arguments.");
  }
  
  /* The input must be a logical (boolean) 1D array.*/
  mrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);
  if (!mxIsLogical(prhs[0])) {
    mexErrMsgIdAndTxt( "MATLAB:lz76:inputNotLogicalArray",
            "Input must be a 1D array of logicals.");
  }
  if ((mrows!=1) && (ncols!=1)) {
    mexErrMsgIdAndTxt( "MATLAB:lz76:inputNot1DArray",
            "Input must be a 1D array of logicals.");
  }
  
  /* Create matrix for the return argument. */
  plhs[0] = mxCreateDoubleMatrix((mwSize) 1, (mwSize) 1, mxREAL);
  
  /* Assign pointers to each input and output. */
  fSS = mxGetLogicals(prhs[0]);
  result = mxGetPr(plhs[0]);

  size_t numel = (mrows > ncols) ? mrows : ncols;
  
  result[0] = LZc(fSS, numel);

}


