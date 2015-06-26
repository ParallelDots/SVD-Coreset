# SVD-Coreset
Calculate a coreset of a matrix in order to find its SVD.

This module helps to calculate a coreset of given data points, which can be used to calculate its SVD/PCA. 
A coreset can be calculated upon a huge matrix, giving a set of a few effective vectors which can then be easily used to find the SVD of the data points, which otherwise would require way more resources.
For more details refer to the paper: http://people.csail.mit.edu/dannyf/big/subspace.pdf . There is a Matlab version somewhere on the MIT website and can be found by a Google Search.
