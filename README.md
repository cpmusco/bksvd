# bksvd -- Block Krylov Singular Value Decomposition

Simple MATLAB code for iterative computing an SVD or PCA via the randomized block Krylov method analyzed in our NIPS 2015 paper,
[Randomized Block Krylov Methods for Stronger and Faster Approximate Singular Value Decomposition](https://papers.nips.cc/paper/5735-randomized-block-krylov-methods-for-stronger-and-faster-approximate-singular-value-decomposition)

## Installataion

Download `bksvd.m`, [add to MATLAB path](https://www.mathworks.com/help/matlab/ref/addpath.html), or include directly in project directory.

## Documentation

`bksvd.m` can be used as a drop-in replacement for MATLAB's `svds` function for computing the Singular Value Decomposition (SVD). It can also be used for Principal Component Analysis (PCA) which performs and SVD on the mean centered data matrix.

### Usage
**Input**
`bksvd(A, k, iter, bsize, center)`

- `A` : matrix to decompose
- `k` : number of singular vectors to compute, default = 6
- `iter` : number of iterations, default = 3, increase for higher accuracy
- `bsize` : block size, must be >= k, default = k
- `center` : set to `true` if A's rows should be mean-centered before the singular value decomposition (e.g. when performing principal component analysis), default = `false`

**Output**
k singular vector/value pairs

- `U` : an orthogonal matrix whose columns are approximate top k left singular vectors for `A`
- `S` : a diagonal matrix whose entries are `A`'s approximate top k singular values
- `V` : an orthogonal matrix whose columns are approximate top k right singular vectors for `A`

`U*S*V'` is a near optimal low-rank approximation for `A`

### Example

**Standard Singular Values Decomposition**:

```
% generate random test matrix
s = 1.5.^[-40:.5:40];
A = randn(10000,161)*diag(s)*randn(161,161);
% compute SVD
[U,S,V] = bksvd(A,10);
```

`bksvd` is typically as accurate as `svds` and often faster:
```
tic; [U,S,V] = svds(A,30); toc;
```
Elapsed time is 1.380471 seconds.
```
norm(A- U*S*V')/norm(A)
```
0.0018

```
tic; [U,S,V] = bksvd(A,30); toc;
```
Elapsed time is 0.062798 seconds.
```
norm(A- U*S*V')/norm(A)
```
0.0018

**Principal Component Analysis**:
For Principal Componenent Analysis, `A`'s rows should be mean centered before performing Singular Value Decomposition. If the center = `true` flag is set, `bksvd` can do this implicitly, without densifying `A`:
```
[U,S,V] = bksvd(A,10,4,10,true);
```
Here `V` contains loading vector for the top 10 principal components. `U*S` can be taken as a dimensionality reduction of `A` to 10 components.

### Parameter Tuning

## Other Implementations

An implementation of `bksvd` is available through [MLPACK](http://mlpack.org/docs/mlpack-git/doxygen.php?doc=classmlpack_1_1svd_1_1RandomizedBlockKrylovSVD.html).

A Python implementation is coming soo.
