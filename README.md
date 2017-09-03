# bksvd -- Block Krylov Singular Value Decomposition

Simple MATLAB code for iterative computing an SVD or PCA via the randomized block Krylov method analyzed in 
[Randomized Block Krylov Methods for Stronger and Faster Approximate Singular Value Decomposition](https://papers.nips.cc/paper/5735-randomized-block-krylov-methods-for-stronger-and-faster-approximate-singular-value-decomposition)

## Installataion

Download `bksvd.m`, [add to MATLAB path](https://www.mathworks.com/help/matlab/ref/addpath.html), or include directly in project directory.

## Documentation

`bksvd` can be used as a drop-in replacement for MATLAB's `svds` function for computing the Singular Value Decomposition (SVD). It can also be used for Principal Component Analysis (PCA) which performs an SVD on the after mean-centering the data matrix.

### Usage
**Input:**
`bksvd(A, k, iter, bsize, center)`

- `A` : matrix to decompose
- `k` : number of singular vectors to compute, default = 6
- `iter` : number of iterations, default = 3, increase for higher accuracy
- `bsize` : block size, must be >= k, default = k
- `center` : set to `true` if A's rows should be mean-centered before the singular value decomposition (e.g. when performing PCA), default = `false`

**Output:**
k singular vector/value pairs

- `U` : an orthogonal matrix whose columns are approximate top k left singular vectors for `A`
- `S` : a diagonal matrix whose entries are `A`'s approximate top k singular values
- `V` : an orthogonal matrix whose columns are approximate top k right singular vectors for `A`

`U*S*V'` is a near optimal low-rank approximation for `A`

### Example

**Standard Singular Values Decomposition**:

```
% generate test matrix
s = 1.5.^[-40:.5:40];
A = randn(10000,161)*diag(s)*randn(161,161);

% compute SVD
[U,S,V] = bksvd(A,10);
```

`bksvd` is typically as accurate as `svds` and often faster:
```
tic; [U,S,V] = svds(A,30); toc;
Elapsed time is 1.380471 seconds.
norm(A- U*S*V')/norm(A)
0.0018
```

```
tic; [U,S,V] = bksvd(A,30); toc;
Elapsed time is 0.062798 seconds.
norm(A- U*S*V')/norm(A)
0.0018
```


**Principal Component Analysis**:

For PCA, `A`'s rows (data points) should be mean-centered before computing the SVD. If the `center` flag is set to `true`, `bksvd` can do this implicitly, without densifying `A`:
```
[U,S,V] = bksvd(A,10,4,10,true);
```
Here `V` contains loading vector for the top 10 principal components. `U*S` can be taken as a dimensionality for the data to 10 components.

### Parameter Tuning

For higher accuracy (at the cost of slower runtime), the number of iterations can be increased from the default of `3`, although for many matrices this is unecessary. Increasing the block size, `bsize`, so that it is > k also increases accuracy. For matrices with quickly decaying singular values, increasing block size can be more effective than increasing iterations. For details, see the [NIPS paper.](https://papers.nips.cc/paper/5735-randomized-block-krylov-methods-for-stronger-and-faster-approximate-singular-value-decomposition)

## Other Implementations

An implementation of `bksvd` is available through [MLPACK](http://mlpack.org/docs/mlpack-git/doxygen.php?doc=classmlpack_1_1svd_1_1RandomizedBlockKrylovSVD.html). We plan to upload a Python implementation soon.
