function [U,S,V] = sisvd(A, k, iter, bsize, center)
%--------------------------------------------------------------------------
% Simple randomized Simultaneous Iteration for truncated SVD
% Computes approximate top singular vectors and corresponding values
% Described in Rokhlin, Szlam, Tygert, 2009 (https://arxiv.org/abs/0809.2274)
%
% usage : 
%
%  input:
%  * A : matrix to decompose
%  * k : number of singular vectors to compute, default = 6
%  * iter : number of iterations, default = 3
%  * bsize : block size, must be >= k, default = k
%  * center : set to true if A's rows should be mean centered before the
%  singular value decomposition (e.g. when performing principal component 
%  analysis), default = false
%
%
%  output:
%  k singular vector/value pairs. 
%  * U : a matrix whose columns are approximate top left singular vectors for A
%  * S : a diagonal matrix whose entries are A's approximate top singular values
%  * V : a matrix whose columns are approximate top right singular vectors for A
%
%  U*S*V' is a near optimal rank-k approximation for A
%--------------------------------------------------------------------------

% Check input arguments and set defaults.
if nargin > 5
    error('sisvd:TooManyInputs','requires at most 5 input arguments');
end
if nargin < 1
    error('sisvd:TooFewInputs','requires at least 1 input arguments');
end
if nargin < 2
    k = 6;
end
k = min(k,min(size(A)));

if nargin < 3
    iter = 3;
end
if nargin < 4
    bsize = k;
end
if nargin < 5
    center = false;
end
if(k < 1 || iter < 1 || bsize < k)
    error('bksvd:BadInput','one or more inputs outside required range');
end

% Calculate row mean if rows should be centered.
u = zeros(1,size(A,2));
if(center)
    u = mean(A);
end
l = ones(size(A,1),1);


% We want to iterate on the smaller dimension of A.
[n, ind] = min(size(A));
tpose = false;
if(ind == 1) 
    tpose = true;
    l = u'; u = ones(1,size(A,1));
    A = A';
end

% Random block initialization.
block = randn(size(A,2),bsize);
[block,R] = qr(block,0);

% Preallocate space for temporary products.
T = zeros(size(A,2),bsize);

% Run power iteration, orthogonalizing at each step using economy size QR.
for i=1:iter
    T = A*block - l*(u*block);
    block = A'*T - u'*(l'*T);
    [block,R] = qr(block,0);
end

% Rayleigh-Ritz postprocessing with economy size dense SVD.
T = A*block - l*(u*block);

[Ut,St,Vt] = svd(T,0);
S = St(1:k,1:k);
if(~tpose)
    U = Ut(:,1:k);
    V = block*Vt(:,1:k);
else
    V = Ut(:,1:k);
    U = block*Vt(:,1:k);
end

end