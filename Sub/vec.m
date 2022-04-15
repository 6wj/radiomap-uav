function [Xv] = vec(X);
%
% [Xv] = vec(X)
%
% Vectorizing (vec) operator.
%
% This operator stacks the columns of the matrix X into a 
% single column vector Xv.
%
% See - Brewer, J. W. &quot;Kronecker Products and Matrix Calculus in System Theory,&quot;
%       IEEE Transactions on Circuits and Systems, Vol. cas-25, no. 9, pp. 772-781.
%

[a,b] = size(X);

Xv = reshape(X,a*b,1);