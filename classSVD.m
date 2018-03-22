function [P,T,L,r,centerX,cX] = classSVD(X,varargin)

%CLASSSVD performs an eigenvector-eigenvalue decomposition of a matrix with more
% rows than columns, based on its singular value decomposition.
%
% Required input: 
%          X : data matrix of size n by p where n > p
%
% Optional input arguments:
%        all : if equal to 1, all eigenvalues and eigenvectors are included
%              in the output. If equal to 0, only the first r eigenvalues and 
%              eigenvectors are included (default). 
%   centered : if equal to 1, X is already centered. If equal to 0, X is
%              first mean-centered (default).
%
% Example:
%  [P,T,L,r,centerX,cX] = classSVD(X)
%  [P,T,L,r,centerX,cX] = classSVD(X,0,1)
%
% The output consist of the following elements:
%          P  : loading matrix containing the eigenvectors of
%                  (X-mean(X))'(X-mean(X)) if centered = 0
%                   X'X                    if centered = 1
%                 The first column of P correspond with the largest
%                 eigenvalue, the second column with the second largest
%                 eigenvalue and so on.
%          T  : scores matrix 
%          L  : vector of eigenvalues in descending order
%          r  : number of non-zero eigenvalues
%    centerX  : mean-centered X
%         cX  : mean of X if centered = 0. 
%               Otherwise a vector of length p with zeros.
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust
%
% Written by Sabine Verboven, Mia Hubert, Wannes Van den Bossche
% Last Update: 25/04/2016

if nargin>2
    all=varargin{1};
    centered=varargin{2};
elseif nargin>1
    all=varargin{1};
    centered=0;
else
    all=0;
    centered=0;
end
    
[n,p]=size(X); 

if n==1
    error('The sample size is 1. No SVD can be performed.')
end

if ~centered
    cX=mean(X);
    centerX= bsxfun(@minus,X,cX);
else
    cX=zeros(1,p);
    centerX=X;
end
[~,S,loadings]=svd(centerX./sqrt(n-1),0); 
eigenvalues=diag(S).^2;
tol = max([n p])*eps(eigenvalues(1));
r=sum(eigenvalues>tol);
if ~all
    L=eigenvalues(1:r);
    P=loadings(:,1:r);
else
    L=eigenvalues;
    P=loadings;
end
T=centerX*P; 

