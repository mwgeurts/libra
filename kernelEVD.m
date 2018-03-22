function [P,T,L,r,centerX,cX] = kernelEVD(X,varargin)

%KERNELEVD performs an eigenvector-eigenvalue decomposition of a matrix with more
% columns than rows than columns, based on its kernel eigenvalue
% decomposition.
%
% Required input: 
%     X : data matrix of size n by p where n < p (else classSVD is invoked)
%
% Optional input arguments:
%        all : if equal to 1, all eigenvalues and eigenvectors are included
%              in the output. If equal to 0, only the first r eigenvalues and 
%              eigenvectors are included (default). 
%   centered : if equal to 1, X is already centered. If equal to 0, X is
%              first mean-centered (default).
%
% Example:
%  [P,T,L,r,centerX,cX] = kernelEVD(X)
%  [P,T,L,r,centerX,cX] = kernelEVD(X,0,1)
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
% Written by Sabine Verboven, Mia Hubert
% Last Update: 17/06/2003, 25/04/2016

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
if n > p 
    [P,T,L,r,centerX,cX]=classSVD(X,all,centered);
else
    if ~centered
        cX=mean(X);
        centerX= bsxfun(@minus,X,cX);
    else
        cX=zeros(1,p);
        centerX=X;
    end
    [P,L]=eig(centerX*centerX'/(n-1));
    [L,I]=greatsort(diag(L));
    P=P(:,I);
    tol=n*eps(L(1));
    r=sum(L>tol);
    
    if ~all
        L=L(1:r);
        loadings=(centerX/sqrt(n-1))'*P(:,1:r)*diag(1./sqrt(L)); 
        %normalizing loadings by dividing by the sqrt(eigenvalues)
    else
        loadings=(centerX/sqrt(n-1))'*P*diag(1./sqrt(L));   
    end
    T=centerX*loadings;
    P=loadings;
end