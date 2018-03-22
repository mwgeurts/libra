function result=twopoints(data,ndirect,seed,centered)

%TWOPOINTS calculates ndirect directions through two randomly chosen data points from data.
% If ndirect is larger than the number of all possible directions, then all
% these combinations are considered.
%
% Required input arguments: 
%     data    : Data matrix
%     ndirect : Number of directions through two random data points that
%               needs to be constructed
%
% Optional input arguments:
%      seed : To define the state of the generator (default=0)
%             (0 sets the generator to its default initial state)
%  centered : If one, the first point is the center and the second point is 
%             the sum of any two random data points (can be twice the same point).
%             (default=0)   
%
%I/O:
%    result=twopoints(x,250,0);
%
% Output arguments:
%   result : matrix containing the ndirect directions (each row is a
%            direction)
% 
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust
% 
% Last modified: 25/04/2016

 if nargin==2
     seed=0;
     centered=0;
 end

 if nargin==3
    centered=0; 
 end
 
[n,p]=size(data);
if centered==1
    nrich1=(n*(n-1)/2) + n;
else
    nrich1=n*(n-1)/2;
end
nrich2=n*(n-1)/2;

ndirect=min(ndirect,nrich1);
B=zeros(ndirect,p);
dirs=zeros(nrich2,2);

perm=[1 1];
k=1;

for ndir=1:nrich2
    k1=2;
    perm(k1)=perm(k1)+1;
    while ~(k1==1 || perm(k1) <=(n-(k+1-k1)))
        k1=k1-1;
        perm(k1)=perm(k1)+1;
        for j=(k1+1):k+1
            perm(j)=perm(j-1)+1;
        end
    end
    dirs(ndir,:)=perm;
end

if centered==1
	if ndirect<=n
        index1=randomset(n,ndirect,seed);
        B(1:ndirect,:)=data(index1,:);
    else
        B(1:n,:)=data(1:n,:);
        index1=randomset(nrich2,(ndirect-n),seed);
        index2=dirs(index1,:);
        B((n+1):ndirect,:)=(data(index2(:,1),:)+data(index2(:,2),:))/2;
	end
    
else
    index1=randomset(nrich2,ndirect,seed);
	index2=dirs(index1,:);
    B(1:ndirect,:)=data(index2(:,1),:)-data(index2(:,2),:);

end
result=B;
