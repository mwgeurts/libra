function [rew,raw]=mcdcov(x,varargin)

%MCDCOV computes the MCD estimator of a multivariate data set.  This
% estimator is given by the subset of h observations with smallest covariance
% determinant.  The MCD location estimate is then the mean of those h points,
% and the MCD scatter estimate is their covariance matrix.  The default value
% of h is roughly 0.75n (where n is the total number of observations), but the
% user may choose each value between n/2 and n. Based on the raw estimates,
% weights are assigned to the observations such that outliers get zero weight.
% The reweighted MCD estimator is then given by the mean and covariance matrix
% of the cases with non-zero weight. To compute the MCD estimator,
% the FASTMCD algorithm is used.
%
% The MCD method is intended for continuous variables, and assumes that
% the number of observations n is at least 5 times the number of variables p.
% If p is too large relative to n, it would be better to first reduce
% p by variable selection or robust principal components (see the functions
% robpca.m and rapca.m).
%
% The MCD method was introduced in:
%
%   Rousseeuw, P.J. (1984), "Least Median of Squares Regression,"
%   Journal of the American Statistical Association, Vol. 79, pp. 871-881.
%
% The MCD is a robust method in the sense that the estimates are not unduly
% influenced by outliers in the data, even if there are many outliers.
% Due to the MCD's robustness, we can detect outliers by their large
% robust distances. The latter are defined like the usual Mahalanobis
% distance, but based on the MCD location estimate and scatter matrix
% (instead of the nonrobust sample mean and covariance matrix).
%
% The FASTMCD algorithm uses several time-saving techniques which
% make it available as a routine tool to analyze data sets with large n,
% and to detect deviating substructures in them. A full description of the
% algorithm can be found in:
%
%   Rousseeuw, P.J. and Van Driessen, K. (1999), "A Fast Algorithm for the
%   Minimum Covariance Determinant Estimator," Technometrics, 41, pp. 212-223.
%
% An important feature of the FASTMCD algorithm is that it allows for exact
% fit situations, i.e. when more than h observations lie on a (hyper)plane.
% Then the program still yields the MCD location and scatter matrix, the latter
% being singular (as it should be), as well as the equation of the hyperplane.
%
%
% Required input argument:
%    x : a vector or matrix whose columns represent variables, and rows represent observations.
%        Missing values (NaN's) and infinite values (Inf's) are allowed, since observations (rows)
%        with missing or infinite values will automatically be excluded from the computations.
%
% Optional input arguments:
%       cor : If non-zero, the robust correlation matrix will be
%             returned. The default value is 0.
%         h : The quantile of observations whose covariance determinant will
%             be minimized.  Any value between n/2 and n may be specified.
%             The default value is 0.75*n.
%     alpha : (1-alpha) measures the fraction of outliers the algorithm should
%             resist. Any value between 0.5 and 1 may be specified. (default = 0.75)
%    ntrial : The number of random trial subsamples that are drawn for
%             large datasets. The default is 500.
%     plots : If equal to one, a menu is shown which allows to draw several plots,
%             such as a distance-distance plot. (default)
%             If 'plots' is equal to zero, all plots are suppressed.
%             See also makeplot.m
%   classic : If equal to one, the classical mean and covariance matrix are computed as well.
%             (default = 0)
%   center  : If equal to one the dataset is considered to be centered, i.e. the
%             MCD location of the data is considered to be the origin.
%
% Input arguments for advanced users:
%     Hsets : Instead of random trial h-subsets (default, Hsets = []), Hsets makes it possible to give certain
%             h-subsets as input. Hsets is a matrix that contains the indices of the observations of one
%             h-subset as a row.
%    factor : If not equal to 0 (default), the consistency factor is adapted. Only useful in case of the
%             kmax approach.
%
% I/O: result=mcdcov(x,'alpha',0.75,'h',h,'ntrial',500)
%  If only one output argument is listed, only the final result ('result')
%  is returned.
%  The user should only give the input arguments that have to change their default value.
%  The name of the input arguments needs to be followed by their value.
%  The order of the input arguments is of no importance.
%
% Examples: [rew,raw]=mcdcov(x);
%           result=mcdcov(x,'h',20,'plots',0);
%           [rew,raw]=mcdcov(x,'alpha',0.8,'cor',0)
%
% The output structure 'raw' contains intermediate results, with the following
% fields :
%
%     raw.center : The raw MCD location of the data.
%        raw.cov : The raw MCD covariance matrix (multiplied by a consistency factor).
%        raw.cor : The raw MCD correlation matrix, if input argument 'cor' was non-zero.
%         raw.wt : Weights based on the estimated raw covariance matrix 'raw.cov' and
%                  the estimated raw location 'raw.center' of the data. These weights determine
%                  which observations are used to compute the final MCD estimates.
%  raw.objective : The determinant of the raw MCD covariance matrix.
%
% The output structure 'rew' contains the final results, namely:
%
%       rew.center : The robust location of the data, obtained after reweighting, if
%                    the raw MCD is not singular.  Otherwise the raw MCD center is
%                    given here.
%          rew.cov : The robust covariance matrix, obtained after reweighting, if the raw MCD
%                    is not singular.  Otherwise the raw MCD covariance matrix is given here.
%          rew.cor : The robust correlation matrix, obtained after reweighting, if
%                    options.cor was non-zero.
%            rew.h : The number of observations that have determined the MCD estimator,
%                    i.e. the value of h.
%    rew. Hsubsets : A structure that contains Hopt and Hfreq:
%                    Hopt  : The subset of h points whose covariance matrix has minimal determinant,
%                            ordered following increasing robust distances.
%                    Hfreq : The subset of h points which are the most frequently selected during the whole
%                            algorithm.
%        rew.alpha : (1-alpha) measures the fraction of outliers the algorithm should
%                    resist.
%           rew.md : The distance of each observation from the classical
%                    center of the data, relative to the classical shape
%                    of the data. Often, outlying points fail to have a
%                    large Mahalanobis distance because of the masking
%                    effect.
%           rew.rd : The distance of each observation to the final,
%                    reweighted MCD center of the data, relative to the
%                    reweighted MCD scatter of the data.  These distances allow
%                    us to easily identify the outliers. If the reweighted MCD
%                    is singular, raw.rd is given here.
%       rew.cutoff : Cutoff values for the robust and mahalanobis distances
%         rew.flag : Flags based on the reweighted covariance matrix and the
%                    reweighted location of the data.  These flags determine which
%                    observations can be considered as outliers. If the reweighted
%                    MCD is singular, raw.wt is given here.
%       rew.method : A character string containing information about the method and
%                    about singular subsamples (if any).
%        rew.plane : In case of an exact fit, rew.plane is a structure that contains
%                    eigvct and eigval:
%                    eigvct: contains the eigenvectors that define the plane. I.e.
%                            the eigenvectors belonging to none-zero
%                            eigenvalues.
%                    eigval: contains the corresponding eigenvalues of
%                            eigvct.
%      rew.classic : If the input argument 'classic' is equal to one, this structure
%                    contains results of the classical analysis: center (sample mean),
%                    cov (sample covariance matrix), md (Mahalanobis distances), class ('COV').
%        rew.class : 'MCDCOV'
%            rew.X : If x is bivariate, same as the x in the call to mcdcov,
%                    without rows containing missing or infinite values.
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at:
%              http://wis.kuleuven.be/stat/robust
%
% Written by Katrien Van Driessen and Bjorn Rombouts
% Revisions by Sanne Engelen, Sabine Verboven, Wannes Van den Bossche
% Last Update: 09/04/2004, 01/08/2007, 25/04/2016

% The FASTMCD algorithm works as follows:
%
%       The dataset contains n cases and p variables.
%       When n < 2*nmini (see below), the algorithm analyzes the dataset as a whole.
%       When n >= 2*nmini (see below), the algorithm uses several subdatasets.
%
%       When the dataset is analyzed as a whole, a trial subsample of p+1 cases
%       is taken, of which the mean and covariance matrix are calculated.
%       The h cases with smallest relative distances are used to calculate
%       the next mean and covariance matrix, and this cycle is repeated csteps1
%       times. For small n we consider all subsets of p+1 out of n, otherwise
%       the algorithm draws 500 random subsets by default.
%       Afterwards, the 10 best solutions (means and corresponding covariance
%       matrices) are used as starting values for the final iterations.
%       These iterations stop when two subsequent determinants become equal.
%       (At most csteps3 iteration steps are taken.) The solution with smallest
%       determinant is retained.
%
%       When the dataset contains more than 2*nmini cases, the algorithm does part
%       of the calculations on (at most) maxgroup nonoverlapping subdatasets, of
%       (roughly) maxobs cases.
%
%       Stage 1: For each trial subsample in each subdataset, csteps1 (see below) iterations are
%       carried out in that subdataset. For each subdataset, the 10 best solutions are
%       stored.
%
%       Stage 2 considers the union of the subdatasets, called the merged set.
%       (If n is large, the merged set is a proper subset of the entire dataset.)
%       In this merged set, each of the 'best solutions' of stage 1 are used as starting
%       values for csteps2 (sse below) iterations. Also here, the 10 best solutions are stored.
%
%       Stage 3 depends on n, the total number of cases in the dataset.
%       If n <= 5000, all 10 preliminary solutions are iterated.
%       If n > 5000, only the best preliminary solution is iterated.
%       The number of iterations decreases to 1 according to n*p (If n*p <= 100,000 we
%       iterate csteps3 (sse below) times, whereas for n*p > 1,000,000 we take only one iteration step).
%

if rem(nargin-1,2)~=0
    error('The number of input arguments should be odd!');
end
% Assigning some input parameters
data = x;
raw.cor = [];
rew.cor = [];
rew.plane = [];
% The maximum value for n (= number of observations) is:
nmax=50000;
% The maximum value for p (= number of variables) is:
pmax=50;
% To change the number of subdatasets and their size, the values of
% maxgroup and nmini can be changed.
maxgroup=5;
nmini=300;
% The number of iteration steps in stages 1,2 and 3 can be changed
% by adapting the parameters csteps1, csteps2, and csteps3.
csteps1=2;
csteps2=2;
csteps3=100;
% dtrial : number of subsamples if not all (p+1)-subsets will be considered.
dtrial=500;

if size(data,1)==1
    data=data';
end

% Observations with missing or infinite values are ommitted.
ok=all(isfinite(data),2);
data=data(ok,:);
xx=data;
[n,p]=size(data);
% Some checks are now performed.
if n==0
    error('All observations have missing or infinite values.')
end
if n > nmax
    error(['The program allows for at most ' int2str(nmax) ' observations.'])
end
if p > pmax
    error(['The program allows for at most ' int2str(pmax) ' variables.'])
end
if n < p
    error('Need at least (number of variables) observations.')
end

%internal variables
hmin=quanf(0.5,n,p);
%Assiging default values
h=quanf(0.75,n,p);
default=struct('alpha',0.75,'h',h,'plots',1,'ntrial',dtrial,'cor',0,'seed',0,'classic',0,'center',0,'Hsets',[],'factor',0);
list=fieldnames(default);
options=default;
IN=length(list);
i=1;
counter=1;

%Reading optional inputarguments
if nargin>2
    %
    % placing inputfields in array of strings
    %
    chklist=cell(1,floor((nargin-1)/2));
    for j=1:nargin-1
        if rem(j,2)~=0
            varargin{j}=strtrim(varargin{j}); %remove white space
            chklist{i}=varargin{j};
            i=i+1;
        end
    end

    %
    % Checking which default parameters have to be changed
    % and keep them in the structure 'options'.
    %
    while counter<=IN
        index=find(strcmp(list(counter,:),chklist));
        if ~isempty(index) % in case of similarity
            for j=1:nargin-2 % searching the index of the accompanying field
                if rem(j,2)~=0 % fieldnames are placed on odd index
                    if strcmp(chklist{index},varargin{j})
                        I=j;  %index of located similarity
                    end
                end
            end
            options.(chklist{index})=varargin{I+1};
%             options=setfield(options,chklist{index},varargin{I+1});
        end
        counter=counter+1;
    end

    dummy=sum(strcmp(chklist,'h')+2*strcmp(chklist,'alpha'));
    switch dummy

        case 1% checking inputvariable h
        % hmin is the minimum number of observations whose covariance determinant
        % will be minimized.            
            if isempty(options.Hsets)
                if options.h < hmin
                    disp(['Warning: The MCD must cover at least ' int2str(hmin) ' observations.'])
                    disp(['The value of h is set equal to ' int2str(hmin)])
                    options.h = hmin;
                elseif options.h > n
                    error('h is greater than the number of non-missings and non-infinites.')
                elseif options.h < p
                    error(['h should be larger than the dimension ' int2str(p) '.'])
                end
            end
            options.alpha=options.h/n;            
        case 2            
            if options.alpha < 0.5
                options.alpha=0.5;
                mess=sprintf(['Attention (mcdcov.m): Alpha should be larger than 0.5. \n',...
                    'It is set to 0.5.']);
                disp(mess)
            end
            if options.alpha > 1
                options.alpha=0.75;
                mess=sprintf(['Attention (mcdcov.m): Alpha should be smaller than 1.\n',...
                    'It is set to 0.75.']);
                disp(mess)
            end
            options.h=quanf(options.alpha,n,p);        
        case 3
            error('Both input arguments alpha and h are provided. Only one is required.')
    end
end

centered=options.center; %data is considered centered if one
h=options.h;  %number of regular datapoints on which estimates are based. h=[alpha*n]
plots=options.plots; %relevant plots are plotted
alfa=options.alpha; %percentage of regular observations
ntrial=options.ntrial; %number of subsets to be taken in the first step
cor=options.cor; %correlation matrix
seed=options.seed; %seed of the random generator
cutoff.rd=sqrt(chi2inv(0.975,p)); %cutoff value for the robust distance
cutoff.md=cutoff.rd; %cutoff value for the mahalanobis distance
Hsets = options.Hsets;
if ~isempty(Hsets)
    Hsets_ind = 1;
else
    Hsets_ind = 0;
end
factor = options.factor;
if factor == 0
    factor_ind = 0;
else
    factor_ind = 1;
end

% Some initializations.
rew.flag=NaN(1,length(ok));
raw.wt=NaN(1,length(ok));
raw.rd=NaN(1,length(ok));
rew.rd=NaN(1,length(ok));
rew.mahalanobis=NaN(1,length(ok));
rew.method=sprintf('\nMinimum Covariance Determinant Estimator.');
correl=NaN;

% weights : weights of the observations that are not excluded from the computations.
%           These are the observations that don't contain missing or infinite values.
% bestobj : best objective value found.
weights=zeros(1,n);
bestobj=inf;

% The classical estimates are computed 
% [Loadings,Scores,Eigs,#nonzero eigs,centeredData,meanData]
[Pcl,~,Lcl,rcl,~,cXcl] = classSVD(data,0,centered);
clasmean=cXcl;
clascov=Pcl*diag(Lcl)*Pcl';

if p < 5
    eps=1e-12;
elseif p <= 8
    eps=1e-14;
else
    eps=1e-16;
end

%h-th order zero statistic
if ~centered
    med=median(data);
else
    med=zeros(1,p);
end
mad=sort(abs(bsxfun(@minus,data,med)));
mad=mad(h,:);
iis=find( (mad < eps)); %all vars with h-th order zero statistic
if ~isempty(iis) 
    ii=iis(1); %smallest variable
    % The h-th order statistic is zero for the ii-th variable. The array plane contains
    % all the observations which have the same value for the ii-th variable.
    plane=find(abs(data(:,ii)-med(ii)) < eps)';
    weights(plane)=1;
    if p==1
        if ~centered 
            meanplane=mean(data(plane,:));
        else
            meanplane=zeros(1,p);
        end
        rew.flag=weights;
        raw.wt=weights;
        [raw.center,rew.center]=deal(meanplane);
        [raw.cov,rew.cov,raw.objective]=deal(0);
        if plots
            rew.method=sprintf('\nUnivariate location and scale estimation.');
            rew.method=char(rew.method,sprintf('%g of the %g observations are identical.',length(plane),n));
        end
    else
        [P,~,L,~,~,meanplane] = classSVD(data(plane,:),0,centered);
        rew.plane.eigvct=P;
        rew.plane.eigval=L;
        covplane=P*diag(L)*P';
        [raw.center,raw.cov,rew.center,rew.cov,raw.objective,raw.wt,rew.flag,...
            rew.method]=displ(3,length(plane),weights,n,p,meanplane,covplane,rew.method,rew.plane,h,iis);
    end
    rew.Hsubsets.Hopt = plane;
    rew.Hsubsets.Hfreq = plane;
    %classical analysis?
    if options.classic==1
        classic.cov=clascov;
        classic.center=clasmean;
        classic.class='COV';
    else
        classic=0;
    end
    %assigning the output
    rewo=rew;rawo=raw;
    rew=struct('center',{rewo.center},'cov',{rewo.cov},'cor',{cor},'h',{h},'Hsubsets',{rewo.Hsubsets},'alpha',{alfa},...
        'flag',{rewo.flag},'plane', {rewo.plane},'method',{rewo.method},'class',{'MCDCOV'},'classic',{classic},'X',{xx});
    raw=struct('center',{rawo.center},'cov',{rawo.cov},'cor',{rawo.cor},'objective',{rawo.objective},...
        'wt',{rawo.wt},'class',{'MCDCOV'},'classic',{classic},'X',{x});
    if size(data,2)~=2
        rew=rmfield(rew,'X');
        raw=rmfield(raw,'X');
    end
    return
end

%exact fit situation
if rcl < p
    % all observations lie on a hyperplane.
    rew.plane.eigvct=Pcl;%./mad';
    rew.plane.eigval=Lcl;
    weights(1:n)=1;
    if cor
        correl=clascov./(sqrt(diag(clascov))*sqrt(diag(clascov))');
        [rew.cor,raw.cor]=deal(correl);
    end
    [raw.center,raw.cov,rew.center,rew.cov,raw.objective,raw.wt,rew.flag,...
        rew.method]=displ(1,n,weights,n,p,clasmean,clascov,rew.method,rew.plane);

    %classical analysis?
    if options.classic==1
        classic.cov=clascov;
        classic.center=clasmean;
        classic.class='COV';
    else
        classic=0;
    end
    rew.Hsubsets.Hopt=1:n;
    rew.Hsubsets.Hfreq=1:n;
    %assigning the output
    rewo=rew;rawo=raw;
    rew=struct('center',{rewo.center},'cov',{rewo.cov},'cor',{rewo.cor},'h',{h},'Hsubsets',{rewo.Hsubsets},'alpha',{alfa},...
        'rd',{rewo.rd},'cutoff',{cutoff},'flag',{rewo.flag},'plane',{rewo.plane},'method',{rewo.method},...
        'class',{'MCDCOV'},'classic',{classic},'X',{xx});
    raw=struct('center',{rawo.center},'cov',{rawo.cov},'cor',{rawo.cor},'objective',{rawo.objective},...
        'cutoff',{cutoff},'wt',{rawo.wt}, 'class',{'MCDCOV'},'classic',{classic},'X',{x});
    if size(data,2)~=2
        rew=rmfield(rew,'X');
        raw=rmfield(raw,'X');
    end
    return
end

% The standardization of the data will now be performed.
data=bsxfun(@rdivide,bsxfun(@minus,data,med),mad);
% The standardized classical estimates are now computed.
clmean=(clasmean-med)./mad;%mean(data);
clcov=bsxfun(@rdivide,bsxfun(@rdivide,clascov,mad),mad');%cov(data);

% The univariate non-classical case is now handled.
if p==1 && h~=n
    [rew.center,rewsca,weights,raw.center,raw.cov,raw.rd,Hopt]=unimcd(data,h,centered);
    rew.Hsubsets.Hopt = Hopt';
    rew.Hsubsets.Hfreq = Hopt';
    [raw.cov,raw.center]=trafo(raw.cov,raw.center,med,mad);
    raw.objective=raw.cov;
    raw.cutoff=cutoff.rd;
    raw.wt=weights;
    rew.cov=rewsca^2;
    mah=(data-rew.center).^2/rew.cov;
    rew.rd=sqrt(mah');
    rew.flag=(rew.rd<=cutoff.rd);
    rew.cutoff=cutoff.rd;
    [rew.cov,rew.center]=trafo(rew.cov,rew.center,med,mad);
    rew.mahalanobis=abs(data'-clmean)/sqrt(clcov);
    %classical analysis
    if options.classic==1
        classic.cov=clascov;
        classic.center=clasmean;
        classic.md=rew.mahalanobis;
        classic.class='COV';
    else
        classic=0;
    end
    %assigning the output
    rewo=rew;rawo=raw;
    rew=struct('center',{rewo.center},'cov',{rewo.cov},'cor',{rewo.cor},'h',{h},'Hsubsets',{rewo.Hsubsets},...
        'alpha',{alfa},'rd',{rewo.rd},'cutoff',{cutoff},'flag',{rewo.flag}, 'plane',{[]},'method',{rewo.method},...
        'class',{'MCDCOV'},'md',{rewo.mahalanobis},'classic',{classic},'X',{xx});
    raw=struct('center',{rawo.center},'cov',{rawo.cov},'cor',{rawo.cor},'objective',{rawo.objective},...
        'rd',{rawo.rd},'cutoff',{cutoff},'wt',{rawo.wt}, 'class',{'MCDCOV'},'classic',{classic},'X',{x});
    if size(data,2)~=2
        rew=rmfield(rew,'X');
        raw=rmfield(raw,'X');
    end
    try
        if plots && options.classic
            makeplot(rew,'classic',1)
        elseif plots
            makeplot(rew)
        end
    catch %output must be given even if plots are interrupted
        %> delete(gcf) to get rid of the menu
        return
    end
    return
end

% The classical case is now handled.
if h==n
    if plots
        msg=sprintf('The MCD estimates based on %g observations are equal to the classical estimates.\n',h);
        rew.method=char(rew.method,msg);
    end
    raw.center=clmean;
    raw.cov=clcov;
    raw.objective=det(clcov);
    mah=mahalanobis(data,clmean,'cov',clcov);
    rew.mahalanobis=sqrt(mah);
    raw.rd=rew.mahalanobis;
    weights= mah <= cutoff.rd^2;
    raw.wt=weights;
	[rew.center,rew.cov]=weightmecov(data,weights,centered);
    if cor
        raw.cor=raw.cov./(sqrt(diag(raw.cov))*sqrt(diag(raw.cov))');
        rew.cor=rew.cov./(sqrt(diag(rew.cov))*sqrt(diag(rew.cov))');
    else
        raw.cor=0;
        rew.cor=0;
    end
    [Pfull,~,r] = covsvd(rew.cov,1);
    Pzero=Pfull(:,r+1:p);
    if r<p    
        [obsinplaneInd,count,~]=obsinplane(data,rew.center,Pzero,centered);    
        if ~centered
            covar = cov(data(obsinplaneInd,:));
            meanvct=mean(data(obsinplaneInd,:));   
        else
            covar=(data(obsinplaneInd,:)'*data(obsinplaneInd,:))/(length(obsinplaneInd)-1);
            meanvct=zeros(1,p);
        end
        [covar,center]=trafo(covar,meanvct,med,mad);          
        [rew.plane.eigvct,rew.plane.eigval,~] = covsvd(covar,0);      
        if cor
            correl=covar./(sqrt(diag(covar))*sqrt(diag(covar))');
        end
        rew.method=displrw(count,n,p,center,covar,rew.method,rew.plane,cor,correl);
        rew.rd=raw.rd;
    else
        mah=mahalanobis(data,rew.center,'cov',rew.cov);
        weights = mah <= cutoff.md^2;
        rew.rd=sqrt(mah);
    end
    [raw.cov,raw.center]=trafo(raw.cov,raw.center,med,mad);
    [rew.cov,rew.center]=trafo(rew.cov,rew.center,med,mad);    
    raw.objective=raw.objective*prod(mad)^2;
    rew.flag=weights;
    %classical analysis?
    if options.classic==1
        classic.cov=clascov;
        classic.center=clasmean;
        classic.md=rew.mahalanobis;
        classic.class='COV';
    else
        classic=0;
    end
    %assigning Hsubsets:
    rew.Hsubsets.Hopt = 1:n;
    rew.Hsubsets.Hfreq = 1:n;
    %assigning the output
    rewo=rew;rawo=raw;
    rew=struct('center',{rewo.center},'cov',{rewo.cov},'cor',{rewo.cor},'h',{h},'Hsubsets',{rewo.Hsubsets},'alpha',{alfa},...
        'rd',{rewo.rd},'cutoff',{cutoff},'flag',{rewo.flag},'plane',{rewo.plane},...
        'method',{rewo.method},'class',{'MCDCOV'},'md',{rewo.mahalanobis},'classic',{classic},'X',{xx});
    raw=struct('center',{rawo.center},'cov',{rawo.cov},'cor',{rawo.cor},'objective',{rawo.objective},...
        'rd',{rawo.rd},'cutoff',{cutoff},'wt',{rawo.wt}, 'class',{'MCDCOV'},'classic',{classic},'X',{x});
    if size(data,2)~=2
        rew=rmfield(rew,'X');
        raw=rmfield(raw,'X');
    end
    try
        if plots && options.classic
            makeplot(rew,'classic',1)
        elseif plots
            makeplot(rew)
        end
    catch %output must be given even if plots are interrupted
        %> delete(gcf) to get rid of the menu
        return
    end
    return
end
percent=h/n;
teller = zeros(1,n+1);

if Hsets_ind
    csteps = csteps1;
    fine = 0;
    part = 0;
    final = 1;
    tottimes = 0;
    nsamp = size(Hsets,1);
    obsingroup = n;
else

    %  If n >= 2*nmini the dataset will be divided into subdatasets.  For n < 2*nmini the set
    %  will be treated as a whole.

    if n >= 2*nmini
        maxobs=maxgroup*nmini;
        if n >= maxobs
            ngroup=maxgroup;
            group(1:maxgroup)=nmini;
        else
            ngroup=floor(n/nmini);
            minquan=floor(n/ngroup);
            group = zeros(1,ngroup);
            group(1)=minquan;
            for s=2:ngroup
                group(s)=minquan+double(rem(n,ngroup)>=s-1);
            end
        end
        part=1;
        adjh=floor(group(1)*percent);
        nsamp=floor(ntrial/ngroup);
        minigr=sum(group);
        obsingroup=fillgroup(n,group,ngroup,seed);
        % obsingroup : i-th row contains the observations of the i-th group.
        % The last row (ngroup+1-th) contains the observations for the 2nd stage
        % of the algorithm.
    else
        [part,group,ngroup,adjh,minigr,obsingroup]=deal(0,n,1,h,n,n);
        replow=[50,22,17,15,14,zeros(1,45)];
        if n < replow(p)
            % All (p+1)-subsets will be considered.
            allSubsets=1;
            if ~centered
                perm=[1:p,p];
                nsamp=nchoosek(n,p+1);
            else
                perm=[1:p-1,p-1];
                nsamp=nchoosek(n,p);
            end
        else
            allSubsets=0;
            nsamp=ntrial;
        end
    end
    % some further initialisations.

    csteps=csteps1;
    % tottimes : the total number of iteration steps.
    % fine     : becomes 1 when the subdatasets are merged.
    % final    : becomes 1 for the final stage of the algorithm.
    [tottimes,fine,final,prevObj]=deal(0);

    if part
        % bmean1 : contains, for the first stage of the algorithm, the means of the ngroup*10
        %          best estimates.
        % bcov1  : analogous to bmean1, but now for the covariance matrices.
        % bobj1  : analogous to bmean1, but now for the objective values.
        % coeff1 : if in the k-th subdataset there are at least adjh observations that lie on
        %          a hyperplane then the coefficients of this plane will be stored in the
        %          k-th column of coeff1.
        Pzeros1=cell(ngroup,1);
        bobj1=inf(ngroup,10);
        bmean1=cell(ngroup,10);
        bP1 = cell(ngroup,10);
        bL1 = cell(ngroup,10);
        [Pzeros1{:}]=deal(NaN);
        [bmean1{:}]=deal(NaN);
        [bL1{:}]=deal(NaN);
        [bP1{:}]=deal(NaN);
    end

    % bmean : contains the means of the ten best estimates obtained in the second stage of the
    %         algorithm.
    % bcov  : analogous to bmean, but now for the covariance matrices.
    % bobj  : analogous to bmean, but now for the objective values.
    % coeff : analogous to coeff1, but now for the merged subdataset.
    % If the data is not split up, the 10 best estimates obtained after csteps1 iterations
    % will be stored in bmean, bcov and bobj.
    
    bobj=inf(1,10);
    bmean=cell(1,10);
%     Pzeros=NaN(p,1); %does not require initialization
    bP = cell(1,10);
    bL = cell(1,10);
    [bmean{:}]=deal(NaN);
    [bP{:}]=deal(NaN);
    [bL{:}] = deal(NaN);
end

seed=0;

while final~=2
    if fine || (~part && final)
        if ~Hsets_ind
            nsamp=10;
        end
        if final
            adjh=h;
            ngroup=1;
            if n*p <= 1e+5
                csteps=csteps3;
            elseif n*p <=1e+6
                csteps=10-(ceil(n*p/1e+5)-2);
            else
                csteps=1;
            end
            if n > 5000

                nsamp=1;
            end
        else
            adjh=floor(minigr*percent);
            csteps=csteps2;
        end
    end

    % found : becomes 1 if we have a singular intermediate MCD estimate.
    %         is set to 0 for every new subgroup
    found=0;

    for k=1:ngroup
        if ~fine
            found=0;
        end
        for i=1:nsamp
            tottimes=tottimes+1;
            % ns becomes 1 if we have a singular trial subsample and if there are at
            % least adjh observations in the subdataset that lie on the concerning hyperplane.
            % In that case we don't have to take C-steps. The determinant is zero which is
            % already the lowest possible value. If ns=1, no C-steps will be taken and we
            % start with the next sample. If we, for the considered subdataset, haven't
            % already found a singular MCD estimate, then the results must be first stored in
            % bmean, bcov, bobj or in bmean1, bcov1 and bobj1.  If we, however, already found
            % a singular result for that subdataset, then the results won't be stored
            % (the hyperplane we just found is probably the same as the one we found earlier.
            % We then let adj be zero. This will guarantee us that the results won't be
            % stored) and we start immediately with the next sample.
            adj=1;
            ns=0;

            % For the second and final stage of the algorithm the array sortdist(1:adjh)
            % contains the indices of the observations corresponding to the adjh observations
            % with minimal relative distances with respect to the best estimates of the
            % previous stage. An exception to this, is when the estimate of the previous
            % stage is singular.  For the second stage we then distinguish two cases :
            %
            % 1. There aren't adjh observations in the merged set that lie on the hyperplane.
            %    The observations on the hyperplane are then extended to adjh observations by
            %    adding the observations of the merged set with smallest orthogonal distances
            %    to that hyperplane.
            % 2. There are adjh or more observations in the merged set that lie on the
            %    hyperplane. We distinguish two cases. We haven't or have already found such
            %    a hyperplane. In the first case we start with a new sample.  But first, we
            %    store the results in bmean1, bcov1 and bobj1. In the second case we
            %    immediately start with a new sample.
            %
            % For the final stage we do the same as 1. above (if we had h or more observations
            % on the hyperplane we would already have found it).
            if ~Hsets_ind

                if final
                    if ~isinf(bobj(i))
                        meanvct=bmean{i};
                        P=bP{i};
                        L = bL{i};
                        if bobj(i)==0
                            [~,~,distToPlane]=obsinplane(data,meanvct,Pzeros,centered);
                            [~,sortdist]=sort(distToPlane);                            
                        else
                            [~,sortdist]=mahal2(bsxfun(@minus,data,meanvct)*P,sqrt(L),part,fine,final,k,obsingroup);
                        end
                    else
                        break
                    end
                elseif fine
                    if ~isinf(bobj1(k,i))
                        meanvct=bmean1{k,i};
                        P=bP1{k,i};
                        L=bL1{k,i};
                        if bobj1(k,i)==0
                            [~,~,distToPlane]=obsinplane(data(obsingroup{end},:),meanvct,Pzeros1{k},centered);
                            [dis,ind]=sort(distToPlane);
                            sortdist=obsingroup{end}(ind);
                            if dis(adjh) < 1e-8
                                if found==0
                                    obj=0;
                                    Pzeros=Pzeros1{k};
                                    found=1;
                                else
                                    adj=0;
                                end
                                ns=1;
                            end
                        else
                            [~,sortdist]=mahal2(bsxfun(@minus,data(obsingroup{end},:),meanvct)*P,sqrt(L),part,fine,final,k,obsingroup);
                        end
                    else
                        break;
                    end
                else %first stage: subgroups
                  	if ~centered
                        psamp = p+1;
                    else
                        psamp = p;
                    end
                    if ~part            %there are no subsets
                        if allSubsets   %all p+1 subsets out of n are considered

                            l=psamp;
                            perm(l)=perm(l)+1;
                            while ~(l==1 || perm(l) <=(n-(psamp-l)))
                                l=l-1;
                                perm(l)=perm(l)+1;
                                for j=(l+1):psamp
                                    perm(j)=perm(j-1)+1;
                                end
                            end
                            index=perm;
                        else            %take a random p+1 subset out of n
                            [index,seed]=randomset(n,psamp,seed);
                        end
                    else                %there are k=1:ngroup subsets               
                        [index,seed]=randomset(group(k),psamp,seed);
                        index=obsingroup{k}(index);
                    end
                    
                    [Pfull,~,Lfull,r,~,meanvct] = classSVD(data(index,:),1,centered);
                    Pzero=Pfull(:,r+1:p);
                    P=Pfull(:,1:r);
                    L=Lfull(1:r);
                    
                    if r==0             %all points collapse into single point
                        ns = 1;
                        rew.center=meanvct;
                        rew.cov=zeros(p,p);
                        rew.flag = zeros(1,n);
                        rew.flag(index) = 1;
                        rew.Hsubsets.Hopt=1:n;
                        rew.Hsubsets.Hfreq=1:n;
                    elseif r < p
                        % The trial subsample is singular.
                        % We distinguish two cases :
                        %
                        % 1. There are adjh or more observations in the subdataset that lie
                        %    on the hyperplane. If the data is not split up, we have adjh=h and thus
                        %    an exact fit. If the data is split up we distinguish two cases.
                        %    We haven't or have already found such a hyperplane.  In the first case
                        %    we check if there are more than h observations in the entire set
                        %    that lie on the hyperplane. If so, we have an exact fit situation.
                        %    If not, we start with a new trial subsample.  But first, the
                        %    results must be stored bmean1, bcov1 and bobj1.  In the second case
                        %    we immediately start with a new trial subsample.
                        %
                        % 2. There aren't adjh observations in the subdataset that lie on the
                        %    hyperplane. We then extend the trial subsample until it isn't singular
                        %    anymore.


                        % The smallest eigenvector belonging to a non-zero eigenvalue contains 
                        %the coefficients of the hyperplane. 
                        
                        if ~part 
                            [obsinplaneInd,count,~]=obsinplane(data,meanvct,Pzero,centered);
                        else
                            [obsinplaneInd,count,~]=obsinplane(data(obsingroup{k},:),meanvct,Pzero,centered);
                        end
                        
                        if count >= adjh
                            if ~part %full dataset, so more than h samples on plane
                                if ~centered
                                    covar = cov(data(obsinplaneInd,:));
                                    meanvct=mean(data(obsinplaneInd,:));  
                                else
                                    covar=(data(obsinplaneInd,:)'*data(obsinplaneInd,:))/(length(obsinplaneInd)-1);
                                    meanvct=zeros(1,p);
                                end                                
                                [covar,center]=trafo(covar,meanvct,med,mad);    
                                [rew.plane.eigvct,rew.plane.eigval,~] = covsvd(covar,0);                             
                                weights(obsinplaneInd)=1;
                                [raw.center,raw.cov,rew.center,rew.cov,raw.objective,...
                                    raw.wt,rew.flag,rew.method]=displ(2,count,weights,n,p,center,covar,...
                                    rew.method,rew.plane);
                                if cor
                                    correl=covar./(sqrt(diag(covar))*sqrt(diag(covar))');      
                                    [rew.cor,raw.cor]=deal(correl);                                    
                                end
                                rew.Hsubsets.Hopt=obsinplaneInd;
                                rew.Hsubsets.Hfreq=obsinplaneInd;
                                return
                            elseif found==0     %no singular sample found in this subgroup so far
                                [obsinplaneInd,count2,~]=obsinplane(data,meanvct,Pzero,centered);                 
                                if count2>=h
                                    if ~centered
                                        covar = cov(data(obsinplaneInd,:));
                                        meanvct=mean(data(obsinplaneInd,:));  
                                    else
                                        covar=(data(obsinplaneInd,:)'*data(obsinplaneInd,:))/(length(obsinplaneInd)-1);
                                        meanvct=zeros(1,p);
                                    end      
                                    [covar,center]=trafo(covar,meanvct,med,mad);    
                                    [rew.plane.eigvct,rew.plane.eigval,~] = covsvd(covar,0);   
                                    weights(obsinplaneInd)=1;
                                    [raw.center,raw.cov,rew.center,rew.cov,raw.objective,...
                                        raw.wt,rew.flag,rew.method]=displ(2,count2,weights,n,p,center,covar,...
                                        rew.method,rew.plane);
                                    if cor
                                        correl=covar./(sqrt(diag(covar))*sqrt(diag(covar))');      
                                        [rew.cor,raw.cor]=deal(correl);                                    
                                    end
                                    rew.Hsubsets.Hopt=obsinplaneInd;
                                    rew.Hsubsets.Hfreq=obsinplaneInd;
                                    return
                                end
                                obj=0;
                                Pzeros1{k}=Pzero;
                                found=1;
                                ns=1;
                            else                %a singular sample is already found in this subgroup
                                ns=1;
                                adj=0;
                            end
                        else %count < adjh
                            covmat = P*diag(L)*P';
                            while det(covmat) < exp(-50*p) %add samples until no longer singular
                                [index1,seed]=addobs(index,n,seed);
                                if ~centered
                                    [covmat,meanvct] = updatecov(data(index,:),covmat,meanvct,data(setdiff(index1,index),:),[],1);
                                else
                                    [covmat,~] = updatecov(data(index,:),covmat,meanvct,data(setdiff(index1,index),:),[],1);
                                end
                                index = index1;
                            end
                        end %count >= adjh
                    end

                    if ~ns              %no subspace plane with count >= adjh has been found: find sorted robust distances
                        if ~part
                            [~,sortdist] = mahal2(bsxfun(@minus,data,meanvct)*P,sqrt(L),part,fine,final,k,obsingroup);
                        else
                            [~,sortdist] = mahal2(bsxfun(@minus,data(obsingroup{k},:),meanvct)*P,sqrt(L),part,fine,final,k,obsingroup);
                        end
                    end
                end
            end

            if ~ns  
            %none singular sample: do C-steps
                for j=1:csteps
                    tottimes=tottimes+1;
                    if j == 1
                        if Hsets_ind
                            obs_in_set = Hsets(i,:);
                        else
                            obs_in_set = sort(sortdist(1:adjh));
                            teller(obs_in_set) = teller(obs_in_set) + 1;
                            teller(end) = teller(end) + 1;
                        end
                    else
                        % The observations correponding to the adjh smallest mahalanobis
                        % distances determine the subset for the next iteration.
                        if ~part
                            [~,sortdist] = mahal2(bsxfun(@minus,data,meanvct)*P,sqrt(L),part,fine,final,k,obsingroup);
                        else
                            if final
                                [~,sortdist] = mahal2(bsxfun(@minus,data,meanvct)*P,sqrt(L),part,fine,final,k,obsingroup);
                            elseif fine
                                [~,sortdist] = mahal2(bsxfun(@minus,data(obsingroup{end},:),meanvct)*P,sqrt(L),part,fine,final,k,obsingroup);
                            else
                                [~,sortdist] = mahal2(bsxfun(@minus,data(obsingroup{k},:),meanvct)*P,sqrt(L),part,fine,final,k,obsingroup);
                            end
                        end
                        % Creation of a H-subset.
                        obs_in_set=sort(sortdist(1:adjh));
                        teller(obs_in_set) = teller(obs_in_set) + 1;
                        teller(end) = teller(end) + 1;
                    end
                    [Pfull,~,Lfull,r,~,meanvct] = classSVD(data(obs_in_set,:),1,centered);
                    Pzero=Pfull(:,r+1:p);
                    P=Pfull(:,1:r);
                    L=Lfull(1:r);
                    
                    if r==0
                        rew.center=meanvct;
                        rew.cov=zeros(p,p);
                        rew.flag=(data==data(obs_in_set(1),:));
                        rew.Hsubsets.Hopt=obs_in_set(1);
                        rew.Hsubsets.Hfreq=obs_in_set(1);
                    elseif r < p
                        obj=0;
                        % The adjh-subset is singular. If adjh=h we have an exact fit situation.
                        % If adjh < h we distinguish two cases :
                        %
                        % 1. We haven't found earlier a singular adjh-subset. We first check if
                        %    in the entire set there are h observations that lie on the hyperplane.
                        %    If so, we have an exact fit situation. If not, we stop taking C-steps
                        %    (the determinant is zero which is the lowest possible value) and
                        %    store the results in the appropriate arrays.  We then begin with
                        %    the next trial subsample.
                        %
                        % 2. We have, for the concerning subdataset, already found a singular
                        %    adjh-subset. We then immediately begin with the next trial subsample.

                        if ~part || final || (fine && n==minigr)  %adjh=h: exact fit
                            [obsinplaneInd,count,~]=obsinplane(data,meanvct,Pzero,centered);  
                            if ~centered
                                covar = cov(data(obsinplaneInd,:));
                                meanvct=mean(data(obsinplaneInd,:));  
                            else
                                covar=(data(obsinplaneInd,:)'*data(obsinplaneInd,:))/(length(obsinplaneInd)-1);
                                meanvct=zeros(1,p);
                            end                              
                            [covar,center]=trafo(covar,meanvct,med,mad);    
                            [rew.plane.eigvct,rew.plane.eigval,~] = covsvd(covar,0); 
                            weights(obsinplaneInd)=1;                            
                            [raw.center,raw.cov,rew.center,rew.cov,raw.objective,...
                                raw.wt,rew.flag,rew.method]=displ(2,count,weights,n,p,center,covar,...
                                rew.method, rew.plane);
                            if cor
                                correl=covar./(sqrt(diag(covar))*sqrt(diag(covar))'); 
                                [rew.cor,raw.cor]=deal(correl);
                            end
                            rew.Hsubsets.Hopt=obsinplaneInd;
                            rew.Hsubsets.Hfreq=obsinplaneInd;
                            return
                        elseif found==0 %till now, no singular solution has been found in the subgroup
                            [obsinplaneInd,count,~]=obsinplane(data,meanvct,Pzero,centered);                        
                            if count >= h
                                if ~centered
                                    covar = cov(data(obsinplaneInd,:));
                                    meanvct=mean(data(obsinplaneInd,:));  
                                else
                                    covar=(data(obsinplaneInd,:)'*data(obsinplaneInd,:))/(length(obsinplaneInd)-1);
                                    meanvct=zeros(1,p);
                                end        
                                [covar,center]=trafo(covar,meanvct,med,mad);    
                                [rew.plane.eigvct,rew.plane.eigval,~] = covsvd(covar,0); 
                                weights(obsinplaneInd)=1;
                                [raw.center,raw.cov,rew.center,rew.cov,raw.objective,...
                                    raw.wt,rew.flag,rew.method]=displ(2,count,weights,n,p,center,covar,...
                                    rew.method,rew.plane);
                                if cor
                                    correl=covar./(sqrt(diag(covar))*sqrt(diag(covar))');      
                                    [rew.cor,raw.cor]=deal(correl);                                    
                                end
                                rew.Hsubsets.Hopt=obsinplaneInd;
                                rew.Hsubsets.Hfreq=obsinplaneInd;
                                return                                
                            end
                            found=1;
                            if ~fine
                                Pzeros1{k}=Pzero;
                            else
                                Pzeros=Pzero;
                            end
                            break;
                        else  %found=1: a singular solution is already found in the subgroup
                            adj=0;
                            break;
                        end
                    else
                        obj=prod(L);
                    end
                    % We stop taking C-steps when two subsequent determinants become equal.
                    % We have then reached convergence.
                    if j >= 2 && obj == prevObj
                        break;
                    end
                    prevObj=obj;

                end % C-steps

            end

            % After each iteration, it has to be checked whether the new solution
            % is better than some previous one.  A distinction is made between the
            % different stages of the algorithm:
            %
            %  - Let us first consider the first (second) stage of the algorithm.
            %    We distinguish two cases if the objective value is lower than the largest
            %    value in bobj1 (bobj) :
            %
            %      1. The new objective value did not yet occur in bobj1 (bobj).  We then store
            %         this value, the corresponding mean and covariance matrix at the right
            %         place in resp. bobj1 (bobj), bmean1 (bmean) and bcov1 (bcov).
            %         The objective value is inserted by shifting the greater determinants
            %         upwards. We perform the same shifting in bmean1 (bmean) and bcov1 (bcov).
            %
            %      2. The new objective value already occurs in bobj1 (bobj). A comparison is
            %         made between the new mean vector and covariance matrix and those
            %         estimates with the same determinant. When for an equal determinant,
            %         the mean vector or covariance matrix do not correspond, the new results
            %         will be stored in bobj1 (bobj), bmean1 (bmean) and bcov1 (bcov).
            %
            %    If the objective value is not lower than the largest value in bobj1 (bobj),
            %    nothing happens.
            %
            %  - For the final stage of the algorithm, only the best solution has to be kept.
            %    We then check if the objective value is lower than the till then lowest value.
            %    If so, we have a new best solution. If not, nothing happens.


            if ~final && adj %&& ~ns
                if fine || ~part
                    if obj < max(bobj)
                        [bmean,bP,bL,bobj]=insertion(bmean,bP,bL,bobj,meanvct,P,L,obj,1,eps);
                    end
                else
                    if obj < max(bobj1(k,:))
                        [bmean1,bP1,bL1,bobj1]=insertion(bmean1,bP1,bL1,bobj1,meanvct,P,L,obj,k,eps);
                    end
                end
            end

            if final && obj< bestobj
                % bestset           : the best subset for the whole data.
                % bestobj           : objective value for this set.
                % initmean, initcov : resp. the mean and covariance matrix of this set.
                bestset=obs_in_set;
                bestobj=obj;
                initmean=meanvct;
                if ~centered
                    initcov = cov(data(bestset,:));                    
                else
                    initcov=(data(bestset,:)'*data(bestset,:))/(length(bestset)-1);
                end                   
                raw.initcov = initcov; 
            end

        end % nsamp
        
    end % ngroup


    if part && ~fine
        fine=1;
    elseif (part && fine && ~final) || (~part && ~final)
        final=1;
    else
        final=2;
    end

end % while loop final ~=2

[P,~,L,~,~,cX] = classSVD(data(bestset,:),0,centered);
mah=mahalanobis(bsxfun(@minus,data,cX)*P,zeros(size(P,2),1),'cov',L);
sortmah=sort(mah);

[~,indbestset] = sort(mah(bestset));
sortbestset = bestset(indbestset);
rew.Hsubsets.Hopt = sortbestset;

if ~factor_ind
    factor = sortmah(h)/chi2inv(h/n,p);
else
    factor = sortmah(h)/chi2inv(h/n,p/2);
end

raw.cov=factor*initcov;
% We express the results in the original units.
[raw.cov,raw.center]=trafo(raw.cov,initmean,med,mad);
raw.objective=bestobj*prod(mad)^2;

if cor
    raw.cor=raw.cov./(sqrt(diag(raw.cov))*sqrt(diag(raw.cov))');
end

% the mahalanobis distances are computed without the factor, therefore we
% have to correct for it now.
mah=mah/factor;
raw.rd=sqrt(mah);
weights=mah<=cutoff.md^2;
raw.wt=weights;
[rew.center,rew.cov]=weightmecov(data,weights,centered);
[trcov,trcenter]=trafo(rew.cov,rew.center,med,mad);

% determination of Hfreq:
[~,indobs] = greatsort(teller(1:(end - 1)));
rew.Hsubsets.Hfreq = indobs(1:(h));
if size(rew.Hsubsets.Hfreq,2) == 1
    rew.Hsubsets.Hfreq = rew.Hsubsets.Hfreq';
end

if cor
    rew.cor=rew.cov./(sqrt(diag(rew.cov))*sqrt(diag(rew.cov))');
end

[Pfull,~,r] = covsvd(rew.cov,1);
Pzero=Pfull(:,r+1:p);
if r<p    
    [obsinplaneInd,count,~]=obsinplane(data,rew.center,Pzero,centered);  
    if ~centered
        covar = cov(data(obsinplaneInd,:));
        meanvct=mean(data(obsinplaneInd,:));  
    else
        covar=(data(obsinplaneInd,:)'*data(obsinplaneInd,:))/(length(obsinplaneInd)-1);
        meanvct=zeros(1,p);
    end           
    [covar,center]=trafo(covar,meanvct,med,mad);          
    [rew.plane.eigvct,rew.plane.eigval,~] = covsvd(covar,0);      
    if cor
        correl=covar./(sqrt(diag(covar))*sqrt(diag(covar))');
    end
    rew.method=displrw(count,n,p,center,covar,rew.method,rew.plane,cor,correl);
    rew.flag=weights;
    rew.rd=raw.rd;
else
    mah=mahalanobis(data,rew.center,'cov',rew.cov);
    rew.flag=(mah <= cutoff.md^2);
    rew.rd=sqrt(mah);
end

rew.mahalanobis=sqrt(mahalanobis(data,clmean,'cov',clcov));
rawo=raw;
reso=rew;
if options.classic==1
    classic.cov=clascov;
    classic.center=clasmean;
    classic.md=rew.mahalanobis;
    classic.flag = (classic.md <= cutoff.md);
    if options.cor==1
        classic.cor=clascov./(sqrt(diag(clascov))*sqrt(diag(clascov))');
    end
    classic.class='COV';
else
    classic=0;
end
%assigning the output
rew=struct('center',{trcenter},'cov',{trcov},'cor',{reso.cor},'h',{h},'Hsubsets',{reso.Hsubsets},'alpha',{alfa},...
    'rd',{reso.rd},'flag',{reso.flag},'md',{reso.mahalanobis},'cutoff',{cutoff},...
    'plane',{reso.plane},'method',{reso.method},'class',{'MCDCOV'},'classic',{classic},'X',{xx});
raw=struct('center',{rawo.center},'cov',{rawo.cov},'cor',{rawo.cor},'objective',{rawo.objective},...
    'rd',{rawo.rd},'cutoff',{cutoff},...
    'wt',{rawo.wt},'class',{'MCDCOV'},'classic',{classic},'X',{xx});

if size(data,2)~=2
    rew=rmfield(rew,'X');
    raw=rmfield(raw,'X');
end

try
    if plots && options.classic
        makeplot(rew,'classic',1)
    elseif plots
        makeplot(rew)
    end
catch %output must be given even if plots are interrupted
    %> delete(gcf) to get rid of the menu
end

%-----------------------------------------------------------------------------------------
function [raw_center,raw_cov,center,covar,raw_objective,raw_wt,mcd_wt,method]=displ(exactfit,...
    count,weights,n,~,center,covar,method,plane,varargin)

% Determines some fields of the output argument REW for the exact fit situation.  It also
% displays and writes the messages concerning the exact fit situation.  If the raw MCD
% covariance matrix is not singular but the reweighted is, then the function displrw is
% called instead of this function.

[raw_center,center]=deal(center);
[raw_cov,covar]=deal(covar);
raw_objective=0;
mcd_wt=weights;
raw_wt=weights;

switch exactfit
    case 1
        msg='The scatter matrix of the data is singular.';
    case 2
        msg='The scatter matrix has become singular during the iterations of the MCD algorithm.';
    case 3
        if length(varargin{2})==1
            fmt=('%g');
        else 
            fmt=[repmat('%g, ',1,length(varargin{2})-1) '%g '];
        end
        msg=sprintf('The %g-th order statistic of the absolute deviation of variable(s) ',varargin{1});
        msg=sprintf([msg fmt],varargin{2});
        msg=sprintf([msg ' is zero. ']);
end

msg=sprintf([msg '\nThere are %g observations in the entire dataset of %g observations that lie on the \n'],count,n);
msg=sprintf([msg '%g dimensional '],size(plane.eigvct,2));
if length(plane.eigvct) < 16
    msg=sprintf([msg 'affine subspace defined by the orthogonal basis vector(s) \n\n']);
    msg=sprintf([msg sprintf([repmat('% 13.4g ',1,size(plane.eigvct,2)) '\n'],plane.eigvct')]);
    msg=sprintf([msg '\n\n and passing through the point a= \n\n']);    
    msg=sprintf([msg sprintf('%g  ',center)]);    
else
    msg=sprintf([msg 'affine subspace plane defined by the orthogonal basis vector(s) in \n']);
    msg=sprintf([msg '''rew.plane.eigvct'' and passing through the point a= ''rew.center''.\n']);
end
method=char(method,msg);
disp(method);    
%-----------------------------------------------------------------------------------------
function method=displrw(count,n,p,center,covar,method,plane,cor,correl)

% Displays and writes messages in the case the reweighted robust covariance matrix
% is singular.

msg=sprintf('The reweighted MCD scatter matrix is singular. \n');
msg=sprintf([msg '\nThere are %g observations in the entire dataset of %g observations that lie on the \n'],count,n);
msg=sprintf([msg '%g dimensional '],size(plane.eigvct,2));
if length(plane.eigvct) < 16
    msg=sprintf([msg 'affine subspace defined by the orthogonal basis vector(s) \n\n']);
    msg=sprintf([msg sprintf([repmat('% 13.4g ',1,size(plane.eigvct,2)) '\n'],plane.eigvct')]);
    msg=sprintf([msg '\n\n and passing through the point a= \n\n']);    
    msg=sprintf([msg sprintf('%g  ',center)]);    
else
    msg=sprintf([msg 'affine subspace plane defined by the orthogonal basis vector(s) in \n']);
    msg=sprintf([msg '''rew.plane.eigvct'' and passing through the point a= ''rew.center''.\n']);
end
msg=sprintf([msg '\n\nTheir covariance matrix equals : \n\n']);
msg=sprintf([msg sprintf([repmat('% 13.4g ',1,p) '\n'],covar)]);
if cor
    msg=sprintf([msg '\n\nand their correlation matrix equals : \n\n']);
    msg=sprintf([msg sprintf([repmat('% 13.4g ',1,p) '\n'],correl)]);
end
method=char(method,msg);
disp(method);  

%------------------------------------------------------------------------------------------
function obsingroup=fillgroup(n,group,ngroup,seed)

% Creates the subdatasets.

obsingroup=cell(1,ngroup+1);

jndex=0;
index=zeros(2,sum(group(1:ngroup)));
for k=1:ngroup
    for m=1:group(k)
        [random,seed]=uniran(seed);
        ran=floor(random*(n-jndex)+1);
        jndex=jndex+1;
        if jndex==1
            index(1,jndex)=ran;
            index(2,jndex)=k;
        else
            index(1,jndex)=ran+jndex-1;
            index(2,jndex)=k;
            ii=find(index(1,1:jndex-1) > ran-1+(1:jndex-1),1);
            if ~isempty(ii)
                index(1,jndex:-1:ii+1)=index(1,jndex-1:-1:ii);
                index(2,jndex:-1:ii+1)=index(2,jndex-1:-1:ii);
                index(1,ii)=ran+ii-1;
                index(2,ii)=k;
            end
        end
    end
    obsingroup{k}=index(1,index(2,:)==k);
    obsingroup{ngroup+1}=[obsingroup{ngroup+1},obsingroup{k}];
end

%-----------------------------------------------------------------------------------------
function [ranset,seed]=randomset(tot,nel,seed)

% This function is called if not all (p+1)-subsets out of n will be considered.
% It randomly draws a subsample of nel cases out of tot.

ranset=zeros(1,nel);
for j=1:nel
    [random,seed]=uniran(seed);
    num=floor(random*tot)+1;
    if j > 1
        while any(ranset==num)
            [random,seed]=uniran(seed);
            num=floor(random*tot)+1;
        end
    end
    ranset(j)=num;
end

%-----------------------------------------------------------------------------------------
function [index,seed]=addobs(index,n,seed)
% Extends a trial subsample with one observation.

indexTot = 1:1:n;
indexDiff = setdiff(indexTot,index); %remaining indices

[random,seed]=uniran(seed);
num=floor( random * length(indexDiff) ) +1; %select index
index(end+1)=indexDiff(num);

%-----------------------------------------------------------------------------------------

function mahsort=mahal(dat,meanvct,covmat,part,fine,final,k,obsingroup,~,~,~,~) %#ok<DEFNU>

% Orders the observations according to the mahalanobis distances.

if ~part || final
    [~,ind]=sort(mahalanobis(dat,meanvct,'cov',covmat));
    mahsort=ind;
elseif fine
    [~,ind]=sort(mahalanobis(dat(obsingroup{end},:),meanvct,'cov',covmat));
    mahsort=obsingroup{end}(ind);
else
    [~,ind]=sort(mahalanobis(dat(obsingroup{k},:),meanvct,'cov',covmat));
    mahsort=obsingroup{k}(ind);
end

%-----------------------------------------------------------------------------------------

function [dis,mahsort]=mahal2(score,sca,part,fine,final,k,obsingroup)

% Orders the observations according to the mahalanobis distances for a diagonal
% covariance matrix and zero mean. sca contains the squareroot of the diagonal elements.

if ~part || final
    [dis,ind]=sort(mahalanobis(score,zeros(size(score,2),1),'cov',sca.^2));
    mahsort=ind;
elseif fine
    [dis,ind]=sort(mahalanobis(score,zeros(size(score,2),1),'cov',sca.^2));
    mahsort=obsingroup{end}(ind);
else
    [dis,ind]=sort(mahalanobis(score,zeros(size(score,2),1),'cov',sca.^2));
    mahsort=obsingroup{k}(ind);
end

%------------------------------------------------------------------------------------------

function [covmat,meanvct]=trafo(covmat,meanvct,med,mad)

% Transforms a mean vector and a covariance matrix to the original units.

covmat=bsxfun(@times,bsxfun(@times,covmat,mad),mad');
meanvct=meanvct.*mad+med;

%-----------------------------------------------------------------------------------------
function [bestmean,bestP,bestL,bobj]=insertion(bestmean,bestP,bestL,bobj,meanvct,P,L,obj,row,~)

% Stores, for the first and second stage of the algorithm, the results in the appropriate
% arrays if it belongs to the 10 best results.

insert=1;

equ=find(obj==bobj(row,:));

for j=equ
    if (meanvct==bestmean{row,j})
        if all(P==bestP{row,j})
            if all(L==bestL{row,j})
                insert=0;
            end
        end
    end
end

if insert
    ins=find(obj < bobj(row,:),1);

    if ins==10
        bestmean{row,ins}=meanvct;
        bestP{row,ins} = P;
        bestL{row,ins} = L;
        bobj(row,ins)=obj;
    else
        [bestmean{row,ins+1:10}]=deal(bestmean{row,ins:9});
        bestmean{row,ins}=meanvct;
        [bestP{row,ins+1:10}] = deal(bestP{row,ins:9});
        bestP{row,ins} = P;
        [bestL{row,ins+1:10}] = deal(bestL{row,ins:9});
        bestL{row,ins} = L;
        bobj(row,ins+1:10)=bobj(row,ins:9);
        bobj(row,ins)=obj;
    end

end

%-----------------------------------------------------------------------------------------
function quan=quanf(alfa,n,rk)
quan=floor(2*floor((n+rk+1)/2)-n+2*(n-floor((n+rk+1)/2))*alfa);
%--------------------------------------------------------------------------

function [P,L,r]=covsvd(cov,varargin)

%COVSVD performs the singular value decomposition of a covariance matrix 

if nargin>1
    all=varargin{1};
else
    all=0;
end

[n,p]=size(cov); 

if n==1
    error('The sample size is 1. No SVD can be performed.')
end

[U,S,~]=svd(cov,0); 
S=diag(S);
tol = max([n p])*eps(S(1));
r=sum(S>tol);
if ~all
    L=S(1:r);
    P=U(:,1:r);
else
    L=S;
    P=U;
end

%---------------------------------------------------------------------------

function[indexOnPlane,countOnPlane,distToPlane]=obsinplane(data,mean,zeroEigvct,centered)
    
    if ~centered
        data=bsxfun(@minus,data,mean);
    end
    
    dists=data*zeroEigvct;
    
    distToPlane=sqrt(sum(dists.^2,2));   
    distToPlane=distToPlane';
    indexOnPlane=find(distToPlane < 1e-8);        
    countOnPlane=length(indexOnPlane);        