function result = robpca(x, varargin)

% ROBPCA is a 'ROBust method for Principal Components Analysis'. 
%  It is resistant to outliers in the data. The robust loadings are computed
%  using projection-pursuit techniques and the MCD method. 
%  Therefore ROBPCA can be applied to both low and high-dimensional data sets.
%  In low dimensions, the MCD method is applied (see mcdcov.m).
% 
% The ROBPCA method is described in
%    Hubert, M., Rousseeuw, P.J., Vanden Branden, K. (2005), 
%    "ROBPCA: a new approach to robust principal components analysis", 
%    Technometrics, 47, 64-79.
%
% To select the number of principal components, a robust PRESS (predicted
% residual sum of squares) curve is drawn, based on a fast algorithm for
% cross-validation. This approach is described in:
%
%    Hubert, M., Engelen, S. (2007),
%    "Fast cross-validation of high-breakdown resampling algorithms for PCA",
%    Computational Statistics and Data Analysis, 51, 5013-5024.
%
% ROBPCA is designed for normally distributed data. If the data are skewed,
% a modification of ROBPCA is available, based on the adjusted outlyingness 
% (see adjustedoutlyingness.m). This method is described in:
%
%    Hubert, M., Rousseeuw, P.J., Verdonck, T. (2009),
%    "Robust PCA for skewed data and its outlier map", 
%    Computational Statistics and Data Analysis, 53, 2264-2274.
% 
% In case of missing data, an ER algorithm is applied according to the 
% approach described in:
% 
%    Serneels, S., Verdonck T. (2008),
%    "Principal component analysis for data containing outliers and
%     missing elements", 
%    Computational Statistics and Data Analysis, 52, 1712–1727. 
%
%
% Required input arguments:
%            x : Data matrix (observations in the rows, variables in the
%                columns)
%
% Optional input arguments:
%            k : Number of principal components to compute. If k is missing, 
%                or k = 0, a scree plot and a press curve are drawn which allows
%                you to select the number of principal components. 
%         kmax : Maximal number of principal components to compute (default = 10).
%                If k is provided, kmax does not need to be specified, unless
%                k is larger than 10.                 
%        alpha : (1-alpha) measures the fraction of outliers the algorithm should
%                resist. Any value between 0.5 and 1 may be specified.
%                (default = 0.75)
%            h : (n-h+1) measures the number of outliers the algorithm should 
%                resist. Any value between n/2 and n may be specified.
%                (default = 0.75*n)
%                Alpha and h may not both be specified.
%          mcd : If equal to one: when the number of variables is sufficiently
%                small, the loadings are computed as the eigenvectors of the
%                MCD covariance matrix, hence the function 'mcdcov.m' is
%                automatically called. The number of principal components is
%                then taken as k = rank(x). (default)
%                If equal to zero, the robpca algorithm is always applied.
%        plots : If equal to one, a scree plot, a press curve and a robust
%                score outlier map are drawn (default). If the input argument
%                'classic' is equal to one, the classical plots are drawn as well.
%                If 'plots' is equal to zero, all plots are suppressed (unless
%                k is missing, then the scree plot and press curve are still drawn).
%                See also makeplot.m
%        labsd : The 'labsd' observations with largest score distance are
%                labeled on the outlier map. (default = 3)
%        labod : The 'labod' observations with largest orthogonal distance are
%                labeled on the outlier map. default = 3)          
%      classic : If equal to one, the classical PCA analysis will be performed
%                (see also cpca.m). (default = 0)
%        scree : If equal to one, a scree plot is drawn. If k is given as input,
%                the value is forced to 0, else the default value is one.
%        press : If equal to one, a plot of robust press-values is drawn. 
%                If k is given as input, the value is forced to 0, else the default
%                value is one.
%                If the input argument 'skew' is equal to one, no plot is drawn.
%    robpcamcd : If equal to one (default), the whole robpca procedure is run
%                (computation of outlyingness and MCD).
%                If equal to zero, the program stops after the computation of
%                the outlyingness. The robust eigenvectors then correspond with
%                the eigenvectors of the covariance matrix of the h observations
%                with smallest outlyingness. This yields the same PCA subspace
%                as the full robpca, but not the same eigenvectors and eigenvalues.
%         skew : If equal to zero, the regular robpca is run. If equal to
%                one, the adjusted robpca algorithm for skewed data is run.
%       center : If equal to one, the dataset is considered to be centered 
%                (the location estimate is the origin). Default is zero.
%      missing : If equal to one, the missing entries in the data matrix
%                are estimated. If equal to zero (default) an error message is 
%                provided in case of missing entries.
%    inimisscw : If equal to one (default), the missing values are initially
%                estimated by the columnwise median. If equal to zero they
%                are estimated by the mean value of both the column- and
%                rowwise median.
%         conv : Convergence criterion in case there are missing entries in
%                the data matrix. In each step the relative change of the
%                imputations must exceed conv. (default = 1e-5)
%      maxiter : Number of imputation iterations in case there are missing
%                entries in the data matrix. (default = 100)
%
% I/O: result=robpca(x,'k',k,'kmax',10,'alpha',0.75,'h',h,'mcd',1,'plots',1, ...
%               'labsd',3,'labod',3,'classic',0);
%  The user should only give the input arguments that have to change their
%  default value.
%  The name of the input arguments needs to be followed by their value.
%  The order of the input arguments is of no importance.
%  
% Examples: 
%    result=robpca(x,'k',3,'alpha',0.65,'plots',0)
%    result=robpca(x,'alpha',0.80,'kmax',15,'labsd',5)
%
% The output of ROBPCA is a structure containing 
% 
%    result.P        : Robust loadings (eigenvectors)
%    result.L        : Robust eigenvalues       
%    result.M        : Robust center of the data
%    result.T        : Robust scores 
%    result.k        : Number of (chosen) principal components
%    result.kmax     : Maximal number of principal components
%    result.alpha    : see interpretation in the list of input arguments
%    result.h        : The quantile h used throughout the algorithm
%    result.Xpred    : Data matrix x with missing values replaced by
%                      imputations.
%    result.Hsubsets : A structure that contains H0, H1 and Hfreq:
%                      H0 : The h-subset that contains the h points with the
%                           smallest outlyingness.  
%                      H1 : The optimal h-subset of mcdcov. 
%                      Hfreq : The subset of h points which are the most
%                              frequently selected during the mcdcov algorithm.
%    result.sd       : Robust score distances within the robust PCA subspace
%    result.od       : Orthogonal distances to the robust PCA subspace 
%    result.cutoff   : Cutoff values for the robust score and orthogonal distances
%    result.flag     : The observations whose score distance is larger than
%                      result.cutoff.sd (==> result.flag.sd)or whose orthogonal
%                      distance is larger than result.cutoff.od (==> result.flag.od)
%                      can be considered as outliers and receive a flag equal
%                      to zero (result.flag.all).
%                      The regular observations receive a flag 1.
%    result.class    : 'ROBPCA' 
%    result.classic  : If the input argument 'classic' is equal to one, this
%                      structure contains results of the classical PCA analysis
%                      (see also cpca.m). 
%
% Short description of the method:
%
% Let n denote the number of observations, and p the number of original variables,
% then ROBPCA finds a robust center (p x 1) of the data M and a loading matrix
% P which is (p x k) dimensional. Its columns are orthogonal and define a new
% coordinate system. The scores (n x k) are the coordinates of the centered
% observations with respect to the loadings: T=(X-M)*P. 
% Note that ROBPCA also yields a robust covariance matrix (often singular) which
% can be computed as
%                         cov=out.P*out.L*out.P'
%
% To select the number of principal components, it is useful to look at the
% scree plot which shows the eigenvalues, and the press curve which displays
% a weighted sum of the squared cross-validated orthogonal distances. 
% The outlier map visualizes the observations by plotting their orthogonal
% distance to the robust PCA subspace versus their robust distances 
% within the PCA subspace. This allows to classify the data points into 4 types:
% regular observations, good leverage points, bad leverage points and 
% orthogonal outliers. 
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust
%
% Written by Mia Hubert, Sabine Verboven, Karlien Vanden Branden, Sanne Engelen,
% Tim Verdonck, Johan Van Kerckhoven, Wannes Van den Bossche
% Last Revision: 25/04/2016

% error if an even number of arguments has been given
if rem(nargin-1,2)~=0
    error('The number of input arguments should be odd!');
end

% default values of the parameters
k = 0;
kmax = -1;
alpha = 0.75;
h = -1;
mcd = 1;
plots = 1;
labsd = 3;
labod = 3;
classic = 0;
scree = 1;
press = 1;
robpcamcd = 1;
cutoff = 0.975;
skew = 0;
centered = 0;
missing = 0;
inimisscw = 1;
conv = 1e-5;
maxiter = 50;
dummy=0;

options = struct('k', k, 'kmax', kmax, 'alpha', alpha, 'h', h, 'mcd', mcd, ...
    'plots', plots, 'labsd', labsd, 'labod', labod, 'classic', classic, ...
    'scree', scree, 'press', press, 'robpcamcd', robpcamcd, 'cutoff', cutoff, ...
    'skew', skew,'center',centered,'missing',missing,'inimisscw',inimisscw, 'conv', conv, 'maxiter', maxiter);
list = fieldnames(options);
IN=length(list);
i=1;
counter=1;

% replace default values with values of provided arguments
if nargin > 2
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
        end
        counter=counter+1;
    end
end

missing=options.missing; % missing values are estimated if one
inimisscw=options.inimisscw; % columnwise initial estimation
centered=options.center; % data is considered centered if one

% initializaion
data = x;
[n, p] = size(data);

% perform first imputation
Mindex = find(isnan(data));

if ~isempty(Mindex)
    if missing
        data = inimiss(data,inimisscw);
    else
        error('The data matrix contains missing values. Solve this or set the parameter ''missing'' to one.');
    end
end

% first step: classical PCA to determine the rank of the data matrix
if n < p
    [P1, T1, L1, r, Xc, clm] = kernelEVD(data,0,centered);
else
    [P1, T1, L1, r, Xc, clm] = classSVD(data,0,centered);
end

if r == 0
    error('All data points collapse!')
end

p1 = size(T1, 2);

% default values of the CPCA dependent parameters if not provided by user
if options.kmax==-1
    kmax = min([10, floor(n/2), r]);
    options.kmax=kmax;
end
if options.h==-1
    h = floor(2*floor((n+kmax+1)/2) - n + 2*(n - floor((n+kmax+1)/2))*alpha);
    options.h=h;
end

% handle provided arguments
if nargin > 2
    
    % check whether h and/or alpha have been provided
    dummy = any(strcmp(chklist, 'h')) + 2 * any(strcmp(chklist, 'alpha'));
    
    switch dummy
        case 0  % neither h or alpha have been provided
            if any(strcmp(chklist, 'kmax'))
                kmax = max(min([floor(options.kmax), floor(n/2), r]), 1);
                options.h = min(floor(2*floor((n+kmax+1)/2) - n + ...
                    2*(n - floor((n+kmax+1)/2))*alpha), n);
                % since h depends on kmax, the default h must be updated if
                % kmax is given by the user
            end
        case 3  % both h and alpha have been provided
            error('Both input arguments alpha and h are provided. Only one is required.')
    end
    % the cases where exactly one of both has been provided is treated
    % below
    
    options.h = floor(options.h);
    options.k = floor(options.k);
    % these must be integers
    kmax = max(min([floor(options.kmax), floor(n/2), r]), 1);
    labod = max(0, min(floor(options.labod), n));
    labsd = max(0, min(floor(options.labsd), n));
    k = max(options.k, 0);
    % these are extracted for later convenience
    
    if k > kmax
        k = kmax;
        sprintf(['Attention (robpca.m): the number of principal components, k = ', ...
            num2str(options.k), '\n  is larger than kmax = ', num2str(kmax), ...
            '; k is set to ', num2str(kmax)]);
    end
    % k must be chosen no larger than kmax
    options.k=k;
    options.kmax=kmax;
    
    if dummy == 1 %checking input variable h
        options.alpha = options.h / n;
        
        if k==0
            if options.h < floor((n+kmax+1)/2)
                options.h = floor((n+kmax+1)/2);
                options.alpha = options.h / n;
                mess=sprintf(['Attention (robpca.m): h should be larger than (n+kmax+1)/2.' ...
                    '\n  It is set to its minimum value ', num2str(options.h)]);        
                disp(mess);
            end
        else
            if options.h < floor((n+k+1)/2)
                options.h = floor((n+k+1)/2);
                options.alpha = options.h / n;
                mess=sprintf(['Attention (robpca.m): h should be larger than (n+k+1)/2.' ...
                    '\n  It is set to its minimum value ', num2str(options.h)]);
                disp(mess);
            end
        end
        % number of observations in the H-set must be large enough
         
        if options.h > n
            options.alpha = 0.75;
            
            if k==0
                options.h = floor(2*floor((n+kmax+1)/2) - n + ...
                    2*(n - floor((n+kmax+1)/2))*options.alpha);
            else
                options.h = floor(2*floor((n+k+1)/2) - n + ...
                    2*(n - floor((n+k+1)/2))*options.alpha);
            end
            
            mess=sprintf(['Attention (robpca.m): h should be smaller than n.' ...
                '\n  It is set to its default value ', num2str(options.h)]);
            disp(mess);
        end
        % number of observations in the H-set cannot be larger than n
        
    elseif dummy == 2 %checking input variable alpha
        if options.alpha < 0.5
            options.alpha = 0.5;
            mess=sprintf(['Attention (robpca.m): Alpha should be larger than 0.5.' ...
                '\n  It is set to 0.5.']);
            disp(mess);
        elseif options.alpha > 1
            options.alpha = 0.75;
            mess=sprintf(['Attention (robpca.m): Alpha should be smaller than 1.' ...
                '\n  It is set to 0.75.']);
            disp(mess);
        end
        
        if k==0
            options.h = floor(2*floor((n+kmax+1)/2) - n + ...
                2*(n - floor((n+kmax+1)/2))*options.alpha);
        else
            options.h = floor(2*floor((n+k+1)/2) - n + ...
                2*(n - floor((n+k+1)/2))*options.alpha);            
        end
    end

    dummyh = any(strcmp(chklist, 'h'));
    
    if ~dummyh
        if k==0
            options.h = min(floor(2*floor((n+kmax+1)/2) - n + ...
                2*(n - floor((n+kmax+1)/2))*options.alpha), n);
        else
            options.h = min(floor(2*floor((n+k+1)/2) - n + ...
                2*(n - floor((n+k+1)/2)) * options.alpha), n);            
        end
    end
    
    if k
        options.scree = 0;
        options.press = 0;
    end
    
    if options.skew
        options.press = 0;
    end
    
    if (options.skew || p1 > min(floor(n/5), kmax)) && options.mcd==1
        options.mcd = 0;
        mess=sprintf('Attention (robpca.m): full robpca is performed, input parameter mcd is set to 0.');
        disp(mess);
    end
    
    % no classic analysis can be performed in case of missing entries
    if ~isempty(Mindex) && options.classic==1
        options.classic = 0;
        mess=sprintf(['Attention (robpca.m): classical analysis is not performed ' ...
            'with missing data.']);
        disp(mess);
    end
    
    if ~isempty(Mindex) && options.mcd==1
        options.mcd = 0;
        mess=sprintf('Attention (robpca.m): only full robpca is available with missing data.');
        disp(mess);
    end
    
    alpha = options.alpha;
    h = options.h;
    scree = options.scree;
    press = options.press;
    plots = options.plots;
    mcd = options.mcd;
    robpcamcd = options.robpcamcd;
    cutoff = options.cutoff;
    skew = options.skew;
    classic = options.classic;
    conv = options.conv;
    maxiter = options.maxiter;
else
    if p1 > min(floor(n/5), kmax)
        options.mcd = 0;
        mcd=0;
        mess=sprintf('Attention (robpca.m): n<5*p, full robpca is performed, parameter mcd is set to 0.');
        disp(mess);
    end

    if ~isempty(Mindex) && options.mcd==1
        options.mcd = 0;
        mcd=0;
        mess=sprintf('Attention (robpca.m): only full robpca is available with missing data.');
        disp(mess);
    end
end

%
% MAIN PART
%

if mcd  % if there are many observations, the mcd covariance matrix
            % can be computed, and the first k eigenvalues/vectors are
            % chosen as output
    rew = mcdcov(T1, 'h', h, 'plots', 0,'center',centered);
    [U, S, ~] = svd(rew.cov, 0);
    L = diag(S);
    
    if k
        options.k = min(k, p1);
    else
        bdwidth = 5;
        topbdwidth = 30;
        set(0, 'Units', 'pixels');
        scnsize = get(0, 'ScreenSize');
        pos1 = [bdwidth, 1/3 * scnsize(4) + bdwidth, scnsize(3)/2 - 2*bdwidth, ...
            scnsize(4)/2 - (topbdwidth + bdwidth)];
        pos2 = [pos1(1) + scnsize(3)/2, pos1(2), pos1(3), pos1(4)];
        
        if press
            outcvMcd = cvMcd(T1, p1, rew, h);
            figure('Position', pos1)
            set(gcf, 'Name', 'PRESS curve', 'NumberTitle', 'off');
            plot(1:p1, outcvMcd.press, 'o-')
            title('MCD')
            xlabel('number of LV')
            ylabel('R-PRESS')
        end
        
        if scree
            figure('Position', pos2)
            screeplot(L, 'MCD');
        end
        
        if scree || press
            cumperc = cumsum(L) ./ sum(L);
            disp(['The cumulative percentage of variance explained by the first ', ...
                num2str(kmax), ' components is:']);
            disp(num2str(cumperc'));
            disp(['How many principal components would you like to retain? Max = ', ...
                num2str(kmax), '.']);
            k = input('');
            options.k = k;
        end
        
        % to close the figures.
        if scree
            close
        end
        
        if press
            close
        end
    end
    
    T = (T1 - repmat(rew.center, size(T1, 1), 1)) * U;
    out.M = clm + rew.center * P1';
    out.L = L(1:options.k)';
    out.P = P1 * U(:, 1:options.k);
    out.T = T(:, 1:options.k);
    out.h = h;
    out.k = options.k;
    out.alpha = alpha;
    out.Hsubsets.H0 = rew.Hsubsets.Hopt;
    out.Hsubsets.H1 = [];
    out.Hsubsets.Hfreq = rew.Hsubsets.Hfreq;
    out.skew = skew;
    
% Full RobPca without missing elements
elseif isempty(Mindex) 
    dataSVD.X=data;
    dataSVD.P1=P1;
    dataSVD.T1=T1;
    dataSVD.L1=L1;
    dataSVD.r=r;
    dataSVD.clm=clm;
    out=robpcafull(dataSVD, options, dummy, robpcamcd, 0);

% Missing elements algorithm
elseif ~isempty(Mindex)

    f = 2 * conv;
    iter = 1;
    
    if k==0
        %with PRESS or SCREE
        dataSVD.X=data;
        dataSVD.P1=P1;
        dataSVD.T1=T1;
        dataSVD.L1=L1;
        dataSVD.r=r;
        dataSVD.clm=clm;    
        out=robpcafull(dataSVD, options, dummy, 0, 0);
        k=out.k;
    end
    
    while iter < maxiter && f > conv
        out = robpcafull(data, options, dummy, 0, 1);
        % no final mcd, no plots
        loadings = out.P;
        scores = out.T;
        
        rew = mcdcov(scores, 'plots', 0,'center',centered);
        TTT = (scores - repmat(rew.center, n, 1));
        TTT = diag( (TTT/rew.cov) * TTT');
        cutofft = chi2inv(0.99, k);
        indout = find(TTT > cutofft);
        
        [rmindex, cmindex] = ind2sub([n p], Mindex);
        missout = rmindex;
        
        for q = 1:length(indout)
            indgood = find(not(missout == indout(q)));
            cmindex = cmindex(indgood);
            missout = missout(indgood);
        end
        
        Xpredicted = scores * loadings' + repmat(out.M, n, 1);        
        Mindexnoout = sub2ind([n p], missout, cmindex);
        f = sum((data(Mindexnoout) - Xpredicted(Mindexnoout)).^2) / ...
            length(Mindexnoout);
        data(Mindex) = Xpredicted(Mindex);
        iter = iter + 1;
    end

    if iter==maxiter
        mess=sprintf('Attention (robpca.m): max. missing value imputation iterations reached.');
        disp(mess);
    end

    %include MCD step at the end
    out=robpcafull(data, options, dummy, robpcamcd, 1);
end

% Classical analysis
if classic
    out.classic.P = P1(:, 1:out.k);
    out.classic.L = L1(1:out.k)';
    out.classic.M = clm;
    out.classic.T = T1(:, 1:out.k);
    out.classic.k = out.k;
    out.classic.Xc = Xc;
end

outpr = out;

% Calculation of the distances and flags
if classic
    outDist = CompDist(data, r, outpr, cutoff, robpcamcd, out.classic);
    
    out.classic.sd = outDist.classic.sd;
    out.classic.od = outDist.classic.od;
    out.classic.cutoff.sd = outDist.classic.cutoff.sd;
    out.classic.cutoff.od = outDist.classic.cutoff.od;
    out.classic.class = outDist.classic.class;
    out.classic.flag = outDist.classic.flag;
else
    outDist = CompDist(data, r, outpr, cutoff, robpcamcd);
end

out.sd = outDist.sd;
out.cutoff.sd = outDist.cutoff.sd;
out.od = outDist.od;
out.cutoff.od = outDist.cutoff.od;
out.flag = outDist.flag;
out.class = outDist.class;
out.classic = outDist.classic;
out.alpha = alpha;
out.kmax = kmax;

if ~isempty(Mindex)
    out.Xpred = data;  % X matrix with missings replaced by last imputations
%     out.Xpredict = out.T * out.P' + repmat(out.M, n, 1); 
%     removed since this should give approx. the same result
else
    out.Xpred=[];
end

result = struct('P', {out.P}, 'L', {out.L}, 'M', {out.M}, 'T', {out.T}, ...
    'k', {out.k}, 'kmax', {out.kmax}, 'alpha', {out.alpha}, 'h', {out.h}, ...
    'Xpred', {out.Xpred}, 'Hsubsets', {out.Hsubsets}, 'sd', ...
    {out.sd}, 'od', {out.od}, 'cutoff', {out.cutoff}, 'flag', ...
    {out.flag}, 'class', {out.class},'classic',{out.classic});

if classic
    result.classic = out.classic;
end

if skew
    result.AO = out.AO;
    result.cutoff.AO = out.cutoff.AO;
end

% Plots
try
    if plots && classic
        makeplot(result, 'classic', 1, 'labsd', labsd, 'labod', labod)
    elseif plots
        makeplot(result, 'labsd', labsd, 'labod', labod)
    end
catch  % output must be given even if plots are interrupted
    %> delete(gcf) to get rid of the menu
end


% ----------------------------------------------------------
% Function doing the actual robpca
function result = robpcafull(data, options, dummy, robpcamcd, missingAlg)
% add mcd perhaps

% initialization
k=options.k;
kmax=options.kmax;
h=options.h;
alpha=options.alpha;
skew=options.skew;
cutoff=options.cutoff;
press=options.press;
scree=options.scree;
centered=options.center;

if isstruct(data)
   Xorig=data.X; 
   P1=data.P1;
   T1=data.T1;
   L1=data.L1;
   r=data.r;
   clm=data.clm;
   n = size(T1,1);
   p = size(P1,1);
else
    % first step: classical PCA to determine the rank of the data matrix
    [n,p] = size(data);
    if n < p
        [P1, T1, L1, r, ~, clm] = kernelEVD(data,0,centered);
    else
        [P1, T1, L1, r, ~, clm] = classSVD(data,0,centered);
    end
end

if r == 0
    error('All data points collapse!')
end

% main part
X = T1;
center = clm;
rot = P1;

niter = 100;
seed = 0;

if h ~= n
    if ~skew
        B = twopoints(X, 250, seed,centered);  % n * ri
        Bnorm = zeros(1,size(B, 1));
        
        for i = 1:size(B, 1)
            Bnorm(i) = norm(B(i, :), 2);
        end
        
        Bnormr = Bnorm(Bnorm > 1.e-12);  % ndirect * 1
        B = B(Bnorm > 1.e-12, :);        % ndirect * n
        A = diag(1 ./ Bnormr) * B;       % ndirect * n
        % projected points in columns
        Y = X * A';                     % n * ndirect
        m = length(Bnormr);
        Z = zeros(n, m);
        
        for i = 1:m
            [tmcdi, smcdi, weights] = unimcd(Y(:, i), h, centered);
            if smcdi < 1.e-12
                r2 = rank(X(weights, :));
                
                if r2 == 1
                    error(['At least', num2str(sum(weights)), ...
                        'observations are identical']);
                end
            else
                Z(:,i)=abs(Y(:,i)-tmcdi)/smcdi;
            end
        end
        
        d = max(Z, [], 2);
    else  % adjusted robpca for skewed data
        outAO = adjustedoutlyingness(X, 'ndir', min(250*p, 2500));
        d = outAO.adjout;
    end
    
    [~, is] = sort(d);
    Xh = X(is(1:h), :);  % Xh contains h (good) points out of Xcentr
    [P2, ~, L2, r2, ~, clmX] = classSVD(Xh,0,centered);
    out.Hsubsets.H0 = is(1:h);
    Tn = (X - repmat(clmX, n, 1)) * P2;
    
else %h==n
    P2 = eye(r);
    Tn = X;
    L2 = L1;
    r2 = r;
    out.Hsubsets.H0 = 1:n;
    clmX = zeros(1, size(X, 2));
end

if size(out.Hsubsets.H0, 2) == 1
    out.Hsubsets.H0 = out.Hsubsets.H0';
end

% dim(P2) = r x r2
L = L2;
if missingAlg; k = min(k, r2);end

if ~missingAlg
    
    kmax = min(kmax, r2);
    % Initialise plots, choice of k
    bdwidth = 5;
    topbdwidth = 30;
    set(0, 'Units', 'pixels');
    scnsize = get(0, 'ScreenSize');
    pos1 = [bdwidth, scnsize(4)/3 + bdwidth, scnsize(3)/2 - 2*bdwidth, ...
        scnsize(4)/2 - (topbdwidth + bdwidth)];
    pos2 = pos1;
    pos2(1) = pos2(1) + scnsize(3)/2;

    % PRESS plot
    if press
        disp('The robust press curve based on cross-validation is now computed.')
        outprMCDkmax = projectMCD(Tn, L, kmax, h, niter, rot, P1, P2, ...
            center, cutoff,centered);

        outprMCDkmax.Hsubsets.H0 = out.Hsubsets.H0;
        outpress = cvRobpca(Xorig, kmax, outprMCDkmax, 0, h);
        figure('Position', pos1)
        set(gcf, 'Name', 'PRESS curve', 'NumberTitle', 'off');
        plot(1:kmax, outpress.press, 'o-')
        title('ROBPCA')
        xlabel('Number of LV')
        ylabel('R-PRESS')
    end

    % scree plot
    if scree
        figure('Position', pos2)
        screeplot(L(1:kmax), 'ROBPCA')
    end

    if scree || press
        cumperc = (cumsum(L(1:kmax)) ./ sum(L))';
        disp(['The cumulative percentage of variance explained by the first ', ...
            num2str(kmax), 'components is:']);
        disp(num2str(cumperc));
        disp(['How many principal components would you like to retain? Max = ', ...
            num2str(kmax), '.'])
        k = input('');
        k = max(min(min(r2, k), kmax), 1);

        % we compute again the robpca results for a specific k value. alpha
        % and h can change again, because until now they were based on the
        % kmax value.
        if dummy == 2
            h = floor(2 * floor((n+k+1)/2) - n + 2 * (n - ...
                floor((n+k+1)/2)) * alpha);
        elseif dummy ~= 1
            h = floor(2 * floor((n+k+1)/2) - n + 2 * (n - ...
                floor((n+k+1)/2)) * alpha);
        end

        % closing the plots
        if scree
            close
        end

        if press
            close
        end
    else
        k = min(min(r2, k), kmax);
    end

end

if k ~= r %extra reweighting step
    XRc = X - repmat(clmX, n, 1);
    Xtilde = XRc * P2(:, 1:k) * P2(:, 1:k)';
    Rdiff = XRc - Xtilde;
    
    odh = zeros(n, 1);
    for i = 1:n
        odh(i, 1) = norm(Rdiff(i, :));
    end
    
    if ~skew
        [m, s] = unimcd(odh .^(2/3), h,centered);
        cutoffodh = sqrt(norminv(cutoff, m, s) .^ 3);
    else  % adjusted robpca for skewed data
        mcodh = mc(odh);
        if mcodh > 0
            cutoffodh = prctilde(odg, 75) + 1.5 * exp(3*mcodh) * iqr(odh);
        else
            cutoffodh = prctilde(odg, 75) + 1.5 * iqr(odh);
        end
        ttup = sort(- odh(odh < cutoffodh));
        cutoffodh = - ttup(1);
    end
    indexset = find(odh <= cutoffodh)';
    [P2, ~, Lh, rh, ~, clmX] = classSVD(X(indexset, :),0,centered);
    
    if k > rh
        k = rh;
    end
end

center = center + clmX * rot';
rot = rot * P2(:, 1:k);
Tn = (X - repmat(clmX, n, 1)) * P2;

% if only the subspace is important, not the PC themselves: do not perform
% MCD anymore
if ~robpcamcd
    out.P = rot;
    out.T = Tn(:, 1:k);
    out.M = center;
    out.L = Lh;
    out.k = k;
    out.h = h;
else
    % projection & mcd
    if ~skew
        outpr = projectMCD(Tn, L, k, h, niter, rot, P1, P2, center, cutoff,centered);

        out.Hsubsets.H1 = outpr.Hsubsets.H1;
        out.Hsubsets.Hfreq = outpr.Hsubsets.Hfreq;
    else
        outpr = projectAO(Tn, k, h, rot, center);
        
        out.AO = outpr.AO;
        out.cutoff.AO = outpr.cutoff.AO;
    end
    
    out.T = outpr.T;
    out.P = outpr.P;
    out.M = outpr.M;
    out.L = outpr.L;
    out.k = k;
    out.h = h;
end

result = struct('P', {out.P}, 'L', {out.L}, 'M', {out.M}, 'T', {out.T}, ...
    'k', {out.k}, 'h', {out.h}, 'Hsubsets', {out.Hsubsets}, 'skew',{skew});

if skew
    result.AO = out.AO;
    result.cutoff.AO = out.cutoff.AO;
end


% ----------------------------------------------------------
function outprMCD = projectMCD(Tn,L,k,h,niter,rot,P1,P2,center,cutoff,centered)

% this function performs the last part of ROBPCA when k is determined.
% input : 
%   Tn : the projected data
%   L  : the matrix of the eigenvalues
%   k  : the number of components
%   h  : lower bound for regular observations
%   niter : the number of iterations
%   rot : the rotation matrix
%   P1, P2: the different eigenvector matrices after each transformation
%   center : the classical center of the data

X2 = Tn(:, 1:k);
n = size(X2, 1);
rot = rot(:, 1:k);
% first apply c-step with h points from first step,i.e. those that
% determine the covariance matrix after the c-steps have converged.
mah = mahalanobis(X2, zeros(size(X2, 2), 1), 'cov', L(1:k));
oldobj = prod(L(1:k));
P4 = eye(k); 
korig = k;

j=1;
while j <= niter
    [~, is] = sort(mah);
    Xh = X2(is(1:h), :);
    [P, ~, L, r3, ~, clmX] = classSVD(Xh,0,centered);
    obj = prod(L);
    X2 = (X2 - repmat(clmX, n, 1)) * P;
    center = center + clmX*rot';
    rot = rot * P;
    mah = mahalanobis(X2, zeros(size(X2, 2), 1), 'cov', diag(L));
    P4 = P4 * P;
    
    if ((r3 == k) && (abs(oldobj - obj) < 1.e-12))
        break;
    else
        oldobj = obj;
        j = j + 1;
        
        if r3 < k
            j = 1;
            k = r3;
        end
    end
end

% dim(P4): k x k0 with k0 <= k but denoted as k
% dim X2: n x k0
% perform mcdcov on X2
[zres,zraw] = mcdcov(X2, 'plots', 0, 'ntrial', 250, 'h', h, 'center', centered);
out.resMCD = zres;

if zraw.objective < obj
    z = zres;
    out.Hsubsets.H1 = zres.Hsubsets.Hopt;
else
    sortmah = sort(mah);
    
    if h == n
        factor = 1;
    else
        factor = sortmah(h) / chi2inv(h/n, k);
    end
    
    mah = mah / factor;
    weights = mah <= chi2inv(cutoff, k);
    [center_noMCD, cov_noMCD] = weightmecov(X2, weights, centered);
    mah = mahalanobis(X2, center_noMCD, 'cov', cov_noMCD);
    z.flag = (mah <= chi2inv(cutoff, k));
    z.center = center_noMCD;
    z.cov = cov_noMCD;
    out.Hsubsets.H1 = is(1:h);
end

covf = z.cov;
centerf = z.center;
[P6, L] = eig(covf);
[L, I]=greatsort(diag(real(L)));
P6 = P6(:, I);
out.T = (X2 - repmat(centerf, n, 1)) * P6;
P = P1 * P2;
out.P = P(:, 1:korig) * P4 * P6;
centerfp = center + centerf*rot';
out.M = centerfp;
out.L = L';
out.k = k;
out.h = h;

% creation of Hfreq
out.Hsubsets.Hfreq = zres.Hsubsets.Hfreq; %used to be zres.Hsubsets.Hfreq(1:h) but causes error in exact fit situation
   
outprMCD = out;


%----------------------------------------------------------------------------
function outprAO = projectAO(Tn, k, h, rot, center)

% this function performs the last part of ROBPCA when k is determined.
% input :
%   Tn : the projected data
%   k  : the number of components
%   h  : lower bound for regular observations
%   rot : the rotation matrix
%   center : the classical center of the data

X2 = Tn(:, 1:k);
[n, p] = size(X2);
outAO = adjustedoutlyingness(X2, 'ndir', min(250*p, 2500));
AO = outAO.adjout;
cutoffAO = outAO.cutoff;
indexset = find(AO <= cutoffAO)';
%SVD
[P6, ~, L6, ~, ~, clmX6] = classSVD(X2(indexset, :));
out.T = (X2 - repmat(clmX6, n, 1)) * P6;
out.P = rot * P6;
out.M = center + clmX6 * rot';
out.L = L6;
out.k = k;
out.h = h;
out.AO = AO;
out.cutoff.AO = cutoffAO;
outprAO = out;


%--------------------------------------------------------------------------------
function outDist = CompDist(data, r, out, cutoff, robpcamcd, classic)

% Calculates the distances.
% input: data : the original data
%           r : the rank of the data
%        out is a structure that contains the results of the PCA.
%        classic: an optional structure:
%               classic.P1
%               classic.T1
%               classic.L1
%               classic.clm
%               classic.Xc

if nargin < 6
    options.classic = 0;
else
    options.classic = 1;
end

n = size(data, 1);
k = out.k;
skew = out.skew;

% Computing distances 
% Robust score distances in robust PCA subspace
if robpcamcd
    if skew == 0
        out.sd = sqrt(mahalanobis(out.T, zeros(size(out.T, 2), 1), 'cov', out.L))';
        out.cutoff.sd = sqrt(chi2inv(cutoff, out.k));
    else
        out.sd = out.AO;
        out.cutoff.sd = out.cutoff.AO;
    end
else
    out.sd = zeros(n, 1);
    out.cutoff.sd = 0;
end

% Orthogonal distances to robust PCA subspace
XRc = data - repmat(out.M, n, 1);
Xtilde = out.T * out.P';
Rdiff = XRc - Xtilde;

for i = 1:n
    out.od(i, 1) = norm(Rdiff(i, :));
end

% Robust cutoff-value for the orthogonal distance
if k ~= r
    if skew == 0
        [m, s] = unimcd(out.od.^(2/3), out.h);
        out.cutoff.od = sqrt(norminv(cutoff, m, s).^3);
    else
        mcod = mc(out.od);
        
        if mcod > 0
            out.cutoff.od = prctile(out.od, 75) + 1.5 * exp(3*mcod) * iqr(out.od);
        else
            out.cutoff.od = prctile(out.od,75) + 1.5 * iqr(out.od);
        end
        
        ttup = sort(-out.od(out.od < out.cutoff.od));
        out.cutoff.od = -ttup(1);
    end
else
    out.cutoff.od = 0;
end

if options.classic == 1
    % Mahalanobis distance in classical PCA subspace
    Tclas = classic.Xc * classic.P(:, 1:out.k);
    out.classic.sd = sqrt(mahalanobis(Tclas, zeros(size(Tclas, 2), 1), ...
        'invcov', 1./classic.L(1:out.k)))';
    % Orthogonal distances to classical PCA subspace
    Xtilde = Tclas * classic.P(:, 1:out.k)';
    Cdiff = classic.Xc - Xtilde;
    
    for i = 1:n
        out.classic.od(i, 1) = norm(Cdiff(i, :));
    end
    
    out.classic.cutoff.sd = sqrt(chi2inv(cutoff, out.k)); % should be defined after od to have the output in the correct order
    
    % Classical cutoff-values
    if k ~= r
        m = mean(out.classic.od.^(2/3));
        s = sqrt(var(out.classic.od.^(2/3)));
        out.classic.cutoff.od = sqrt(norminv(cutoff, m, s)^3); 
    else
        out.classic.cutoff.od = 0;
    end
    
    out.classic.cutoff.sd = sqrt(chi2inv(cutoff, out.k));
    out.classic.flag.od = (out.classic.od <= out.classic.cutoff.od);
    out.classic.flag.sd = (out.classic.sd <= out.classic.cutoff.sd);
    out.classic.flag.all = (out.classic.flag.od) & (out.classic.flag.sd);
    out.classic.class = 'CPCA';
else
    out.classic = 0;
end   

out.class = 'ROBPCA';

if k ~= r
    out.flag.od = (out.od <= out.cutoff.od);
    out.flag.sd = (out.sd <= out.cutoff.sd);
    out.flag.all = (out.flag.od) & (out.flag.sd);
else
    out.flag.od = (out.od <= out.cutoff.od);
    out.flag.sd = (out.sd <= out.cutoff.sd);
    out.flag.all = (out.sd <= out.cutoff.sd);
end

outDist = out;
%-----------------------------------------------------------------
function [X]=inimiss(X,varargin)

%**************************************************************************
%function [Xn,y]=inimiss(X)
%AIM:         initial estimation of missing elements
%INPUT:       X-data matrix X (objects x variables) with missing elements
%             cw: if 1 uses columnwise median to replace missing values (=default)
%                 if 0 uses mean of column- and rowwise medians 
%OUTPUT:      X-data matrix with initial estimates of missing elements
%**************************************************************************

if nargin>1
    cw=varargin{1};
else
    cw=1;
end

[m,~]=size(X);

mX=nanmedian(X);
stX=nanmedian(abs(X-repmat(mX,m,1))); %nanmad
X=(X-repmat(mX,m,1))./repmat(stX,m,1);

[i,j]=find(isnan(X));
meanc=nanmedian(X,1);

if cw
	for k=1:length(i),
        X(i(k),j(k))=meanc(j(k));
	end 
else
    meanr=nanmedian(X,2);
    for k=1:length(i),
        X(i(k),j(k))=(meanr(i(k))+meanc(j(k)))/2;
    end
end
X=X.*repmat(stX,m,1)+repmat(mX,m,1);
