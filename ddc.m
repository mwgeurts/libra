function result = ddc(X,varargin)

%DDC is a robust method to detect anamolous data cells. For the method to work 
% well, the variables (columns) of the data should be numerical and take on more  
% than a few values. Therefore the method temporarily sets aside columns that  
% do not satisfy this condition, such as binary dummy variables, and  
% performs its computations on the remaining columns. It can handle data  
% with missing values, but the method also sets aside columns and rows with more 
% than a specified amount of missings (typically 20%). The method works best when 
% the remaining variables are approximately gaussian in their center, that is, 
% apart from any outlying values. 
% The method starts by standardizing the data and flagging the cells that 
% stand out in their column. Next, each data cell is estimated based on 
% the unflagged cells in the same row whose column is correlated with 
% the column in question. Finally, a cell for which the observed value  
% differs much from its estimated value is considered anomalous. The method then 
% produces an imputed data matrix, which is like the original one except that 
% cellwise outliers and missing values are replaced by their estimated values.  
% The method can also flag an entire row if it contains too much anomalous behavior.
%
% The DDC method was introduced in:
%   Rousseeuw, P.J., Van den Bossche, W. (2016),
%   Detecting anomalous data cells, arXiv:1601.07251v1

% Note that the function CELLMAP displays a ddc result as a heatmap image.
%
% Required input argument:
%    X : a matrix whose columns represent variables, and rows represent observations.
%
% Optional input arguments:
%       fracNaN :   only consider columns and rows with fewer NaNs (missing
%                   values) than this fraction (percentage).
%   numDiscrete :   a column that takes on numDiscrete or fewer values will
%                   be considered discrete and not used in the analysis.
%     precScale :   only consider columns whose scale is > precScale.
%                   Here scale is measured by the median absolute deviation.
%       tolProb : 	tolerance probability, with default 0.99, which
%                   determines the cutoff values for flagging outliers.
%                   Used in all the steps of the algorithm.
%       corrlim : 	when tring to estimate z_ij from other variables h, we
%                   will only use variables h with abs(corr(j,h)) >= corrlim.
%                   Variables j without any correlated variables h satisfying
%                   this are considered standalone, and treated on their own.
%    combinRule : 	the operation to combine estimates of z_ij coming from
%                   other variables h: can be 'wmean', 'wmedian', 'nanmean', 
%                   'nanmedian'.
%   includeSelf : 	whether or not the combination rule will include the
%                   variable j itself.
%     rowdetect : 	whether the rule for flagging rows is to be applied.
%                   (default = false)
% returnBigXimp : 	if true (default=false), the imputed data matrix Ximp 
%                   in the output will include the rows and columns that 
%                   were not part of the analysis (and can still contain NaNs).
%      dispInfo :   whether additional information should be supplied 
%                   (false = no info, true = info) with default false.
%
% I/O: result=ddc(X,'fracNaN',0.2,'numDiscrete',3,'precScale',1e-12,...
%           'tolProb',0.99,'corrlim',0.5,'combinRule','wmean','includeSelf',true,...
%           'rowdetect',true,'returnBigXimp',false,'dispInfo',false);
%  The user should only give the input arguments that have to change their
%  default value.
%  The name of the input arguments needs to be followed by their value.
%  The order of the input arguments is of no importance.
%  
% Examples: 
%    result=ddc(X,'fracNaN',0.15,'numDiscrete',5,'dispInfo',true)
%    result=ddc(X,'rowdetect',false,'returnBigXimp',true)
%
% The output structure 'result' contains the final results, namely:
%
%    result.colInAnalysis : The column numbers that remain in the analysis 
%                           after preprocessing.
%    result.rowInAnalysis : The row numbers that remain in the analysis 
%                           after preprocessing.
%  result.namesNotNumeric : The names of the non-numeric columns.
%  result.namesCaseNumber : The names of the columns that were identical to
%                           the case number (number of the row).
%       result.namesNAcol : The names of the columns with over
%                           'fracNaN' * 100% NaNs.
%       result.namesNArow : The names of the rows with over
%                           'fracNaN' * 100% NaNs.
%    result.namesDiscrete : The names of the columns with 'numDiscrete' or 
%                           less different values. 
%   result.namesZeroScale : The names of the columns with zero or tiny 
%                           median absolute deviation.
%     result.namesgoodcol : The names of the good, retained columns.
%     result.namesgoodrow : The names of the good, retained rows.
%             result.remX : The remaining part of the X matrix which consist of 
%                           columns colInAnalysis and rows rowInAnalysis.
%                result.Z : The robustly standardized data.
%                result.k : The amount of correlated non-self neigbours that 
%                           has been searched for by the DDC algorithm.  
%            result.ngbrs : The matrix with the k+1 most correlated neigbours 
%                           for each variable in descending order of absolute value.
%                           I.e. row i indicates the k+1 most correlated  
%                           neigbours of variable i; k+1 since the first value 
%                           is the variable number itself which always has 
%                           perfect correlation. 
%          result.robcors : The matrix with the k+1 highest correlation
%                           values (in descending order of absolute value). 
%                           corresponding with the variables in .ngbrs. 
%        result.robslopes : The matrix with the k+1 robust slope values 
%                           corresponding with the variables in .ngbrs. 
%             result.Xest : The estimated data matrix.
%         result.stdResid : The standardized residuals.
%               result.Ti : Vector that indicates the outlyingness of each
%                           row.
%         result.indcells : The index of the outlying cells.
%          result.indrows : The index of the outlying rows.
%           result.indall : The index of all outlying cells (includes the 
%                           cells of row outliers)
%             result.Ximp : The imputed data matrix
%
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at:
%              http://wis.kuleuven.be/stat/robust
%
% Written by Wannes Van den Bossche
% Revisions by Wannes Van den Bossche
% Last Update: 25/04/2016

% Check optional input arguments (are supplied in pairs)
if rem(nargin-1,2)~=0
    error('The number of input arguments should be odd!');
end 

%Assiging default values
default=struct('fracNaN',0.2,'numDiscrete',3,'precScale',1e-12,...
    'tolProb',0.99,'corrlim',0.5,'combinRule','wmean',...
    'includeSelf',true,'rowdetect',true,'returnBigXimp',false,'dispInfo',false);
list=fieldnames(default);
options=default;
IN=length(list);
i=1;
counter=1;

%Reading optional inputarguments
if nargin>2
    % placing inputfields in array of strings
    chklist=cell(1,floor((nargin-1)/2));
    for j=1:nargin-1
        if rem(j,2)~=0
            varargin{j}=strtrim(varargin{j}); %remove white space
            chklist{i}=varargin{j};
            i=i+1;
        end
    end
    % Checking which default parameters have to be changed
    % and keep them in the structure 'options'.
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
end %end reading optional input arguments

fracNaN = options.fracNaN;
numDiscrete = options.numDiscrete;
precScale = options.precScale;
dispInfo = options.dispInfo;
returnBigXimp = options.returnBigXimp;
DDCpars = options;

% Check the data set and set aside columns and rows that do
% not satisfy the conditions:
if verLessThan('matlab','8.2')
    outCheck = checkDataSetPre2013b(X,fracNaN,numDiscrete,precScale,dispInfo);
else
    outCheck = checkDataSet(X,fracNaN,numDiscrete,precScale,dispInfo);
end


% Carry out the actual DetectDeviatingCells algorithm on
% the remaining dataset out1$remX :
outDDC = DDCcore(outCheck.remX,DDCpars);

if ~verLessThan('matlab','8.2')
    if istable(X)
        outDDC.Z = array2table(outDDC.Z,'Rownames',outCheck.namesgoodrow,...
            'VariableNames',outCheck.namesgoodcol);
        outDDC.ngbrs = array2table(outDDC.ngbrs,'Rownames',outCheck.namesgoodcol,...
            'VariableNames',sprintfc('N%d ',1:size(outDDC.ngbrs,2)));
        outDDC.robcors = array2table(outDDC.robcors,'Rownames',outCheck.namesgoodcol,...
            'VariableNames',sprintfc('N%d ',1:size(outDDC.robcors,2)));
        outDDC.robslopes = array2table(outDDC.robslopes,'Rownames',outCheck.namesgoodcol,...
            'VariableNames',sprintfc('N%d ',1:size(outDDC.robslopes,2)));
        outDDC.Xest = array2table(outDDC.Xest,'Rownames',outCheck.namesgoodrow,...
            'VariableNames',outCheck.namesgoodcol);
        outDDC.stdResid = array2table(outDDC.stdResid,'Rownames',outCheck.namesgoodrow,...
            'VariableNames',outCheck.namesgoodcol);
        outDDC.Ti = array2table(outDDC.Ti,'Rownames',outCheck.namesgoodrow,...
            'VariableNames',{'Ti'});
        outDDC.Ximp = array2table(outDDC.Ximp,...
            'Rownames',outCheck.namesgoodrow,'VariableNames',outCheck.namesgoodcol);
        outCheck.remX = array2table(outCheck.remX,...
            'Rownames',outCheck.namesgoodrow,'VariableNames',outCheck.namesgoodcol);
    end
end

if(returnBigXimp)
    Ximp = X;
    Ximp(outCheck.rowInAnalysis,outCheck.colInAnalysis) = outDDC.Ximp;
    outDDC.Ximp = Ximp;
end

result = mergestruct(outCheck,outDDC);

end

function [result]=checkDataSet(X,fracNaN,numDiscrete, precScale, dispInfo)

narginchk(1, 5);
%set default input arguments
switch nargin
    case 1
        fracNaN = 0.2;
        numDiscrete = 3;
        precScale = 1e-12;
    case 2
        numDiscrete = 3;
        precScale = 1e-12;
    case 3
        precScale = 1e-12;
end

if(iscell(X) || isstruct(X)) 
    error('The input data must be a matrix or a data table.')
end
n = size(X,1);
if (n < 3) 
    error(' The input data must have at least 3 rows (cases)')  
end
d = size(X,2);
if(d < 2) 
    error(' The input data must have at least 2 columns (variables)')
end

msg=sprintf('\n');
msg=sprintf([msg 'The input data has %g rows and %g columns. \n'],n,d);

if dispInfo; disp(msg); end; 
    
if(~istable(X))
   namesgoodrow = sprintfc('%g',1:n);
   X = array2table(X,'RowNames',namesgoodrow);
   namesgoodcol = X.Properties.VariableNames;
else
    if (isempty(X.Properties.RowNames))
        namesgoodrow = sprintfc('%g',1:n);
        X.Properties.RowNames = namesgoodrow;
    else
        namesgoodrow = X.Properties.RowNames';
    end
   namesgoodcol = X.Properties.VariableNames;
end

% 1. Deselect categorical column(s)

remX = X; %remX will be the remaining part of the data
colInAnalysis = varfun(@isnumeric,remX);colInAnalysis=colInAnalysis{:,:};
numgoodcol = sum(colInAnalysis);
vecNotNumeric = (colInAnalysis == 0);
numbadcol = sum(vecNotNumeric); %can be 0
namesNotNumeric={};

if(numbadcol > 0) 
    msg=sprintf('\n');
    msg=sprintf([msg 'The input data contained %g non-numeric columns (variables). \n'],numbadcol);
    msg=sprintf([msg 'Their column names are: \n']);
    msg=sprintf([msg '\n']);
    
    namesNotNumeric = remX.Properties.VariableNames(vecNotNumeric);    
    msg=sprintf([msg dispNames(namesNotNumeric)]);
    msg=sprintf([msg '\n']);
    
    if(numgoodcol > 1)
        msg=sprintf([msg 'These columns will be ignored in the analysis. \n']);
        msg=sprintf([msg 'We continue with the remaining %g numeric columns: \n'],numgoodcol);
        remX = remX(:,colInAnalysis);
    else 
      if dispInfo; disp(msg); end; 
      if(numgoodcol == 0); error(' No numeric columns remain, so we stop.');end
      if(numgoodcol == 1); error(' Only 1 numeric column remains, so we stop.');end
    end   
    msg=sprintf([msg '\n']);
    namesgoodcol = remX.Properties.VariableNames;
    msg=sprintf([msg dispNames(namesgoodcol)]);
    
    if dispInfo; disp(msg); end; 
end

% Turn data into a matrix
remX = table2array(remX);
% The row names are namesgoodrow and the column names are namesgoodcol.

% 2. Deselect column(s) containing only the case number

caseNumber = (1:n)';
dists = mean(abs(bsxfun(@minus,remX,caseNumber)));
dists(~isfinite(dists)) = 1;    %takes out NaN, Inf, -Inf
vecbadcol  = (dists == 0);
numbadcol  = sum(vecbadcol);    %can be 0
goodcol    = (vecbadcol == 0);
numgoodcol = sum(goodcol);
namesCaseNumber = {};

if(numbadcol > 0) 
    msg=sprintf('\n');
    msg=sprintf([msg 'The input data contained %g columns that were '...
        'identical to the case number \n'],numbadcol);
    msg=sprintf([msg ' (number of the row).\n']);
    msg=sprintf([msg 'Their column names are: \n']);
    msg=sprintf([msg '\n']);
    
    namesCaseNumber = namesgoodcol(vecbadcol);    
    msg=sprintf([msg dispNames(namesCaseNumber)]);
    msg=sprintf([msg '\n']);
    
    if(numgoodcol > 1)
        msg=sprintf([msg 'These columns will be ignored in the analysis. \n']);
        msg=sprintf([msg 'We continue with the remaining %g numeric columns: \n'],numgoodcol);
        remX = remX(:,goodcol);
        namesgoodcol = namesgoodcol(goodcol);
    else 
      if dispInfo; disp(msg); end; 
      if(numgoodcol == 0); error(' No columns remain, so we stop.');end
      if(numgoodcol == 1); error(' Only 1 column remains, so we stop.');end
    end   
    % Update array colInAnalysis using goodcol:
    colInAnalysis(colInAnalysis==1) = goodcol; % overwrites correctly    
    msg=sprintf([msg '\n']);
    msg=sprintf([msg dispNames(namesgoodcol)]);
    
    if dispInfo; disp(msg); end;
end

 

% 3. Deselect variables with over fracNA% of missing values
%    (e.g. fracNA=0.20). Then update the vector colInAnalysis.

remX(~isfinite(remX)) = NaN; % sets NaN, Inf, -Inf all to NaN
acceptNaN = n*fracNaN;
NaNcounts   = sum(isnan(remX));
goodcol    = (NaNcounts < acceptNaN);
numgoodcol = sum(goodcol);
vecNaNcol   = (goodcol == 0);
numNaNcol   = sum(vecNaNcol);
namesNaNcol = {};

if(numNaNcol > 0)
    msg=sprintf('\n');
    msg=sprintf([msg 'The data contained %g columns with over '...
            '%0.1f percent of NaNs. \n'],numNaNcol,100*fracNaN);
    msg=sprintf([msg 'Their column names are: \n']);
    msg=sprintf([msg '\n']);
    namesNaNcol = namesgoodcol(vecNaNcol);
    msg=sprintf([msg dispNames(namesNaNcol)]);
    msg=sprintf([msg '\n']);

    if(numgoodcol > 1)
        msg=sprintf([msg 'These columns will be ignored in the analysis. \n']);
        msg=sprintf([msg 'We continue with the remaining %g numeric columns: \n'],numgoodcol);
        remX = remX(:,goodcol);
        namesgoodcol = namesgoodcol(goodcol);
    else 
      if dispInfo; disp(msg); end; 
      if(numgoodcol == 0); error(' No columns remain, so we stop.');end
      if(numgoodcol == 1); error(' Only 1 column remains, so we stop.');end
    end   
    % Update array colInAnalysis using goodcol:
    colInAnalysis(colInAnalysis==1) = goodcol; % overwrites correctly    
    msg=sprintf([msg '\n']);
    msg=sprintf([msg dispNames(namesgoodcol)]);
    
    if dispInfo; disp(msg); end; 
end

% 4. Deselect rows with too many NAs.
%    Create the vector rowInAnalysis.

acceptNaN = size(remX,2)*fracNaN;
NaNcounts   = sum(isnan(remX),2);
goodrow    = (NaNcounts < acceptNaN);
numgoodrow  = sum(goodrow);
vecNaNrow   = (goodrow == 0);
numNaNrow    = sum(vecNaNrow );
rowInAnalysis = goodrow;        % in case we need to remove more rows later.
namesNaNrow = {};

if(numNaNrow > 0)
    msg=sprintf('\n');
    msg=sprintf([msg 'The data contained %g rows with over '...
            '%0.1f percent of NaNs. \n'],numNaNrow,100*fracNaN);
    msg=sprintf([msg 'Their row names are: \n']);
    msg=sprintf([msg '\n']);
    namesNaNrow = namesgoodrow(vecNaNrow);
    msg=sprintf([msg dispNames(namesNaNrow)]);
    msg=sprintf([msg '\n']);

    if(numgoodrow > 2)
        msg=sprintf([msg 'These rows will be ignored in the analysis. \n']);
        msg=sprintf([msg 'We continue with the remaining %g rows: \n'],numgoodrow);
        remX = remX(goodrow,:);
        namesgoodrow = namesgoodrow(goodrow);
    else 
      if dispInfo; disp(msg); end; 
      if(numgoodrow == 0); error(' No rows remain, so we stop.');end
      if(numgoodrow == 1); error(' Only 1 row remains, so we stop.');end
      if(numgoodrow == 2); error(' Only 2 rows remain, so we stop.');end
    end    
    msg=sprintf([msg '\n']);
    msg=sprintf([msg dispNames(namesgoodrow)]);
    
    if dispInfo; disp(msg); end; 

end

%   5. Deselect discrete variables, loosely defined as variables that
%      take on numDiscrete or fewer values, such as binary variables.
  
difference  = diff(sort([remX;max(remX,[],1)+1],1),[],1);
difference(difference~=0)=1;
difference(isnan(difference))= 0; % only counts non-NaNs
valueCount = sum(difference);
goodcol = (valueCount > numDiscrete);
numgoodcol  = sum(goodcol);
vecbadcol   = (goodcol == 0);
numbadcol   = sum(vecbadcol);
namesDiscrete = {};

if(numbadcol > 0)
    msg=sprintf('\n');
    msg=sprintf([msg 'The data contained %g discrete columns with '...
            '%g or fewer values. \n'],numbadcol,numDiscrete);
    msg=sprintf([msg 'Their column names are: \n']);
    msg=sprintf([msg '\n']);
    namesDiscrete = namesgoodcol(vecbadcol);
    msg=sprintf([msg dispNames(namesDiscrete)]);
%     msg=sprintf([msg repmat('%s ',1,numbadcol) '\n'],namesDiscrete{:}); 
    msg=sprintf([msg '\n']);

    if(numgoodcol > 1)
        msg=sprintf([msg 'These columns will be ignored in the analysis. \n']);
        msg=sprintf([msg 'We continue with the remaining %g numeric columns: \n'],numgoodcol);
        remX = remX(:,goodcol);
        namesgoodcol = namesgoodcol(goodcol);
    else 
      if dispInfo; disp(msg); end; 
      if(numgoodcol == 0); error(' No columns remain, so we stop.');end
      if(numgoodcol == 1); error(' Only 1 column remains, so we stop.');end
    end   
    % Update array colInAnalysis using goodcol:
    colInAnalysis(colInAnalysis==1) = goodcol; % overwrites correctly    
    msg=sprintf([msg '\n']);
    msg=sprintf([msg dispNames(namesgoodcol)]);
    
    if dispInfo; disp(msg); end; 

end

% 6. Deselect columns for which the median absolute deviation is
%    zero. This is equivalent to saying that 50% or more of its
%    values are equal.

colScale    = 1.4826*mad(remX,1,1);
goodcol     = (colScale > precScale);
numgoodcol  = sum(goodcol);
vecbadcol   = (goodcol == 0);
numbadcol   = sum(vecbadcol); % can be 0
namesZeroScale = {};

if(numbadcol > 0) 
    
    msg=sprintf('\n');
    msg=sprintf([msg 'The data contained %g columns with zero'...
            ' or tiny median absolute deviation. \n'],numbadcol);
    msg=sprintf([msg 'Their column names are: \n']);
    msg=sprintf([msg '\n']);
    namesZeroScale = namesgoodcol(vecbadcol);
    msg=sprintf([msg dispNames(namesZeroScale)]); 
    msg=sprintf([msg '\n']);

    if(numgoodcol > 1)
        msg=sprintf([msg 'These columns will be ignored in the analysis. \n']);
        msg=sprintf([msg 'We continue with the remaining %g numeric columns: \n'],numgoodcol);
        remX = remX(:,goodcol);
        namesgoodcol = namesgoodcol(goodcol);
    else 
      if dispInfo; disp(msg); end; 
      if(numgoodcol == 0); error(' No columns remain, so we stop.');end
      if(numgoodcol == 1); error(' Only 1 column remains, so we stop.');end
    end   
    % Update array colInAnalysis using goodcol:
    colInAnalysis(colInAnalysis==1) = goodcol; % overwrites correctly    
    msg=sprintf([msg '\n']);
    msg=sprintf([msg dispNames(namesgoodcol)]);
    
    if dispInfo; disp(msg); end; 

end

% check whether we have reduced the size of X
if(size(remX,1) < n || size(remX,2) < d)
    msg=sprintf('\n');
    msg=sprintf([msg 'The final data set we will analyze has %g rows '...
        'and %g columns. \n'],size(remX,1),size(remX,2));
    
    if dispInfo; disp(msg); end;     
end

  result=struct('colInAnalysis',{find(colInAnalysis)},...
              'rowInAnalysis',{find(rowInAnalysis)},...
              'namesNotNumeric',{namesNotNumeric},...
              'namesCaseNumber',{namesCaseNumber},...
              'namesNaNcol',{namesNaNcol},...
              'namesNaNrow',{namesNaNrow},...
              'namesDiscrete',{namesDiscrete},...
              'namesZeroScale',{namesZeroScale},...
              'namesgoodcol',{namesgoodcol},...
              'namesgoodrow',{namesgoodrow},...
              'remX',{remX});

end

function [result]=checkDataSetPre2013b(X,fracNaN,numDiscrete, precScale, dispInfo)

narginchk(1, 5);
%set default input arguments
switch nargin
    case 1
        fracNaN = 0.2;
        numDiscrete = 3;
        precScale = 1e-12;
    case 2
        numDiscrete = 3;
        precScale = 1e-12;
    case 3
        precScale = 1e-12;
end

if(iscell(X) || isstruct(X)) 
    error('The input data must be a matrix or a data table.')
end
n = size(X,1);p = size(X,2);
if (n < 3) 
    error(' The input data must have at least 3 rows (cases)')  
end
d = size(X,2);
if(d < 2) 
    error(' The input data must have at least 2 columns (variables)')
end

msg=sprintf('\n');
msg=sprintf([msg 'The input data has %g rows and %g columns. \n'],n,d);

if dispInfo; disp(msg); end; 
    
namesgoodrow = sprintfc('%g',1:n);
namesgoodcol = sprintfc('X%g',1:p);

% 1. Deselect categorical column(s)

remX = X; %remX will be the remaining part of the data

colInAnalysis = zeros(p,1);
for j=1:p
    colInAnalysis(j) = isnumeric(remX(:,j));
end

numgoodcol = sum(colInAnalysis);
vecNotNumeric = (colInAnalysis == 0);
numbadcol = sum(vecNotNumeric); %can be 0
namesNotNumeric={};

if(numbadcol > 0) 
    msg=sprintf('\n');
    msg=sprintf([msg 'The input data contained %g non-numeric columns (variables). \n'],numbadcol);
    msg=sprintf([msg 'Their column names are: \n']);
    msg=sprintf([msg '\n']);
    
    namesNotNumeric = namesgoodcol(vecNotNumeric);    
    msg=sprintf([msg dispNames(namesNotNumeric)]);
    msg=sprintf([msg '\n']);
    
    if(numgoodcol > 1)
        msg=sprintf([msg 'These columns will be ignored in the analysis. \n']);
        msg=sprintf([msg 'We continue with the remaining %g numeric columns: \n'],numgoodcol);
        remX = remX(:,colInAnalysis);
    else 
      if dispInfo; disp(msg); end; 
      if(numgoodcol == 0); error(' No numeric columns remain, so we stop.');end
      if(numgoodcol == 1); error(' Only 1 numeric column remains, so we stop.');end
    end   
    msg=sprintf([msg '\n']);
    namesgoodcol = namesgoodcol(colInAnalysis == 1);
    msg=sprintf([msg dispNames(namesgoodcol)]);
    
    if dispInfo; disp(msg); end; 
end

% The row names are namesgoodrow and the column names are namesgoodcol.

% 2. Deselect column(s) containing only the case number

caseNumber = (1:n)';
dists = mean(abs(bsxfun(@minus,remX,caseNumber)));
dists(~isfinite(dists)) = 1;    %takes out NaN, Inf, -Inf
vecbadcol  = (dists == 0);
numbadcol  = sum(vecbadcol);    %can be 0
goodcol    = (vecbadcol == 0);
numgoodcol = sum(goodcol);
namesCaseNumber = {};

if(numbadcol > 0) 
    msg=sprintf('\n');
    msg=sprintf([msg 'The input data contained %g columns that were '...
        'identical to the case number \n'],numbadcol);
    msg=sprintf([msg ' (number of the row).\n']);
    msg=sprintf([msg 'Their column names are: \n']);
    msg=sprintf([msg '\n']);
    
    namesCaseNumber = namesgoodcol(vecbadcol);    
    msg=sprintf([msg dispNames(namesCaseNumber)]);
    msg=sprintf([msg '\n']);
    
    if(numgoodcol > 1)
        msg=sprintf([msg 'These columns will be ignored in the analysis. \n']);
        msg=sprintf([msg 'We continue with the remaining %g numeric columns: \n'],numgoodcol);
        remX = remX(:,goodcol);
        namesgoodcol = namesgoodcol(goodcol);
    else 
      if dispInfo; disp(msg); end; 
      if(numgoodcol == 0); error(' No columns remain, so we stop.');end
      if(numgoodcol == 1); error(' Only 1 column remains, so we stop.');end
    end   
    % Update array colInAnalysis using goodcol:
    colInAnalysis(colInAnalysis==1) = goodcol; % overwrites correctly    
    msg=sprintf([msg '\n']);
    msg=sprintf([msg dispNames(namesgoodcol)]);
    
    if dispInfo; disp(msg); end;
end

 

% 3. Deselect variables with over fracNA% of missing values
%    (e.g. fracNA=0.20). Then update the vector colInAnalysis.

remX(~isfinite(remX)) = NaN; % sets NaN, Inf, -Inf all to NaN
acceptNaN = n*fracNaN;
NaNcounts   = sum(isnan(remX));
goodcol    = (NaNcounts < acceptNaN);
numgoodcol = sum(goodcol);
vecNaNcol   = (goodcol == 0);
numNaNcol   = sum(vecNaNcol);
namesNaNcol = {};

if(numNaNcol > 0)
    msg=sprintf('\n');
    msg=sprintf([msg 'The data contained %g columns with over '...
            '%0.1f percent of NaNs. \n'],numNaNcol,100*fracNaN);
    msg=sprintf([msg 'Their column names are: \n']);
    msg=sprintf([msg '\n']);
    namesNaNcol = namesgoodcol(vecNaNcol);
    msg=sprintf([msg dispNames(namesNaNcol)]);
    msg=sprintf([msg '\n']);

    if(numgoodcol > 1)
        msg=sprintf([msg 'These columns will be ignored in the analysis. \n']);
        msg=sprintf([msg 'We continue with the remaining %g numeric columns: \n'],numgoodcol);
        remX = remX(:,goodcol);
        namesgoodcol = namesgoodcol(goodcol);
    else 
      if dispInfo; disp(msg); end; 
      if(numgoodcol == 0); error(' No columns remain, so we stop.');end
      if(numgoodcol == 1); error(' Only 1 column remains, so we stop.');end
    end   
    % Update array colInAnalysis using goodcol:
    colInAnalysis(colInAnalysis==1) = goodcol; % overwrites correctly    
    msg=sprintf([msg '\n']);
    msg=sprintf([msg dispNames(namesgoodcol)]);
    
    if dispInfo; disp(msg); end; 
end

% 4. Deselect rows with too many NAs.
%    Create the vector rowInAnalysis.

acceptNaN = size(remX,2)*fracNaN;
NaNcounts   = sum(isnan(remX),2);
goodrow    = (NaNcounts < acceptNaN);
numgoodrow  = sum(goodrow);
vecNaNrow   = (goodrow == 0);
numNaNrow    = sum(vecNaNrow );
rowInAnalysis = goodrow;        % in case we need to remove more rows later.
namesNaNrow = {};

if(numNaNrow > 0)
    msg=sprintf('\n');
    msg=sprintf([msg 'The data contained %g rows with over '...
            '%0.1f percent of NaNs. \n'],numNaNrow,100*fracNaN);
    msg=sprintf([msg 'Their row names are: \n']);
    msg=sprintf([msg '\n']);
    namesNaNrow = namesgoodrow(vecNaNrow);
    msg=sprintf([msg dispNames(namesNaNrow)]);
    msg=sprintf([msg '\n']);

    if(numgoodrow > 2)
        msg=sprintf([msg 'These rows will be ignored in the analysis. \n']);
        msg=sprintf([msg 'We continue with the remaining %g rows: \n'],numgoodrow);
        remX = remX(goodrow,:);
        namesgoodrow = namesgoodrow(goodrow);
    else 
      if dispInfo; disp(msg); end; 
      if(numgoodrow == 0); error(' No rows remain, so we stop.');end
      if(numgoodrow == 1); error(' Only 1 row remains, so we stop.');end
      if(numgoodrow == 2); error(' Only 2 rows remain, so we stop.');end
    end    
    msg=sprintf([msg '\n']);
    msg=sprintf([msg dispNames(namesgoodrow)]);
    
    if dispInfo; disp(msg); end; 

end

%   5. Deselect discrete variables, loosely defined as variables that
%      take on numDiscrete or fewer values, such as binary variables.
  
difference  = diff(sort([remX;max(remX,[],1)+1],1),[],1);
difference(difference~=0)=1;
difference(isnan(difference))= 0; % only counts non-NaNs
valueCount = sum(difference);
goodcol = (valueCount > numDiscrete);
numgoodcol  = sum(goodcol);
vecbadcol   = (goodcol == 0);
numbadcol   = sum(vecbadcol);
namesDiscrete = {};

if(numbadcol > 0)
    msg=sprintf('\n');
    msg=sprintf([msg 'The data contained %g discrete columns with '...
            '%g or fewer values. \n'],numbadcol,numDiscrete);
    msg=sprintf([msg 'Their column names are: \n']);
    msg=sprintf([msg '\n']);
    namesDiscrete = namesgoodcol(vecbadcol);
    msg=sprintf([msg dispNames(namesDiscrete)]);
%     msg=sprintf([msg repmat('%s ',1,numbadcol) '\n'],namesDiscrete{:}); 
    msg=sprintf([msg '\n']);

    if(numgoodcol > 1)
        msg=sprintf([msg 'These columns will be ignored in the analysis. \n']);
        msg=sprintf([msg 'We continue with the remaining %g numeric columns: \n'],numgoodcol);
        remX = remX(:,goodcol);
        namesgoodcol = namesgoodcol(goodcol);
    else 
      if dispInfo; disp(msg); end; 
      if(numgoodcol == 0); error(' No columns remain, so we stop.');end
      if(numgoodcol == 1); error(' Only 1 column remains, so we stop.');end
    end   
    % Update array colInAnalysis using goodcol:
    colInAnalysis(colInAnalysis==1) = goodcol; % overwrites correctly    
    msg=sprintf([msg '\n']);
    msg=sprintf([msg dispNames(namesgoodcol)]);
    
    if dispInfo; disp(msg); end; 

end

% 6. Deselect columns for which the median absolute deviation is
%    zero. This is equivalent to saying that 50% or more of its
%    values are equal.

colScale    = 1.4826*mad(remX,1,1);
goodcol     = (colScale > precScale);
numgoodcol  = sum(goodcol);
vecbadcol   = (goodcol == 0);
numbadcol   = sum(vecbadcol); % can be 0
namesZeroScale = {};

if(numbadcol > 0) 
    
    msg=sprintf('\n');
    msg=sprintf([msg 'The data contained %g columns with zero'...
            ' or tiny median absolute deviation. \n'],numbadcol);
    msg=sprintf([msg 'Their column names are: \n']);
    msg=sprintf([msg '\n']);
    namesZeroScale = namesgoodcol(vecbadcol);
    msg=sprintf([msg dispNames(namesZeroScale)]); 
    msg=sprintf([msg '\n']);

    if(numgoodcol > 1)
        msg=sprintf([msg 'These columns will be ignored in the analysis. \n']);
        msg=sprintf([msg 'We continue with the remaining %g numeric columns: \n'],numgoodcol);
        remX = remX(:,goodcol);
        namesgoodcol = namesgoodcol(goodcol);
    else 
      if dispInfo; disp(msg); end; 
      if(numgoodcol == 0); error(' No columns remain, so we stop.');end
      if(numgoodcol == 1); error(' Only 1 column remains, so we stop.');end
    end   
    % Update array colInAnalysis using goodcol:
    colInAnalysis(colInAnalysis==1) = goodcol; % overwrites correctly    
    msg=sprintf([msg '\n']);
    msg=sprintf([msg dispNames(namesgoodcol)]);
    
    if dispInfo; disp(msg); end; 

end

% check whether we have reduced the size of X
if(size(remX,1) < n || size(remX,2) < d)
    msg=sprintf('\n');
    msg=sprintf([msg 'The final data set we will analyze has %g rows '...
        'and %g columns. \n'],size(remX,1),size(remX,2));
    
    if dispInfo; disp(msg); end;     
end

  result=struct('colInAnalysis',{find(colInAnalysis)},...
              'rowInAnalysis',{find(rowInAnalysis)},...
              'namesNotNumeric',{namesNotNumeric},...
              'namesCaseNumber',{namesCaseNumber},...
              'namesNaNcol',{namesNaNcol},...
              'namesNaNrow',{namesNaNrow},...
              'namesDiscrete',{namesDiscrete},...
              'namesZeroScale',{namesZeroScale},...
              'namesgoodcol',{namesgoodcol},...
              'namesgoodrow',{namesgoodrow},...
              'remX',{remX});

end

function [result]=DDCcore(X,DDCpars)
    
narginchk(1, 2);
%set default input arguments
if (nargin==1)
    DDCpars.precScale = 1e-12;
    DDCpars.tolProb = 0.99;
    DDCpars.corrlim = 0.5;
    DDCpars.combinRule = 'wmean';
    DDCpars.includeSelf = true;
    DDCpars.rowdetect = true;
    DDCpars.dispInfo = false;
end

% Retrieve parameters from the list:

precScale = DDCpars.precScale;
if(DDCpars.tolProb < 0.5); error('tolProb must be >= 0.5');end;
if(DDCpars.tolProb >= 1);error('tolProb must be < 1.0');end;
qCorr = chi2inv(DDCpars.tolProb,2);
% = cutoff used for bivariate correlation.
qRegr = sqrt(chi2inv(DDCpars.tolProb,1));
% = cutoff used for bivariate slope.
qCell = sqrt(chi2inv(DDCpars.tolProb,1));
% = cutoff used for flagging cells.
qRow  = sqrt(chi2inv(DDCpars.tolProb,1)) ;
% = cutoff used for flagging rows.
includeSelf = DDCpars.includeSelf;
corrlim     = DDCpars.corrlim;
combinRule  = DDCpars.combinRule;
rowdetect   = DDCpars.rowdetect;
dispInfo   = DDCpars.dispInfo;

% Now keeping fixed:
robLoc     = @loc1StepM;    % DDCpars.robLoc
robScale   = @scale1StepM;  % DDCpars.robScale  
robSlope   = @slopeMedWLS;  % DDCpars.robSlope  
robCorr    = @corrGKWLS;    % DDCpars.robCorr  
uniDetect  = @limitFilt;    % DDCpars.uniDetect
% robLoc    : robust univariate location estimator.
% robScale  : robust univariate scale estimator.
% robSlope  : robust slope estimator.
% robCorr   : robust correlation estimator.
% uniDetect : 'limitFilt' compares absolute value to a cutoff.
%             'equiGYfilt' means use permutation-equivariant version
%             of the Gervini-Yohai outlier identifier.

% Turn data into a matrix
if(iscell(X) || isstruct(X)) 
    error('The input data must be a matrix or a data table.')
end

if ~verLessThan('matlab','8.2')
    if(istable(X));X=table2array(X);end
end

n = size(X,1);
if (n < 3) 
    error(' The input data must have at least 3 rows (cases)')  
end
d = size(X,2);
if(d < 2) 
    error(' The input data must have at least 2 columns (variables)')
end

msg=sprintf('\n');
msg=sprintf([msg 'The computation started at: %s '],datestr(now));

if dispInfo; disp(msg); end; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     STEP 1: STANDARDIZE DATA      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Robust standardization
locX   = zeros(1,d);
scaleX   = zeros(1,d);
Z = X;
for j= 1:d
    locX(j) = robLoc(X(:,j),[],precScale);
    Z(:,j) = X(:,j) - locX(j);
    scaleX(j) = robScale(Z(:,j),[],[],precScale);
    Z(:,j) = Z(:,j) / scaleX(j);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     STEP 2: UNIVARIATE ANALYSIS   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

% Univariate analysis by column: replace outliers by NAs
U=Z;
for j= 1:d
    U(:,j) = uniDetect(Z(:,j),qCell);
end 
UniIndex = find(isnan(U));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     STEP 3: CALCULATE CORRELATIONS AND SLOPES     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      

k = min(d-1,1000);  % -1 since we do not have to compute the
                    % correlation of a column with itself

% For each column j of U, find the k columns h != j of U that
% it has the highest absolute correlation robCorr with:
U  = [1:d ; U];
bc = zeros(2*k,d);

for j=1:d
    bc(:,j) = kBestCorr(U(:,j), U, robCorr, k, qCorr, precScale);
end

U = U(2:end,:);
ngbrs   = bc(1:k,:)';         % identities of these columns
robcors = bc((k+1):(2*k),:)'; % the correlations

corrweight = abs(robcors);   % should have no NAs
if (corrlim > 0);corrweight(corrweight < corrlim) = 0;end
ngb0 = ngbrs;
ngb0(corrweight == 0) = 0;
U    = [U ; ngb0'];
robslopes = zeros(k,d);
for j=1:d
    robslopes(:,j) = compSlopes(U(:,j), U, robSlope, k, qRegr, precScale);
end
robslopes = robslopes';
U = U(1:n,:);

colStandalone = find(sum(corrweight,2)==0);
colConnected  = find(~ismember((1:d),colStandalone));

indexStandalone = col2cell(colStandalone,n);
indexStandalone = indexStandalone(ismember(indexStandalone,UniIndex));
% = list of flagged cells in standalone variables.

if(includeSelf)
    % if you want to include column j in its own prediction:
    ngbrs      = [(1:d)' ngbrs];
    robcors    = [ones(d,1) robcors];
    corrweight = [ones(d,1) corrweight];
    robslopes  = [ones(d,1) robslopes];
end

numiter = 1; % increase this if you want to iterate
for iter=1:numiter

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     STEP 4 : ESTIMATE CELLS      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     

    Zest = U; % These values will remain for standalone columns.  

    % Estimation for connected variables:
    U = [1:d ;U]; %#ok<AGROW> OK FOR NUMITER=1
    for j=colConnected
        Zest(:,j) = predictCol(U(:,j),U,ngbrs,...
                                    corrweight,robslopes,combinRule); 
    end
    U = U(2:end,:);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     STEP 5 : DESHRINKAGE         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Deshrinkage: rescale Zest[,j] using robSlope of Z[,j] on Zest[,j]
    Zest = [1:d ; Zest]; %#ok<AGROW> OK FOR NUMITER=1
    for j=colConnected
        Zest(:,j) = deShrink(Zest(:,j),Z,robSlope,qRegr,precScale); 
    end
    Zest = Zest(2:end,:);

    % Finally, all NaNs are replaced by zeroes:
    Zest(isnan(Zest)) = 0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     STEP 6 : FLAGGING CELLS      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      

    % Compute cell residuals:
    Zres = Z-Zest; % original minus estimated
    Zres(:,colStandalone) = Z(:,colStandalone); 
    scalest = zeros(1,length(colConnected));colIt = 1;
    for j= colConnected
        scalest(colIt) = robScale(Zres(:,j),[],[],precScale);
        colIt = colIt + 1;
    end
    Zres(:,colConnected) = bsxfun(@rdivide,Zres(:,colConnected),scalest);
    % where the generic R-function scale() scaled these residuals.
    % We don't have to scale the standalone columns, as they
    % were already standardized in the beginning.

    % Next, flag outlying cells by their large residuals:
    indcells = find(abs(Zres) > qCell); % does not flag the NaNs as cells
    U(indcells) = NaN;

end % ends the iteration
clear U; % we no longer use it.

indcells = setdiff(indcells,col2cell(colStandalone,n));
% are the indices of outlying cells in connected variables only
indcells = unique(sort([indcells; indexStandalone]));
% are the indices of both types of outlying cells

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     STEP 7 : FLAGGING ROWS       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      

Ti      = [];
indrows = [];
indall  = indcells;

if(rowdetect == true)
    Ti = zeros(n,1);
    for i = 1:n
        Ti(i,:) = compT(Zres(i,:));    %Does not contain NaNs
    end
    % calculate the test value (outlyingness) of each row:  
    TiC = bsxfun(@minus,Ti,median(Ti)); %centered data
    mad = 1.4826 * median(abs(TiC));
    Ti = bsxfun(@rdivide,TiC,mad);
    indrows = find(isnan(uniDetect(Ti,qRow)) & (Ti>0));
    indall  = unique([indcells; row2cell(indrows,n,d) ]); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     STEP 8: UNSTANDARDIZE AND IMPUTE      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute Xest (storing in the existing matrix Zest to save space)
Zest = bsxfun(@plus,bsxfun(@times,Zest,scaleX),locX);

% compute Ximp (storing in the existing matrix X to save space)
X(indall)  = Zest(indall);  % imputes the outlying cells of X
unavail    = find(isnan(X));
X(unavail) = Zest(unavail); % imputes the missing values of X


msg=sprintf('The computation ended at: %s \n',datestr(now));

if dispInfo; disp(msg); end;

result=struct('Z',{Z},...
          'k',{k},...
          'ngbrs',{ngbrs},...
          'robcors',{robcors},...
          'robslopes',{robslopes},...
          'Xest',{Zest},...
          'stdResid',{Zres},...
          'indcells',{indcells},...
          'Ti',{Ti},...
          'indrows',{indrows},...
          'indall',{indall},...
          'Ximp',{X});

end

% FUNCTIONS TO TRANSFORM INDICES IN A MATRIX:

function [cindex] = col2cell(colNrs,n)
% Transforms column indices to cellwise indices.
% Here colNrs is a vector with column numbers between 1 and d.

narginchk(2, 2);

if isempty(colNrs)
    cindex=[];
else
    cindex = bsxfun(@plus,bsxfun(@times,(colNrs'-1)*n,ones(n,1)),(1:n)');
    cindex = cindex(:);
end
% contains the cell numbers
end

function [rindex] = row2cell(rowNrs,n,d)
% Transforms row indices to cellwise indices.
% Here rowNrs is a vector with row numbers between 1 and n.

narginchk(3, 3);

if isempty(rowNrs)
    rindex=[];
else
    rindex = bsxfun(@plus,bsxfun(@times,rowNrs',ones(d,1)),(0:n:n*(d-1))');
    rindex = rindex(:);
end

end

% UNIVARIATE ESTIMATORS OF LOCATION AND SCALE:

function [mu] = loc1StepM(x,c1,precScale)

% Computes the first step of an algorithm for
% a location M-estimator using the biweight.
% Note that c1 is expressed in unnormalized MAD units.
% In the usual units it is thus c1/qnorm(3/4). 

narginchk(1, 3);

%set default input arguments
switch nargin
    case 1
        c1=3;
        precScale = 1e-12;
    case 2
        precScale = 1e-12;
end

if isempty(c1); c1=3;end
if isempty(precScale); precScale=1e-12;end

x = x(~isnan(x));    % we always take out NaNs
medx = median(x);
ax = abs(x - medx);
denom = c1 * median(ax);
if(denom > precScale)
  ax = ax/denom;
  w = 1 - ax .* ax;
  w = ((abs(w) + w)/2).^2;
  mu = sum(x .* w)/sum(w);
else mu = medx;
end

end

function [sigma] = scale1StepM(x,c2,delta,precScale)

% Computes the first step of an algorithm for
% a scale M-estimator using the Huber rho. 
% Assumes that the univariate data in x has already
% been centered, e.g. by subtracting loc1StepM.
% Note that c2 is expressed in unnormalized MAD units.
% In the usual units it is thus c2/qnorm(3/4).
% If you change c2 you must also change delta.

narginchk(1, 4);

%set default input arguments
switch nargin
    case 1
        c2=2.5;
        delta=0.844472;
        precScale = 1e-12;
    case 2
        delta=0.844472;
        precScale = 1e-12;
    case 3
        precScale = 1e-12;
end

if isempty(c2); c2=2.5;end
if isempty(delta); delta=0.844472;end
if isempty(precScale); precScale=1e-12;end

x = x(~isnan(x));    % we always take out NaNs
n = length(x);
sigma0 = median(abs(x));

if (c2*sigma0 < precScale)
    sigma = sigma0;
else
  x = x/sigma0;
  rho = x.^2;
  rho(rho > c2^2) = c2^2;
  sigma = (sigma0 * sqrt(sum(rho)/(n*delta)));
end

end

function [mu] = wmean(x,w,nanrm) %#ok<DEFNU>

% Computes the weighted mean of vector x using the weights w.
% nanrm=true indicates whether NaN values in x should be ignored (= default).

narginchk(2, 3);

%set default input arguments
if nargin == 2
    nanrm = true;
end

if (nanrm == 1)
    w = w(~isnan(x));
    x = x(~isnan(x));
end

w = w/sum(w);
mu = w * x';

end

% wmedian with interpolation does not correspond to R-version
function [mu] = wmedian(x,w,interpolate,nanrm) %#ok<DEFNU>

% Computes the weighted mean of vector x using the weights w.
% nanrm=true indicates whether NaN values in x should be ignored (= default).
% If true (=default), linear interpolation is used to get a consistent estimate of the 
% weighted median.

narginchk(2, 4);

%set default input arguments
switch nargin
    case 2
        interpolate=true;
        nanrm = true;
    case 3
        nanrm = true;
end

if nargin == 2
    nanrm = 1;
end

if (nanrm == 1)
    w = w(~isnan(x));
    x = x(~isnan(x));
end

[sx,order] = sort(x);
sw = w(order);

midpoint = sum(sw)/2;
if any(w > midpoint)
    mu = x(w == max(w));
else
    csumw = cumsum(sw);
    j = find(csumw<=midpoint,1,'last');
    if csumw(j) == midpoint
         mu = mean(sx([j j+1]));
    else
         if interpolate==false; mu = sx(j+1);
         else mu = sx(j)+ ( (sx(j+1)-sx(j)) / (csumw(j+1)-csumw(j))) * (midpoint - csumw(j) );
         end
    end
end

end

% UNIVARIATE OUTLIER DETETCTION:

function [vout] = limitFilt(v,qCut)

% Detects outliers and sets them to NaN.
% Assumes that the data have already been standardized.

narginchk(2, 2);

vout = v;
vout(abs(v) > qCut) = NaN;

end

function [vout] = rawEquiGYfilt(v,qCut)
    % raw version of equiGYfilt, has no iteration
    
    narginchk(2, 2);
    
    v2 = v.^2;
    n = length(v2);
    u = sort(v2);
    i0 = find(u < qCut); % initial inliers
    n0 = 0;
    if (~isempty(i0))
        i0 = flip(i0);
        i0 = i0(1);
        dn = max(max(chi2cdf(u(i0:n),1) - ( (i0:n)' - 1)/n, 0));
        n0 = round(dn * n);
    end
    vnan = v;
    if (n0 > 0)
        cutoff = u(n - n0 + 1);
        vnan(v2 >= cutoff) = NaN;
    end
    vout=vnan;
end

function [vout] = equiGYfilt(v,qCut,miter)  %#ok<DEFNU>

% Detects outliers and sets them to NA. This is a
% permutation-equivariant version of the GY filter.
% Assumes that the data have already been standardized.

narginchk(2, 3);
    
%set default input arguments
if(nargin==2);miter = 30;end

qCut = qCut.^2;
converge = 0;
iter = 0;  
n = length(v);
observed = ~isnan(v);
vobs = v(observed);
nobs = length(vobs);
id = 1:nobs;
while (converge == 0 && iter < miter)
  iter = iter + 1;
  vobs = rawEquiGYfilt(vobs,qCut);
  id = id(~isnan(vobs));
  if (~any(isnan(vobs)));converge = 1;end
  vobs(isnan(vobs))=[];
end
vobsout = NaN(nobs,1);
vobsout(id) = vobs;
vout = NaN(n,1);
vout(observed) = vobsout;

end

% FUNCTIONS FOR ROBUST CORRELATION AND SLOPE:

function [corr]= corrGKWLS(xcol,qCorr,colj,precScale,varargin)

% Computes a robust correlation between the columns xcol and
% colj using the formula of Gnanadesikan-Kettenring (GK),
% followed by a Weighted Pearson correlation.
% qCorr is a quantile used to determine the weights.
%
% Assumes colj is a vector with same length as xcol
% and that all normalizations have already happened.  

%set default input arguments
if(nargin<=3);precScale=1e-12;end

corr = (scale1StepM(xcol+colj,[],[],precScale,varargin{:})^2 - ...
          scale1StepM(xcol-colj,[],[],precScale,varargin{:})^2)/4;
% This corr should not be NaN since data columns with too many NaNs
% have already been taken out. But just in case someone increases
% the allowed fraction fracNaN too much, we add the precaution:  
if(isnan(corr)) 
    corr = 0;
else
    corr = min(0.99,max(-0.99,corr));
    % Compute reweighted Pearson correlation:
    corrMatInv = diag(ones(1,2)); 
    corrMatInv([2,3]) = -corr; % replace off-diagonal elements
    corrMatInv = corrMatInv/abs(1-corr^2); % divide by determinant
    xtemp = [xcol colj];
    RDi2  = sum((xtemp * corrMatInv) .* xtemp,2);
    rowSel = (RDi2 < qCorr);
    rowSel(isnan(rowSel))=0;
    xcol = xcol(rowSel);
    colj = colj(rowSel);
    corr = sum(xcol.*colj)/sqrt(sum(xcol.^2)*sum(colj.^2));        
end

end

function [kbest]= kBestCorr(colj,U,robCorr,k,qCorr,precScale)

% For a given column colj this computes the k highest absolute 
% correlations (obtained by robCorr) between colj and the columns 
% of the matrix U, and also reports which k columns are selected.
%
% Assumes that the first entry of colj is the number of the column.
% Assumes the remainder of colj is a vector with same length as the 
% remainder of the columns of U, and that all normalizations have
% already happened.
% Assumes k <= d-1 which is ensured before calling this function.
    
narginchk(4, 6);

%set default input arguments
switch nargin
    case 4
        robCorr=@corrGKWLS;
        precScale = 1e-12;
    case 5
        precScale = 1e-12;
end

indj = colj(1);
colj = colj(2:end);    % take out first entry
U    = U(2:end,:);     % take out first row
n    = size(U,1);
p    = size(U,2);
if(~(length(colj) == n)); error(' colj and U do not match');end
if(k > (size(U,2)-1)); error(' k is larger than ncol(U)-1');end
allCorrs = zeros(1,p);
for j = 1:p
    allCorrs(j) = robCorr(U(:,j),qCorr,colj,precScale);
end
[~,selected] = sort(abs(allCorrs),'descend');
notSelf  = ~(selected == indj);  
selected = selected(notSelf); % removes correlation of j with j
% selected has length d-1
selected = selected(1:k); % always works since k <= d-1
% We stack indices and correlations into a single vector to make it
% easier to use apply() later:
kbest=[selected allCorrs(selected)];

end

function [slope]= slopeMedWLS(xcol,qRegr,colj,precScale)

% Computes the slope of a robust regression without intercept
% of the column colj on the column xcol.
% The computation starts by Median regression, followed by
% weighted least squares (WLS) in which the weights are
% determined by the quantile qRegr.
%
% Assumes that colj is a vector with the same length as xcol
% and that both columns are already centered.
   
narginchk(3, 4);

%set default input arguments
if(nargin<=3);precScale=1e-12;end

slopes = colj./xcol;
rawb = median(slopes(~isnan(slopes))); % raw slope
% This raw slope is not capped.
if(isnan(rawb)) 
    slope = 0;
else
    % Now compute weighted LS slope:
    r = colj - rawb.*xcol; % raw residuals
    cutoff = qRegr*scale1StepM(r,[],[],precScale); 
    % cutoff can be zero, which is okay.
    rowSel = ~(abs(r) > cutoff); % selects the inliers
    rowSel(isnan(r))=0;  
    yw = colj(rowSel);
    xw = xcol(rowSel);
    slope = sum(xw.*yw)/(sum(xw.^2)); % slope of colj ~ xcol + 0
    if(isnan(slope)); slope = 0; end
end

end
    
function [slopes]= compSlopes(colj,U,robSlope,k,qRegr,precScale)
    
% For a given column colj this computes the slopes (obtained by
% robSlope) of regressing colj on each of k given columns of 
% the matrix U.
%
% Assumes that the final k entries of colj are the indices of 
% those k columns.
% Assumes the remainder of colj is a vector with same length as the 
% remainder of the columns of U, and that all of these columns are
% already centered.

narginchk(4, 6);

%set default input arguments
switch nargin
    case 4
        robSlope=@slopeMedWLS;
        precScale = 1e-12;
    case 5
        precScale = 1e-12;
end

n    = size(U,1) - k;
ngbr = colj(n+1:end); % has length k, may contain zeroes
colj = colj(1:n);
U    = U(1:n,:);
slopes = zeros(k,1); % initialize
ngbr = ngbr(~(ngbr == 0)); % keep only the necessary neighbors
if(~isempty(ngbr))
  if(length(ngbr)==1)
    b = robSlope(U(:,ngbr),qRegr,colj,precScale);
  else
    b = zeros(length(ngbr),1);
    
    for j=1:length(ngbr) 
        b(j) = robSlope(U(:,ngbr(j)),qRegr,colj,precScale);
    end
  end
  slopes(1:length(b)) = b;
end

end

% ESTIMATING VALUES IN A COLUMN:

function [estcol]= predictCol(colj,U,ngbrs,corrweight,robslopes,combinRule)
    
% Predicts the values in column colj using the set 'ngbrs' of
% columns of U, by applying the combination rule 'combinRule' whose
% inputs are the weights in 'corrweight' and the slopes in 'robslopes'.
%
% Assumes that the first entry of colj is the number of the column.
% Assumes the remainder of colj is a vector with same length as the 
% remainder of the columns of U, and that all of these columns are
% already centered. 
      
narginchk(6, 6);

j    = colj(1);
U    = U(2:end,:);   % take out first row
n = size(U,1);
contributors = (corrweight(j,:) > 0);
if(length(contributors) < 1)
  estcol = 0;
else
  ngb1    = ngbrs(j,contributors);
  slopes1 = robslopes(j,contributors);
  corrwt1 = corrweight(j,contributors);
  ZestAllh = bsxfun(@times,U(:,ngb1),slopes1); % is n by k matrix
  % Predicts column j from each column h, using slope(Z_j ~ Z_h).
  % This array has the estimates from k variables.    
  if ~any(strcmp(combinRule,{'wmean' 'wmedian' 'mean' 'median'}))
      error('Unknown combination rule. Please check function name.');
  else
      fhandle = str2func(combinRule);
      estcol = zeros(n,1);
      for i =1:n
        estcol(i) = fhandle(ZestAllh(i,:),corrwt1);
      end
  end
  
end

end

% DESHRINKAGE:

function [colj]= deShrink(colj,Z,robSlope,qRegr,precScale)
    
% Deshrinks the column colj by multiplying it by the robSlope
% of column j of the matrix Z on colj.
%
% Assumes that the first entry of colj is the number of the column.
% Assumes the remainder of colj is a vector with same length as the 
% columns of Z, and that both columns were already centered.
   
narginchk(4, 5);

%set default input arguments
if(nargin<=4);precScale=1e-12;end

j    = colj(1);
colj = colj(2:end); % the regressor (from Zest)
zj   = Z(:,j);    % column with response (from Z)
a    = robSlope(colj,qRegr,zj,precScale);
colj = a * [0; colj];

end

% CALCULATE TEST VALUE T

function [Tval] = compT(rowi)
    
narginchk(1, 1);

pval = chi2cdf(rowi.^2,1);
pval = pval(~isnan(pval));
Tval = mean(pval) - 0.5;

end

% MERGE STRUCTURES

function sout = mergestruct(varargin)
     %MERGESTRUCT Merge structures with unique fields.

     %   Copyright 2009 The MathWorks, Inc.

     % Start with collecting fieldnames, checking implicitly 
     % that inputs are structures.
     fn = [];
     for k = 1:nargin
        try
            fn = [fn ; fieldnames(varargin{k})]; %#ok<AGROW>
        catch MEstruct
            throw(MEstruct)
       end
    end

    % Make sure the field names are unique.
    if length(fn) ~= length(unique(fn))
        error('mergestruct:FieldsNotUnique',...
            'Field names must be unique');
    end

    % Now concatenate the data from each struct.  Can't use 
    % structfun since input structs may not be scalar.
    c = [];
    for k = 1:nargin
        try
           c = [c ; struct2cell(varargin{k})]; %#ok<AGROW>
        catch MEdata
            throw(MEdata);
        end
    end

    % Construct the output.
    sout = cell2struct(c, fn, 1);
end

% AUXILIARY FUNCTION TO DISPLAY NAMES IN COMMAND WINDOW

function msgStr = dispNames(names)

    numNames = length(names);
    maxNames = max(cellfun('length', names));  %reserved space for each string
    nrPerRow = floor(80/maxNames);             %determine amount of columns
    dispIts = ceil(numNames/nrPerRow);         %amount of rows
    msgStr = '';
    strForm = sprintf('%%-%ds ',maxNames);
    % fill msgStr
    for i = 1:dispIts              
        namesDone = (i-1)*nrPerRow;
        if i==dispIts; reps = numNames-namesDone;
        else reps = nrPerRow; end
        msgStr=sprintf([msgStr repmat(strForm,1,reps) '\n'],names{namesDone+1:namesDone+reps});
    end
end