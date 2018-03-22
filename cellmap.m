function [hImage, hText, hXText] = cellmap(D, R, indCells, indRows, varargin)
%
% CELLMAP displays a ddc result as a heatmap image.
%
% USAGE:
% [hImage, hText, hTick] = cellmap(D, R, indCells, indRows, labelsx, labelsy, varargin)
%
% Required input argument:
%
%            D : the data matrix
%            R : The matrix of cell residuals
%     indCells : vector which contains indices of outlying cells
%     indRows  : vector which contains indices of outlying rows
%
% Optional input arguments:
% 
%     labelsx  : can be either a numeric vector or cell array of strings that
%                represent the columns of the matrix. If either is not specified 
%                or empty, no labels will be drawn.
%     labelsy  : can be either a numeric vector or cell array of strings that
%                represent the rows of the matrix. If either is not specified 
%                or empty, no labels will be drawn.
%        title : title of the plot
%       xtitle : title for the x-axis
%       ytitle : title for the y-axis
%     showVals : can either be: 1 (or true) in which case the values of the 
%                data matrix D will be displayed in each square, or a format 
%                string in which case the decimal matrix values will be displayed 
%                formatted according to the string specified. (e.g.'%0.1f')   
%   xblocksize : amount of blocks to combine in the x direction
%   yblocksize : amount of blocks to combine in the y direction
%    autoLabel : automatically combines labels of blocks (default=false)
%    TickAngle : Angle of rotation of tick labels on x-axis. (Default: 0)
%     FontSize : The initial fontSize of the text labels on the image. As
%                the image size is scaled the fontSize is shrunk appropriately.      
%     ColorBar : Display colorbar. The corresponding value parameter should
%                be either logical 1 or 0 or a cell array of any additional parameters
%                you wish to pass to the colorbar function (such as location)
% TickFontSize : Font size of the X and Y tick labels. Default value is the 
%                default axes font size, usually 10. Set to a lower value 
%                if many tick labels are being displayed. 
% TickTexInterpreter : Set to 1 or true to render tick labels using a TEX
%                interpreter. For example, '_b' and '^o' would be rendered 
%                as subscript b and the degree symbol with the TEX interpreter. 
%                This parameter is only available in MATLAB R2014b and above 
%                (Default: false)
%
% OUTPUTS:
% * hImage: handle to the image object
% * hText : handle to the text objects (empty if no text labels are drawn)
% * hTick : handle to the X-tick label text objects if tick angle is not 0
%           (empty otherwise)
%
% EXAMPLE:
% 
% cellmap(R, D, indCells, indRows, 'labelsx', labelsx, 'labelsy', labelsy...
%         'FontSize', 10, 'textFormat', %0.1f, 'TickAngle',45)
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at:
%              http://wis.kuleuven.be/stat/robust
%
% Written by Wannes Van den Bossche
% Revisions by Wannes Van den Bossche
% Last Update: 25/04/2016

% Handle missing inputs
if nargin < 1, error('Cellmap requires at least one input argument'); end
if nargin < 2, D = []; end
if nargin < 3, indCells = []; end
if nargin < 4, indRows = []; end

% Parse parameter/value inputs
pars = parseInputs(varargin{:});

if ~verLessThan('matlab','8.2')
    if (istable(R));R=table2array(R);end;
    if (istable(D));D=table2array(D);end;
end

n = size(R,1);
d = size(R,2);

% check input arguments 
blockMap = false;
if (pars.xblocksize > 1 || pars.yblocksize > 1 )
    blockMap = true;
    if (pars.xblocksize > d);error('Input argument xblocksize cannot be larger than d');end;
    if (pars.yblocksize > n);error('Input argument yblocksize cannot be larger than n');end;
    if (pars.showVals);warning(['The option showVals=TRUE cannot be combined with' ... 
                                        'xblocksize or yblocksize greater than 1']);end;
    pars.showVals = false;

    if (~all(size(R) == size(D)) && ~blockMap);error('Dimensions of D and R must match');end;
    if (pars.autoLabel && ~isempty(pars.labelsx) && length(pars.labelsx)~= d); warning('Dimension of labelsx does not match with d');end;
    if (pars.autoLabel && ~isempty(pars.labelsy) && length(pars.labelsy)~= n); warning('Dimension of labelsy does not match with n');end;
end

% create X: Matrix indicating type of cell
% 0=yellow, 1=blue, 2=red, 3=black, 4=white
X = zeros(n,d);
X(indRows,:) = 3;
pcells = intersect(indCells,find(R>=0));
ncells = intersect(indCells,find(R< 0));
X(ncells) = 1;
X(pcells) = 2;
if (~blockMap);X(isnan(D)) = 4;end;

if (blockMap) 

    n = floor(n/pars.yblocksize);
    d = floor(d/pars.xblocksize);

    % if autolabel==false, labels{x,y} will be used for the blocks.
    if (pars.autoLabel) % automatically combine labels for blocks
        labx = pars.labelsx;
        laby = pars.labelsy;
        pars.labelsx = cell(1,d);
        for ind = 1:d
          pars.labelsx(ind) = cellstr(strjoin([labx((1+((ind-1)*pars.xblocksize))) '-' ...
                               labx((ind*pars.xblocksize))]));
        end
        pars.labelsy = cell(1,n);
        for ind = 1:n
          pars.labelsy(ind) = cellstr(strjoin([laby((1+((ind-1)*pars.yblocksize))) '-' ...
                               laby((ind*pars.yblocksize))]));
        end
    end

    X = squeeze(X ,n, d, pars.xblocksize, pars.yblocksize);
    textmat=[];
    colormapStr = 'ddcblock';
    ColorLevels = 64;
else

    % Transform data matrix D in text
    if ischar(pars.showVals) % If a format string, convert to text with specified format
        formatText = pars.showVals;
        textmat = mat2cstr(D,formatText);
    elseif pars.showVals==1  % If 1, convert to text with default format
        formatText = '%0.1f';
        textmat = mat2cstr(D,formatText);
    else                     % Do not show data
        textmat=[];
    end
    colormapStr = 'ddc';
    ColorLevels = [];
end

[hImage, hText, hXText] = heatmap(X, pars.labelsx, pars.labelsy, textmat,'Colormap',colormapStr,...
    'TickAngle',pars.TickAngle,'FontSize',pars.FontSize,'ColorBar',pars.ColorBar,...
    'TickFontSize',pars.TickFontSize,'TickTexInterpreter',pars.TickTexInterpreter,...
    'Title',pars.Title,'xTitle',pars.xTitle,'yTitle',pars.yTitle,'ColorLevels',ColorLevels);


end

% Parse PV inputs & return structure of parameters
function param = parseInputs(varargin) 

p = inputParser;
if ~verLessThan('matlab','8.2')
    p.addParameter('TickAngle',90);
    p.addParameter('FontSize',[]);
    p.addParameter('showVals',[]);
    p.addParameter('ColorBar',[]);
    p.addParameter('TickFontSize',[]);
    p.addParameter('TickTexInterpreter',false);
    p.addParameter('Title','');
    p.addParameter('xTitle','');
    p.addParameter('yTitle','');
    p.addParameter('labelsx',[]);
    p.addParameter('labelsy',[]);
    p.addParameter('autoLabel',false);
    p.addParameter('xblocksize',1);
    p.addParameter('yblocksize',1);
else
    p.addParamValue('TickAngle',90);
    p.addParamValue('FontSize',[]);
    p.addParamValue('showVals',false);
    p.addParamValue('ColorBar',false);
    p.addParamValue('TickFontSize',[]);
    p.addParamValue('TickTexInterpreter',false);
    p.addParamValue('Title','');
    p.addParamValue('xTitle','');
    p.addParamValue('yTitle','');
    p.addParamValue('labelsx',[]);
    p.addParamValue('labelsy',[]);
    p.addParameter('autoLabel',false);
    p.addParameter('xblocksize',1);
    p.addParameter('yblocksize',1);
end
p.parse(varargin{:});

param = p.Results;

% Add a few other parameters
param.Colormap = 'ddc';
 
end

% Function to format data matrix into cell array of strings
function txtX = mat2cstr(X,format)
    
    txtX = cell(size(X));
    for j=1:size(X,2)
        Xsel = X(~isnan(X(:,j)),j);
        if all(Xsel==floor(Xsel))   %integer values
            txtX(:,j) = sprintfc('%i',X(:,j));  
        else                            %decimal values
            txtX(:,j) = sprintfc(format,X(:,j));
        end
    end
    
end

function Xgrad = squeeze(Xin ,n, d, xblocksize, yblocksize)
% function for use in cellMap By Block
X = zeros(n,d);
Xgrad = X;
for i = 1:n
    for j = 1:d
        Xsel = Xin((1+((i-1)*yblocksize)):(i*yblocksize),...
                   (1+((j-1)*xblocksize)):(j*xblocksize));
        seltable = sum(hist(Xsel,0:3),2);
        seltable = seltable(2:end);
        if (sum(seltable) > 0)
          [cntmax,indmax] = max(seltable);
          gradmax = cntmax / (xblocksize*yblocksize);
        else
          indmax = 0;
          gradmax = 1;
        end
        X(i,j) = indmax;
        Xgrad(i,j) = indmax+gradmax;
        if gradmax==1; Xgrad(i,j) = Xgrad(i,j)-0.0001;end
    end
end

end