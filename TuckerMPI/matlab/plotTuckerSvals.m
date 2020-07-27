function plotTuckerSvals(datadir, prefix, nd, myeps, figid)
%PLOTTUCKERSVALS Plot TuckerMPI mode-n singular values to predict cutoffs.
%
%   PLOTTUCKEREVALS(DIR,PREFIX,NDIMS,EPSILON) displays the singular values
%   from the ST-HOSVD code in the TuckerMPI library in such a way that it
%   can be used to predict the cutoff values for different choices of
%   EPSILON. DIR Is the data where the singular value files live, PREFIX is
%   the file prefix (i.e., each file is named <PREFIX>_mode_X.txt), and
%   NDIMS is the number of dimensions. We assume the singular values for
%   each mode are ordered greatest to least, i.e., 
%
%     lambda_1 >= lambda_2 >= ... >= lambda_{N_k}
%
%   where N_k is the size of mode k. The plot shows the cummulative sum,
%   i.e., (i, sum(lambda(i:end)). It also plots the cutoff based on
%   EPSILON, which is ||X||^2 EPSILON^2 / NDIMS. The value of ||X||^2 is
%   the maximum sum of the singular values across all modes.
%
%   PLOTTUCKEREVALS(DIR,PREFIX,NDIMS,EPSILON,FIGID) also specifies the
%   figure to plot. If none is specified, then a new figure is used.
%
%   Note that because the singular values are computed after truncation,
%   the compression ratios may be somewhat speculative but do provide an
%   accurate bound for any EPSILON less than or equal to the values used in
%   the ST-HOSVD computation.
%
%   Examples:
%   plotTuckerSvals('C:\Users\tgkolda\Documents\tmp','sv',3,0.1,1)

%% Check inputs
if nargin < 4
    error('Need 4 input arguments');
end
   
%% Read the data
svals = cell(nd,1);
cumsevals = cell(nd,1);
lgnd = cell(nd,1);
for n = 1:nd
    fname = sprintf('%s_mode_%d.txt',prefix,n-1);
    svals{n} = importdata(fullfile(datadir,fname));
    cumsevals{n} = cumsum(svals{n},'reverse');
    lgnd{n} = sprintf('Mode %d',n);
end

%% Assume max sum is the tensor norm (not sure order that modes were processed)

% Extract first elements
Xnormsqr = max(cellfun(@(x) x(1), cumsevals));

% Get tensor size
tsz = cellfun(@length, cumsevals);

% Compute sum cutoff
cutoff = Xnormsqr * myeps^2 / nd;

% Determine cutoff
rsz = cellfun(@(x) find(x >= cutoff,1,'last'), cumsevals);


%% Plot 
if ~exist('figid','var')
    figure;
else
    figure(figid);
end
clf;

for n = 1:nd
    semilogy(cumsevals{n})
    hold on;
end
xl = xlim;
plot(xl, [cutoff cutoff], 'k--');
lgnd{end+1} = 'Cutoff';
legend(lgnd);
reduction = prod(tsz)/prod(rsz);
gtitle = sprintf('\\epsilon = %.1e: %s to %s (%.0fX)', ...
    myeps, size2str(tsz), size2str(rsz), reduction);
title(gtitle);

function s = size2str(sz)
%SIZE2STR Convert size to a string that can be printed.

if isempty(sz)
    s = sprintf('[empty]');
    return;
end

if numel(sz) == 1
    s = sprintf('%d',sz);
else
    s = [sprintf('%d x ',sz(1:end-1)) sprintf('%d', sz(end)) ];
end


