function processTuckerTimers(csv_filename, nmodes,fignum)
%processTuckerTimers Plots the runtime of various segments of the Tucker
%code

% Read a CSV
A = importdata(csv_filename);

% Strip out commands that never got executed
i=1;
while(i <= size(A.data,2))
    maxent = max(A.data(:,i));
    if(maxent <= 0)
        A.colheaders(i) = [];
        A.data(:,i) = [];
    else
        i = i+1;
    end
end

% Plot each mode separately
runtimeSoFar=zeros(size(A.data,1),1);
for i=1:nmodes
    % Find all entries that contain this mode number
    fig = figure(fignum+i);
    modeSubs = arrayfun(@(x) strfind(x,int2str(i-1)), A.colheaders);
    modeSubs = arrayfun(@(x) ~cellfun('isempty',x), modeSubs, 'UniformOutput', false);
    modeSubs = find(cell2mat(modeSubs));
    nthings = size(modeSubs,2);
    
    % Determine which entry is Gram(i-1)
    gramSubs = arrayfun(@(x) strfind(x,strcat('Gram(',int2str(i-1))), A.colheaders);
    gramSubs = arrayfun(@(x) ~cellfun('isempty',x), gramSubs, 'UniformOutput', false);
    gramSubs = find(cell2mat(gramSubs));
    
    % Determine which entry is Eigensolve(i-1)
    eigSubs = arrayfun(@(x) strfind(x,strcat('Eigensolve(',int2str(i-1))), A.colheaders);
    eigSubs = arrayfun(@(x) ~cellfun('isempty',x), eigSubs, 'UniformOutput', false);
    eigSubs = find(cell2mat(eigSubs));
    
    % Determine which entry is TTM(i-1)
    ttmSubs = arrayfun(@(x) strfind(x,strcat('TTM(',int2str(i-1))), A.colheaders);
    ttmSubs = arrayfun(@(x) ~cellfun('isempty',x), ttmSubs, 'UniformOutput', false);
    ttmSubs = find(cell2mat(ttmSubs));
    
    colormap(fig,[[1 1 1]; jet])
    bar([runtimeSoFar A.data(:,[gramSubs+1:eigSubs ttmSubs+1:gramSubs+nthings-1])],'stacked')
    hold on
    plot(A.data(:,gramSubs)+runtimeSoFar,'k')
    plot(A.data(:,ttmSubs)+A.data(:,gramSubs)+A.data(:,eigSubs)+runtimeSoFar,'k')
    hold off
    axis tight
    legend(['total so far' A.colheaders([gramSubs+1:eigSubs ttmSubs+1:gramSubs+nthings-1 gramSubs ttmSubs])]);
    runtimeSoFar = runtimeSoFar + A.data(:,gramSubs) + A.data(:,eigSubs) + A.data(:,ttmSubs);
    
    runtimeSoFar = runtimeSoFar - min(runtimeSoFar);
end

% Only plot the higher-level routines
fig = figure(fignum+nmodes+1);
out = arrayfun(@(x) strfind(x,' '), A.colheaders);
out = arrayfun(@(x) cellfun('isempty',x), out, 'UniformOutput', false);
out = find(cell2mat(out));
A.data = A.data(:,out);
A.colheaders = A.colheaders(out);

out = arrayfun(@(x) strfind(x,'Eigensolve'), A.colheaders);
out = arrayfun(@(x) cellfun('isempty',x), out, 'UniformOutput', false);
out = find(cell2mat(out));
    
%colors = distinguishable_colors(size(out,2),'w',@(x) colorspace('RGB->Lab',x));
%set(groot,'defaultAxesColorOrder',colors)
nthings = size(out,2);
bar(A.data(:,out(1:nthings-1)),'stacked')
hold on
plot(A.data(:,out(nthings)),'k')
hold off
legend(A.colheaders(out))
axis tight
colormap(fig,jet)
end