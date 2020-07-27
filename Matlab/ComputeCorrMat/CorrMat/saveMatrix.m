function X = saveMatrix()
%Generate Random Matrix X
%rng(1);
%X = randn(10,4);

%Readin
T = readtable('timeSeries.csv');
X = table2array(T);
X = X';
%Save to binary file
fileID = fopen('matrix.bin','w');
fwrite(fileID,X,'double');
fclose(fileID);
end