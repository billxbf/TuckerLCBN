function X = saveMatrix()
%Generate Random Matrix X
rng(1);
X = randn(100,100);

%Save to binary file
fileID = fopen('matrix.bin','w');
fwrite(fileID,X,'double');
fclose(fileID);
end