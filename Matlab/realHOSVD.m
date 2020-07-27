vox = 1;

if vox
    a = extractData('../data/Original_Data/voxelBased1/1back', 18225);
    b = extractData('../data/Original_Data/voxelBased1/2back', 18225);
    c = extractData('../data/Original_Data/voxelBased1/rest', 18225);
else 
    a = extractData('../data/Original_Data/regionBased1/1back', 268);
    b = extractData('../data/Original_Data/regionBased1/2back', 268);
    c = extractData('../data/Original_Data/regionBased1/rest', 268);
    a(268,:,:) = [];
    a(131,:,:) = [];
    b(268,:,:) = [];v
    b(131,:,:) = [];
    c(268,:,:) = [];
    c(131,:,:) = [];
end

A = cat(3,a,b);
A = cat(3,A,c);

assert(nnz(A)>0);
tic;
W = 60;
[U,~,T4,etime1] = fast_HOSVD(A,W);
etime = toc;
%save_RightSingVec(A, W, U)
%writematrix(U1, 'voxelSingVec_full.csv')
%writematrix(S1, 'voxelSingVal_full.csv')
%writematrix(D,'voxelDistance_full.csv')
