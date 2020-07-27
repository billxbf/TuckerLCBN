A = extractData('regionBased', 268);

%Preprocess by removing TxP slices where there
%exits a all-0 fiber along T.
[a, T, P] = size(A);
delList = zeros(1,a);
for i = 1:a
    for j = 1:P
        if sum(A(i,:,j)) == 0
            delList(i) = 1;
            break 
        end
    end
end
A(logical(delList),:,:) = [];

W = 60;
[U,T4, etime1] = corr_HOSVD(A,W);