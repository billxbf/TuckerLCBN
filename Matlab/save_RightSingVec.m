function save_RightSingVec(A, W, U)
    [a, T, P] = size(A);
    N = T-W+1;
    U1 = U(:,1)';
    U2 = U(:,2)';
    % Correlation tensor: a x a x #window x subj
    for i = 1:N
        V1 = zeros(a);
        V2 = zeros(a);
        for j = 1:P            
            Aij = A(:,i:i+W-1,j);
            % Shift & scale window
            Aij = Aij - mean(Aij,2);
            Aij = bsxfun(@rdivide,Aij,vecnorm(Aij')');           
            C = Aij * Aij';
            V1 = V1 + C * U1(j);
            V2 = V2 + C * U2(j);
        end
        fileID = fopen('TMP/V1window'+string(i)+'.bin','w');
        fwrite(fileID,V1,'double');
        fclose(fileID);
        fileID = fopen('TMP/V2window'+string(i)+'.bin','w');
        fwrite(fileID,V2,'double');
        fclose(fileID);

    end

end