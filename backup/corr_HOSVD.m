function [U,T4,etime] = corr_HOSVD(A, W)
    tic
    [a, T, P] = size(A);
    N = T-W+1;
    
    % Correlation tensor: a x a x #window x subj
    C = zeros(a,a,N, P);
    for j = 1:P
        for i = 1:N
            
            Aij = A(:,i:i+W-1,j);
            
            % Shift & scale window
            Aij = Aij - mean(Aij,2);
            Aij = bsxfun(@rdivide,Aij,vecnorm(Aij')');
            
            C(:,:,i,j) = Aij * Aij';
            
            % Without centering window, for speedup
%             C(:,:,i,j) = A(:,i:i+W-1,j) * A(:,i:i+W-1,j)';
        end
    end

    % T(4)
    T4 = reshape(permute(C, [4,1,2,3]), [P, a*a*N]);
    
    % HOSVD
    [U,E] = eig(T4 * T4');
    
    % put the eigenvalues in descending order
    [E,ind] = sort(diag(E));
    U = U(:,ind);
    etime = toc;
    %save('corr_HOSVDMat')
end