function [U,S,V,etime] = corr_3D_HOSVD(A, W)
    tic
    [a, T, P] = size(A);
    N = T-W+1;
    
    % Correlation tensor: a x a x #window x subj
    C = zeros(a*(a-1)/2,N, P);
    mask = tril(true(a,a),-1);
    
    for j = 1:P
        for i = 1:N
            
            Aij = A(:,i:i+W-1,j);
            
            % Shift & scale window
            Aij = Aij - mean(Aij,2);
            Aij = bsxfun(@rdivide,Aij,vecnorm(Aij')');
            
            M = Aij * Aij';
            C(:,i,j) = M(mask);
            
            % Without centering window, for speedup
%             C(:,:,i,j) = A(:,i:i+W-1,j) * A(:,i:i+W-1,j)';
        end
    end

    % T(3)
    T3 = reshape(permute(C, [3,1,2]), [P, a*(a-1)/2*N]);
    
    % HOSVD
    [U,E] = eig(T3 * T3');
    
    % put the eigenvalues in descending order
    [S,ind] = sort(diag(E),'desc');
    U = U(:,ind);
    V = diag(S)\U'*T3;
    
    etime = toc;
    %save('corr_HOSVDMat')
end