function D = DistanceMatrix_NP2(A,W)
    [a, T, P] = size(A);
    N = T-W+1;
    D = zeros(N*P);
    for i = 1:N
        for j = 1:i
            [i,j]
            for a = 1:P
                if i == j
                    blim = a;
                else
                    blim = P;
                end
                
                W1 = A(:,i:i+W-1, a);
                W1W1 = W1'*W1; % Store W1'W1 to avoid redundancy ?
                W1W2s = W1'*reshape(A(:,j:j+W-1,1:blim), [size(A(:,j:j+W-1,1:blim),1), W*blim]); 
                for b = 1:blim
                    W2W2 =  A(:,i:i+W-1, b)' * A(:,i:i+W-1, b);
                    D((i-1)*P+a, (j-1)*P+b) = sqrt(norm(W1W1,'fro')^2 - 2*norm(W1W2s(:,(b-1)*W+1:b*W), 'fro')^2 + norm(W2W2,'fro')^2);
                                     
                end
            end
        end
    end
    


end