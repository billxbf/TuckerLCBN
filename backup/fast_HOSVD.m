% Suna Guo
% 05.08.2019

% Fast Tucker decomposition for analyzing brain imagining data
% (dynamic correlation with sliding window method).
% Exploits special structure of correlaiton tensor, construct T(4)*T(4)'
% for participant mode HOSVD directly from the activity tensor
% without explicitly constructing the correlation tensor.
% Returns the left singular vectors of T(4).

% Includes algorithm for T(3) commented out at the bottom.

% *Note: 
% Tradeoff between running time and storage: 
% Center the windows right before putting the result in T(4)T(4)',
% or start from a tensor with all centered windows.

% Input: 
% A: activity tensor (area x time x subject)
% W: window size

% Output:
% U: Left singular vectors of T(4) (participant x (area x area x time))
%    from eig(T(4)*T(4)').

function [U,S,etime] = fast_HOSVD(A, W)
    tic
    [~, T, P] = size(A);    % a: #areas; T: total time points; N: #subj
    N = T-W+1;    
    
    % ***** G = T(4)T(4)' for eig decomp *****
    G = zeros(P, P);
    
    % Construct T(4)T(4)'
    for n = 1:N
        for i = 1:P
            
            % Shift & scale windows (need to do this more efficiently)
            Ai = A(:,n:n+W-1,i);
            Ai = Ai - mean(Ai,2);
            rnrms = vecnorm(Ai')';
            Ai = bsxfun(@rdivide,Ai,rnrms);
             
            for j = 1:i
                
                % Each entry is sum of F-norm^2 of participant i&j at window n
                
                % Shift & scale windows (need to do this more efficiently)
                Aj = A(:,n:n+W-1,j);
                Aj = Aj - mean(Aj,2);
                rnrms = vecnorm(Aj')';
                Aj = bsxfun(@rdivide,Aj,rnrms);
                
                % Put in T(4)T(4)'
                G(i,j) = G(i,j) + norm(Ai'*Aj, 'fro')^2;
                
                % Without centering window, for speedup
                %                 G(i,j) = G(i,j) + norm(A(:,t:t+W-1,i)'*A(:,t:t+W-1,j), 'fro')^2;
            end
        end
    end
    
    G = G + tril(G, -1)';
    [UG,EG] = eig(G);
    
    % put the eigenvalues in descending order
    [S,ind] = sort(diag(EG),'desc');
    UG = UG(:,ind);
    
    % Assign U with correct mat
    U = UG;
    etime = toc;
    %save('fast_HOSVDMat')
end

    
%     % ***** C = T(3)T(3)' for eig decomp *****
%     C = zeros(N, N);
%     
%     % Construct T(3)T(3)'
%     for k = 1:N   % For each participant to sum together
%         for i = 1:N
%             for j = 1:N
%                 C(i,j) = C(i,j) + norm(A(:,i:i+W-1,k)'*A(:,j:j+W-1,k), 'fro');
%             end
%         end
%     end
%     [UC,~] = eig(C);


%%%%%%%%%%%%%%%%%%