function [R] = sparse_cholesky(A)
n = length(A);              
R = A;
for k = 1:n-1
    R(k,k+1:n)     = 0;                                   % Setting all values of of k row in front of k column to 0
    R(k+1:n,k+1:n) = R(k+1:n,k+1:n) - ...
                    (R(k+1:n,k)*(R(k+1:n,k))')/(R(k,k));  % Computing submatrix
    R(k,k)         = (R(k,k))^0.5;                        % Reassigning entry k,k to its sq root
    R(k+1:n,k)     = R(k+1:n,k)./R(k,k);                  % Building v vector divided by sq root of first entry
end
R(end,end) = (R(end,end))^0.5;                            % Taking square root of last entry
R = R';    