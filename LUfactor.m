function [ L U ] = LUfactor( A )
% LU factorization
% A must be a matrix of matrices

Temp = size(A);
n = Temp(1); % Number of scalar points in sub-matrices
N = Temp(3); % Number of sub-matrices in A matrix


A2 = (A);
A = (A2);

% Matrices of split tridiagonal
T1 = eye(n,n);
L = zeros(n,n,N,N);
for i = 1:N
    L(:,:,i,i) = (T1); % L matrix set up, with I matrices
end
L(:,:,2,1) = (A(:,:,2,1)/A(:,:,1,1)); % T1 = B1/A0 = B1/D0


U = zeros(n,n,N,N); % U matrix set up, all zeros
U(:,:,1,1) = (A(:,:,1,1)); % Delta(0,0) = A0
for i = 1:N-1
    U(:,:,i,i+1) = (A(:,:,i,i+1)); % C0 through CJ
end


for i = 2:N-1
    U(:,:,i,i)= (A(:,:,i,i)) - (L(:,:,i,i-1))*(A(:,:,i-1,i)); % Filling out D values, 1-(J-1), starting at the second row L(:,:,i,i-1),Z   A(:,:,i,i-1)/A(:,:,i-1,i-1),Z
    L(:,:,i+1,i)= (A(:,:,i+1,i)/U(:,:,i,i)); % Filling out T values, T2-TJ started at the third row    
end
%whos L

U(:,:,N,N) = A(:,:,N,N)-L(:,:,N,N-1)*A(:,:,N-1,N); % DJ



end
