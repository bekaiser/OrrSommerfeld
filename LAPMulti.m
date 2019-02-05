function [ C ] = LAPMulti( A,B)
% Linear Algebraic Product (LAP) of two multidimensional matrics of equal 
% dimensions
% A = Matrix of matrices (all square)
% B = Matrix of matrices (all square) OR(!) vector of vectors

Asize = size(A);
nA = Asize(1); % Minor rows in A
mA = nA; % Square minor matrix, minor cols
NA = Asize(3); % Major rows in A
MA = NA; % Square major matrix, major cols

% A2 = (A);
% A = (A2);
% B2 = (B);
% B = (B2);


if NA ~= MA
    fprintf('A is not square')
end

Bsize = size(B);
Temp = length(Bsize);

if Temp == 3 % B is a vector  
    %fprintf('B is a vector') 
    nB = Bsize(1); % Minor rows in B
    mB = Bsize(2); % Minor cols in B (should be 1)
    if mB ~= 1
        fprintf('Problem with vector dimensions, check it')
    end
    NB = Bsize(3); % Major rows in B (must equal cols of A)
    MB = 1; % One major column in B (must be 1)

    % Insert matrix dimension check: NB = MA 
    N = NB;
    Ci = zeros(nB,mB,N,1);
    for j = 1:N % For each row
        k = 1; % For each column
            E = zeros(nB,mB,1,N);
            Itr = j;
            Itc = k;
            for i = 1:N
                E(:,:,1,i) = (A(:,:,j,i)*B(:,:,i,k));
            end
            Ci(:,:,j,k) = (SumMultiArray(E));
        
    end
    C = (Ci);
    
    

elseif Temp == 4 % B is a matrix
    %fprintf('B is a matrix')    
    nB = Bsize(1); % Minor rows in B
    mB = Bsize(2); % Minor cols in B
    NB = Bsize(3); % Major rows in B
    MB = Bsize(4); % Major cols in B

    % Insert matrix dimension check: NB = MB = NA = MB
    N = MB;
    Ci = zeros(nB,mB,N,N);
    for j = 1:N % For each row
        for k = 1:N % For each column
            E = zeros(nB,mB,1,N);
            Itr = j;
            Itc = k;
            for i = 1:N
                E(:,:,1,i) = (A(:,:,j,i)*B(:,:,i,k));
            end
            Ci(:,:,j,k) = (SumMultiArray(E));
        end
    end
    C = (Ci);


else  
    fprintf('Error in LAPMatMat due to column sizing, check it')
end



end