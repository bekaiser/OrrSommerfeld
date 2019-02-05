function [ C ] = SumMultiArray( A )
% Sum of Multidimensional Array across rows
% Input "A" Must be a single major row of a multidimensional array
% i.e. input: A(1:n,1:m,1,1:N)

% N = Ni; % major dimension
% M = 1;


S = size(A);
n = S(1);
m = S(2);
M = S(3); % Should equal 1
N = S(4);

A2 = (A);
A = (A2);

if M ~= 1
    fprintf('Error in SumMultiArray function, check it')
end

Cii = zeros(n,m);
for k = 1:n
    for j = 1:m

    Ci = zeros(1,N);
        for i = 1:N
            Ci(i) = (A(k,j,1,i)); % A(loc,loc,1,vector to be summed)
        end
        
    Cii(k,j) = (sum(Ci));
    end
end
C = (Cii);
% Output: C(1:n,1:m,1,1)

end

