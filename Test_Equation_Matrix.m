% Equation Matrix Set up

clear all
clc
close all


ar = 1.93; % alpha = real component of the small disturbance wavenumber
wr = 0.336; % omega = real component of the small disturbance frequency
Re = 50000; % Reynolds number based on channel height h - Does this affect 
N = 100; % N = gridpoints. This parameter determines the location of all 
n = 4;
Damp = 0.1;
j = 1;

[ dy G ] = GLgrid(N); % Function that creates a Gauss-Lobatto grid

Uy = zeros(1,N-1); % Set up
for i = 1:N-1
    Uy(i) = (1-4*( (G(i+1)-G(i))/2-1/2)^2);
end
Uydp = -8; % U(y)''

% length(Uy)
% length(dy)


A = zeros(n,n,N,N); % A Matrix populated with zeros only.
% 4 Sub rows, 4 Sub columns, N rows of matrices, N columns of matrices
% "A" matrices corresponds to a {phi,s,f,g} unknown vector
% Inner B Matrices

%%%%%%%%%%%%%%%%%%%%%

% A Matrices
A(:,:,1,1) = [ 1,0,0,0; 0,1,0,0; -2,0,-dy(1), 0; 0,-2,0,-dy(1) ]; % BC
for i = 1:N-2
    A(:,:,i+1,i+1) = [ -(ar(j))^2*dy(i), -dy(i), 2, 0; Re*ar(j)*Uydp*dy(i)*sqrt(-1), (-dy(i))*((ar(j)^2)+(sqrt(-1)).*Re*(ar(j)*Uy(i)-wr(j))) , 0, 2; -2, 0, -dy(i+1), 0;  0, -2, 0, -dy(i+1) ];
end
A(:,:,N,N) = [ -(ar(j))^2*dy(N-1), -dy(N-1), 2, 0;  (Re*ar(j)*Uydp*dy(N-1)).*(sqrt(-1)), (-dy(N-1))*( (ar(j))^2+(sqrt(-1))*Re*(ar(j)*Uy(N-1)-wr(j) )  ), 0, 2; 0, 0, 1, 0;  1, 0, 0, 0 ];

format long
A(:,:,1,1);
A(:,:,N-1,N-1);

% B Matrices
for i = 1:N-1
    A(:,:,i+1,i) = [ -(ar(j))^2*dy(i), -dy(i), -2, 0; Re*ar(j)*Uydp*dy(i)*sqrt(-1), (-dy(i))*((ar(j)^2)+(sqrt(-1)).*Re*(ar(j)*Uy(i)-wr(j))),0,-2;  0,0,0,0;  0,0,0,0 ];
end
A(:,:,N,N-1) = [ -(ar(j))^2*dy(N-1), -dy(N-1), -2, 0; Re*ar(j)*Uydp*dy(N-1)*(sqrt(-1)), (-dy(N-1))*((ar(j)^2)+(sqrt(-1)).*Re*(ar(j)*Uy(N-1)-wr(j))),0,-2;  0,0,0,0;  0,0,0,0 ];

format long
%A(:,:,1,1)
A(:,:,N-1,N-2)
A(:,:,N,N-1)

% Inner C Matrices
for i = 1:(N-1)
    A(:,:,i,i+1) = ([ 0,0,0,0; 0,0,0,0; 2,0,-dy(i),0; 0,2,0,-dy(i)]);
end 

format long
%A(:,:,1,1)
A(:,:,1,2);
A(:,:,N-1,N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
