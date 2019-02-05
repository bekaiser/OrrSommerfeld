function [ dy G1 ] = StructGrid( n )
% Generates a grid of Gauss Lobatto distribution

% N = Number of grid points
N = n;

% Generating the Gauss-Lobatto grid points
G1 = linspace(0,1,n); % Location of y, starting at 0
% and ending at 1.

% Vector dy (center points = N - 1)
dy = abs(diff(G1)); % y difference from cell bottom to cell top


end

