function [ dy G1 ] = GLgrid( n )
% Generates a grid of Gauss Lobatto distribution

% N = Number of grid points
N = n;

% Generating the Gauss-Lobatto grid points
G1i(1:N) = (0.5+cos((0:N-1)*pi/(N-1)) /2); % Location of y, starting at 0
% and ending at 1.

G1 = [ 0 fliplr(G1i(1:N-1)) ];

% Vector dy (center points = N - 1)
dy = abs(diff(G1)); % y difference from cell bottom to cell top

% Test plot
% close all
% Y1 = linspace(0,1,N-1);
% plot(dy,Y1)
% title('Gauss Lobatto Grid')
% xlabel('Cell Center Height (dy)')
% ylabel('Channel Height')

end

