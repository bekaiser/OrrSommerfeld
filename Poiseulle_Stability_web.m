% Numerical Analysis of Stability in a Plane Poiseulle Flow
% Orr Sommerfeld Equation
% 05/03/13
% Bryan Kaiser

clear all
clc
close all


% Neutral Stability Curve Loop set up

% 1) User Specified values -------------------------------
I = 15; % Number of iterations
% Vectors of important variables per iteration:
ar = zeros(1,I); % alpha real
dar = zeros(1,I);
wr = zeros(1,I); % omega real
dwr = zeros(1,I);
T1 = zeros(1,I);
T2 = zeros(1,I);
Numa = zeros(1,I);
Numw = zeros(1,I);
ar(1) = 2.040697; % alpha = real component of the small disturbance wavenumber
wr(1) = 0.80806; % omega = real component of the small disturbance frequency
Re = 7695.559 % Reynolds number based on channel height h
N = 2500; % N = gridpoints
e1 = 11; % Test 10)
e2 = 15; % Test 14)

% 2a) Grid generation ---------------------------------------------------------
[ dy G ] = GLgrid(N); % Function that creates a Gauss-Lobatto grid 
% G = locations of all y's, starting at 0 and ending at 1.
% dy = y difference from cell bottom to cell top

% 2b) Grid Plot Option
% plot(1:(N-1),dy)
% xlabel('Gridpoint')
% ylabel('dy')
% g = fliplr(G)
% plot(1:(N),g)
% xlabel('Gridpoint')
% ylabel('height')

% 3) Plane Poiseulle channel flow velocity U(y) profile ----------------------
% U(y) = 1.5 - 6(y-1/2)^2 ,  U(y)'' = -12 
Uy = zeros(1,N-1); % Set up
for i = 1:N-1
    Uy(i) = (1.5-6*( (G(i+1)+G(i))/2-1/2)^2);
end
Uydp = -12; % U(y)''

tic
% 4) Beginning to iterate alpha & omega
for j = 1:(I-1) % Number of iterations
     
% 5) Equation matrix "A" set up ----------------------------------------------
n = 4;
A = zeros(n,n,N,N); % A Matrix populated with zeros only.

% A Matrices
A(:,:,1,1) = [ 1,0,0,0; 0,1,0,0; -2,0,-dy(1), 0; 0,-2,0,-dy(1) ]; % BC
for i = 1:N-2
    A(:,:,i+1,i+1) = [ -(ar(j))^2*dy(i), -dy(i), 2, 0; Re*ar(j)*Uydp*dy(i)*sqrt(-1), (-dy(i))*((ar(j)^2)+(sqrt(-1)).*Re*(ar(j)*Uy(i)-wr(j))) , 0, 2; -2, 0, -dy(i+1), 0;  0, -2, 0, -dy(i+1) ];
end
A(:,:,N,N) = [ -(ar(j))^2*dy(N-1), -dy(N-1), 2, 0;  (Re*ar(j)*Uydp*dy(N-1)).*(sqrt(-1)), (-dy(N-1))*( (ar(j))^2+(sqrt(-1))*Re*(ar(j)*Uy(N-1)-wr(j) )  ), 0, 2; 0, 0, 1, 0;  1, 0, 0, 0 ];

% B Matrices
for i = 1:N-1
    A(:,:,i+1,i) = [ -(ar(j))^2*dy(i), -dy(i), -2, 0; Re*ar(j)*Uydp*dy(i)*sqrt(-1), (-dy(i))*((ar(j)^2)+(sqrt(-1)).*Re*(ar(j)*Uy(i)-wr(j))),0,-2;  0,0,0,0;  0,0,0,0 ];
end
A(:,:,N,N-1) = [ -(ar(j))^2*dy(N-1), -dy(N-1), -2, 0; Re*ar(j)*Uydp*dy(N-1)*(sqrt(-1)), (-dy(N-1))*((ar(j)^2)+(sqrt(-1)).*Re*(ar(j)*Uy(N-1)-wr(j))),0,-2;  0,0,0,0;  0,0,0,0 ];

% Inner C Matrices
for i = 1:(N-1)
    A(:,:,i,i+1) = ([ 0,0,0,0; 0,0,0,0; 2,0,-dy(i),0; 0,2,0,-dy(i)]);
end 

% 6) R Matrix Setup -------------------------------------------------------
R = zeros(4,1,N,1);
R(2,1,1,1) = 1; % The top of R is y = 0, where s = 1, and  the bottom is y = h.
Rs = R; % Saved value of original R, for later in the loop.

% 7) L U Factorization of A matrix ----------------------------------------
[ L U ] = LUfactor(A); 

% 8) Putting L & U into vector form ---------------------------------------
S = size(L);
n = S(1);
N = S(3);
T = zeros(n,n,1,N-1);
for i = 1:N-1
    T(:,:,1,i) = (L(:,:,i+1,i)); % Vector of Lower diagonal of L
end
Ts = (T); % Saved value of original T.
D = zeros(n,n,1,N);
for i = 1:N
    D(:,:,1,i) = (U(:,:,i,i)); % Vector of Lower diagonal of U
end
Ds = (D); % Saved value of original D.
C = zeros(n,n,1,N-1);
for i = 1:N-1
    C(:,:,1,i) = (U(:,:,i,i+1)); % Vector of Upper diagonal of U
end
Cs = (C); % Saved value of original C.

% 9) Solving for delta (d) unknown column of vectors ---------------------
% Creating the omega (w) vector of secondary unknown vectors 
% Downwards sweep:
w = zeros(4,1,N,1);
w(2,1,1,1) = 1; % w0 = R0 (Notes 2/14)
for i = 1:N-1
    w(:,1,i+1,1) = (R(:,1,i+1,:))-(T(:,:,1,i)*w(:,1,i,:));
end
% Creating the delta (d) vector of primary unknown vectors 
% Upwards sweep:
Li = [1:N];
L = fliplr(Li);
d = zeros(4,1,N,1);
d(:,1,N,1) = ((D(:,:,1,N))\w(:,1,N,1)); % dJ = (DJ^(-1))*wJ (Notes 2/14)
for i = 1:N-1
    d(:,1,L(i+1),:) = (D(:,:,1,L(i+1))\w(:,1,L(i+1),:))-(D(:,:,1,L(i+1))\(C(:,:,1,L(i+1))*d(:,1,L(i),:))); 
end
% The above loop: dj = (Dj^(-1))*(wj - cj*Dj+1), (Notes 2/14)
% The target variable f(y=h) => 0 :
T1(j) = (real(d(3,1,1))); % fr(y=h) real part of missing BC
T2(j) = (imag(d(3,1,1))); % fi(y=h) imag part of missing BC

% 10) Check to be inserted: A*d = R ? -------------------------------------
%[ Rcheck ] = LAPMulti( A,d );
% IsAns2 = find(Rcheck >= 1*10^(-e1))

% 11) Partial Equation Matrices (dA/da and dA/dw)--------------------------
% Equation matrix dAda set up
% Partial derivative of A with respect to alpha = dAda
n = 4;
dAda = zeros(n,n,N,N); % dAda Matrix populated with zeros only. 
% with respect to alpha, in vector form (see line 43).
% Inner B Matrices of dAda
for i = 1:(N-1)
    dAda(:,:,i+1,i) = ([ -2*(ar(j))*dy(i),0,0,0; Re*Uydp*dy(i).*(sqrt(-1)), -dy(i)*(2*(ar(j))+sqrt(-1).*Re*(Uy(i))) ,0,0; 0,0,0,0; 0,0,0,0]);
end
% Inner A Matrices of dAda
for i = 1:(N-1)
    dAda(:,:,i+1,i+1) = ([ -2*(ar(j))*dy(i),0,0,0; sqrt(-1).*Re*Uydp*dy(i), -dy(i)*(2*(ar(j))+sqrt(-1).*Re*(Uy(i))) ,0,0; 0,0,0,0; 0,0,0,0]);
end

% Equation matrix dAdw set up
% Partial derivative of A with respect to omega = dAdw
n = 4;
dAdw = zeros(n,n,N,N); % dAdw Matrix populated with zeros only.
% Inner B Matrices of dAda
for i = 1:(N-1)
    dAdw(:,:,i+1,i) = ([0,0,0,0; 0,dy(i)*sqrt(-1)*Re,0,0; 0,0,0,0; 0,0,0,0]);
end
% Inner A Matrices of dAda
for i = 1:(N-1)
    dAdw(:,:,i+1,i+1) = ([0,0,0,0; 0,dy(i)*sqrt(-1).*Re,0,0; 0,0,0,0; 0,0,0,0]);
end

% 12) New R vector column (for dAda and dAdw)------------------------------
% Since A*(dd/dw) = -(dA/dw)*d and A*(dd/da) = -(dA/da)*d:
% Ra = -(dA/da)*d 
% Rw = -(dA/dw)*d 
% and A*(dd/da) = Ra and A*(dd/dw) = Rw.
Ra = LAPMulti( -dAda, d ); 
Rw = LAPMulti( -dAdw, d ); 

% 13) Solving for delta (d) unknown column of vectors for Ra & Rw  --------
% Creating the omega (w) vector of secondary unknown vectors 
T = (Ts);
D = (Ds); % Checked, saved values work properly.
C = (Cs);
% Downwards sweep (matrix) up physically
w = zeros(4,1,N,1);
w(:,:,1,1) = (Ra(:,:,1,1)); % W0 = Ra0
for i = 1:N-1
    w(:,1,i+1,1) = (Ra(:,1,i+1,1))-(T(:,:,1,i)*w(:,1,i,:));
end
% Creating the delta (d) vector of primary unknown vectors 
% Upwards sweep
Li = [1:N];
L = fliplr(Li);
ddda = zeros(4,1,N,1);
ddda(:,1,N,1) = (D(:,:,1,N)\w(:,1,N,1)); % Final d value 
for i = 1:N-1
    ddda(:,1,L(i+1),:) = (D(:,:,1,L(i+1))\(w(:,1,L(i+1),:)))-(D(:,:,1,L(i+1))\(C(:,:,1,L(i+1))*ddda(:,1,L(i),:)));
end
% Creating the omega (w) vector of secondary unknown vectors 
% Downwards sweep
w = zeros(4,1,N,1);
w(:,:,1,1) = Rw(:,:,1,1); % W0 = Rw0
for i = 1:N-1
    w(:,1,i+1,1) = (Rw(:,1,i+1,:))-(T(:,:,1,i)*w(:,1,i,:));
end
% Creating the delta (d) vector of primary unknown vectors 
% Upwards sweep
Li = [1:N];
L = fliplr(Li);
dddw = zeros(4,1,N,1);
dddw(:,1,N,1) = (D(:,:,1,N)\w(:,1,N,1)); % Final d value 
for i = 1:N-1
    dddw(:,1,L(i+1),:) = (D(:,:,1,L(i+1)) \w(:,1,L(i+1),:))-(D(:,:,1,L(i+1))\(C(:,:,1,L(i+1))*dddw(:,1,L(i),:)));
end

% 14) Check to be inserted: A*dd/dw = -dAdw*d and A*dd/da = Ra ? <--------------
% Rwtest = vpa(LAPMulti( A, dddw ));
% RwAns = find(Rw >= 1*10^(-e2));
% IsRwtestsame = find(Rwtest >= 1*10^(-e2));
% Before = Rw(:,:,N-1,1);
% After = Rwtest(:,:,N-1,1);
% RSsum = sum(Rw(RwAns(:)));
% RSsumtest = sum(Rwtest(IsRwtestsame(:)));
% Shouldbeone = RSsum/RSsumtest
% Very close for N = 20!

%%% CODE VERIFIED CORRECT UP TO THIS POINT %%%

% 15) Target variables at df/dw(y=0)and df/da(y=0) ------------------------
dT1da = (real(ddda(3,:,1))); % dfr/da(y=0)
dT2da = (imag(ddda(3,:,1))); % dfi/da(y=0)
dT1dw = (real(dddw(3,:,1))); % dfr/dw(y=0)
dT2dw = (imag(dddw(3,:,1))); % dfi/dw(y=0)

% 16) Calulation of next guess for alpha and omega ------------------------
Den = (dT1da*dT2dw-dT2da*dT1dw);
Numa(j) = (T2(j)*dT1dw-T1(j)*dT2dw); 
dar(j) = (Numa(j)/Den); % Change in alpha for next iteration
dara = abs(dar);
Numw(j) = (T1(j)*dT2da-T2(j)*dT1da); 
dwr(j) = (Numw(j)/Den); % Change in omega for next iteration
dwra = abs(dwr);
if dara(j) <= 10^(-8)
    if dwr(j) <= 10^(-8)
        formatSpec = 'At iteration %2.0f, convergence requirement met \n';
        fprintf(formatSpec,j)
    end
end

% 17) The next Newton Iteration:
ar(j+1) = (ar(j) + dar(j));
wr(j+1) = (wr(j) + dwr(j));
%pause

end
toc

% 18) Alpha Results
dara = abs(dar);
dwra = abs(dwr);
Convaloc = find(dara <= 10^(-8)); % Convergence requirement
Convwloc = find(dwra <= 10^(-8)); % Convergence requirement
Conva = dar(Convaloc(1)); % Converged requirement value
Convw = dwr(Convwloc(1)); % Converged requirement value

formatSpec = 'For iteration %2.0f and Re = %2.0f, alpha %4.3f and omega %4.3f \n';
fprintf(formatSpec,I(end),Re, ar(end),wr(end))
formatSpec = 'At iteration %2.0f, alpha convergence and omega convergence below 10^(-8) \n';
fprintf(formatSpec,Convaloc(1))

format Long
Result_Vector = [ 20; ar(end); wr(end); Convaloc(1) ]

% 19) Iteration Plots
figure
set(gca,'FontName','Times New Roman','Fontsize',16)
p1 = plot(1:I,dar);
xlabel('Iteration')
ylabel('d alpha')
figure
set(gca,'FontName','Times New Roman','Fontsize',16)
p2 = plot(1:I,dwr);
xlabel('Iteration')
ylabel('d omega')


% 20) Eigenvalue Plots 
% eta = G; % Dell's notation
% del = d; % Dell's notation
% for j=1:N
%     phi(j)=del(1,1,j,1);
%     s(j)=del(2,1,j,1);
%     f(j)=del(3,1,j,1);
%     g(j)=del(4,1,j,1);
% end
% subplot(1,4,1)
% set(gca,'FontName','Times New Roman','LineWidth',2,'Fontsize',14)
% p = plot(real(phi),eta,'r',imag(phi),eta,'c',abs(phi),eta);
% set(p,'LineWidth',2)
% grid on;
% set(gca,'FontName','Times New Roman','LineWidth',2,'Fontsize',14)
% xlabel('\phi');
% ylabel('y');
% subplot(1,4,2)
% p1 = plot(real(s),eta,'r',imag(s),eta,'c',abs(s),eta);
% set(p1,'LineWidth',2)
% grid on;
% set(gca,'FontName','Times New Roman','LineWidth',2,'Fontsize',14)
% xlabel('s');
% ylabel('y')
% subplot(1,4,3)
% p2 = plot(real(f),eta,'r',imag(f),eta,'c',abs(f),eta);
% set(p2,'LineWidth',2)
% grid on;
% set(gca,'FontName','Times New Roman','LineWidth',2,'Fontsize',14)
% xlabel('f');
% ylabel('y');
% subplot(1,4,4)
% p3 = plot(real(g),eta,'r',imag(g),eta,'c',abs(g),eta);
% set(p3,'LineWidth',2)
% grid on;
% set(gca,'FontName','Times New Roman','LineWidth',2,'Fontsize',14)
% xlabel('g');
% ylabel('y');
% legend('Real','Imaginary','Mag.','Location','NorthEastOutside');
