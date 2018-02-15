function V = calc_potHA3(phi,r)
% solves AU = b
%
%

h = r(2) - r(1); % Stepping distance 
N=length(r);
n = abs(phi).^2; % Electron density
UN = 0; % BC at r_max

% A
A = zeros(N,N);
for i = 2:N-1
   A(i,i) = -2;
   A(i,i-1) = 1;
   A(i,i+1) = 1;   
end
A(1,1) = -2;
A(N,N) = -2;
A(N,N-1) = 1;
A(1,2) = 1;

A = A/h^2;

% b
b = zeros(N,1);
for i = 1:N
   b(i) = -4*pi*r(i)*n(i);
end
b(N) = b(N) - UN/h^2;
% Solve system of equations
U_0 = A\b;
U = U_0 + r/max(r);
V = U./r;   % V_sH
end