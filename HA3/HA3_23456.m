%% 2
clc, close all
% r
r_min = 0;
r_max = 10;
N = 1000;
r = linspace(r_min,r_max,N+2);
r = r(2:end-1)';


phi=1/sqrt(pi)*exp(-r); % Hydrogen ground-state wave function

V = calc_potHA3(phi,r); % Calculate potential

plot(r,V,'cyan')
hold on
V_H = 1./r - (1+1./r).*exp(-2*r);
plot(r,V_H,'--r') 

xlabel('$r$ [a$_0$]','Interpreter','Latex','fontsize',16)
ylabel('Hydrogen potential [E$_H$] ','Interpreter','Latex','fontsize',16)
title('Hartree potential for H$_{1s}$','Interpreter','Latex','fontsize',16)

%% 3
clc, close all
% r
r_min = 0;
r_max = 10;
N = 1000;
r = linspace(r_min,r_max,N+2);
r = r(2:end-1)';
Z=1;    % Hydrogen

V=zeros(N); % No external potential
H=hamiltonianHA3(V,r,Z);

[F,lambda] = eig(H);
u = F(:,1); % Eigenfunction to radial Kohm-Sham
phi=u/sqrt(trapz(r,u.^2)); % Normalized wave function
E=lambda(1,1) % Ground state energy

plot(r,-phi,'blue')
xlabel('$r$ [a$_0$]','Interpreter','Latex','fontsize',16)
ylabel('$\phi_{1s}$ [a$_0^{-3/2}$] ','Interpreter','Latex','fontsize',16)
title('Wave function for H$_{1s}$','Interpreter','Latex','fontsize',16)

%% 4
clc, close all
% r
r_min = 0;
r_max = 10;
N = 1000;
r = linspace(r_min,r_max,N+2);
r = r(2:end-1)';

phi=ones(N,1)/N;   % Initial guess for wave function
E1 = 0; % Set energies to enter while-loop
E2 = 3;
tol=1e-7;   % Convergence tolerance
Z=2;    % Helium
i = 0;  % Number of iterations
% Self consistency loop
while abs(E2-E1) > tol
    E1 = E2;
    V = calc_potHA3(phi,r); % V_sH
    
    A = hamiltonianHA3(V,r,Z);    
    [F,lambda] = eig(A);
    
    u = F(:,1); % Eigenfunction to radial Kohm-Sham
    phi = u./(sqrt(4*pi)*r);      % Wave function  
    phi = phi / sqrt(trapz(r,4*pi*r.^2.*phi.^2));  % Normalization
    E2 = 2*lambda(1,1) - trapz(r,4*pi*r.^2.*V.*abs(phi).^2)  % Calculate energy  

    i = i+1;
end
i
E = E2

n = abs(phi).^2;
n=2*n/trapz(r,4*pi*r.^2.*n);
plot(r,-phi,'g')
xlim([0 2])
hold on
plot(r,-phi,'blue')
xlabel('$r$ [a$_0$]','Interpreter','Latex','fontsize',16)
ylabel('$\phi_{1s}$ [a$_0^{-3/2}$] ','Interpreter','Latex','fontsize',16)
title('Wave function for He$_{1s}$','Interpreter','Latex','fontsize',16)

%% 5
clc, %close all
% r
r_min = 0;
r_max = 10;
N = 2000;
r = linspace(r_min,r_max,N+2);
r = r(2:end-1)';

phi=ones(N,1)/N;   % Initial guess for wave function
E1 = 0; % Set energies to enter while-loop
E2 = 3;
tol=1e-5;   % Convergence tolerance
Z=2;    % Helium
i = 0;  % Number of iterations
% Self consistency loop
while abs(E2-E1) > tol
    E1 = E2;
    V_H = calc_potHA3(phi,r);   % V_sH
    n=2*abs(phi).^2;    % Electron density 2*n_s
    
    eps_x = -3/4*(3*n/pi).^(1/3);   % Exchange energy
    deps_x = -1/4*(3*n/pi).^(1/3);  % Derivative of exchange energy 
    V_x = eps_x + deps_x;   % Exchange potential
    V = 2*V_H + V_x;    %  Total potential
    
    A = hamiltonianHA3(V,r,Z);    
    [F,lambda] = eig(A);

    u = F(:,1); % Eigenfunction to radial Kohm-Sham
    phi = u./(sqrt(4*pi)*r);    % Wave function
    phi = phi / sqrt(trapz(r,4*pi*r.^2.*phi.^2)); % Normalization
    E2 = 2*lambda(1,1) - 2*trapz(u.^2.*(V_H + V_x - eps_x))   % Calculate energy
    i=i+1;
end
i
E = E2
plot(r,-phi,'blue')
xlabel('$r$ [a$_0$]','Interpreter','Latex','fontsize',16)
ylabel('$\phi_{1s}$ [a$_0^{-3/2}$] ','Interpreter','Latex','fontsize',16)
title('Wave function for He$_{1s}$','Interpreter','Latex','fontsize',16)
plot(r,phi,'r-.')

xlim([0 2])

%% 6
clc, %close all
% r
r_min = 0;
r_max = 10;
N = 2000;
r = linspace(r_min,r_max,N+2);
r = r(2:end-1)';
Z = 2;

% Parameters for correlation terms
A = 0.0311;
B = -0.048;
C = 0.0020;
D = -0.0116;
gamma = -0.1423;
beta1 = 1.0529;
beta2 = 0.03334;



phi=ones(N,1)/N;   % Initial guess for wave function
E1 = 0; % Set energies to enter while-loop
E2 = 3;
tol=1e-7;   % Convergence tolerance
Z=2;    % Helium
i = 0;  % Number of iterations
eps_c = zeros(N,1); % Correlation energy
V_c = zeros(N,1);   % Correlation potential
% Self consistency loop
while abs(E2-E1) > tol
    E1 = E2;
    V_H = calc_potHA3(phi,r); % V_sH
    n=2*abs(phi).^2;    % Electron density 2*n_s
    eps_x = -3/4*(3*n/pi).^(1/3);   % Exchange energy
    deps_x = -1/4*(3*n/pi).^(1/3);  % Derivative of exchange energy 
    V_x = eps_x + deps_x;   % Exchange potential
    r_s = (3./(4*pi*n)).^(1/3);
    
    % Correlation energy/potential
    for j = 1:length(r_s)
        if r_s(j) < 1
            eps_c(j) = A*log(r_s(j)) + B + C*r_s(j)*log(r_s(j)) + D*r_s(j);
            V_c(j) = A*log(r_s(j)) + B - A/3 + 2/3*C*r_s(j)*log(r_s(j)) + (2*D-C)*r_s(j)/3;
        else
            eps_c(j) = gamma/(1 + beta1*sqrt(r_s(j)) + beta2*r_s(j));
            V_c(j) = eps_c(j)*(1+7/6*beta1*sqrt(r_s(j))+beta2*r_s(j))/(1+beta1*sqrt(r_s(j))+beta2*r_s(j));           
        end
    end
    
    V = 2*V_H + V_x + V_c; %  Total potential
    
    H = hamiltonianHA3(V,r,Z);    
    [F,lambda] = eig(H);

    u = F(:,1);
    u=u/sqrt(trapz(r,u.^2)); % Eigenfunction to radial Kohm-Sham
    phi = u./(sqrt(4*pi)*r); % Wave function
    E2 = 2*lambda(1,1) - 2*trapz(r,u.^2.*(V_H + V_x + 2*V_c - eps_x - 2*eps_c)); % Calculate energy
    i=i+1;
end
i
E = E2

plot(r,-phi,'r')
xlim([0 2])
%% Grid points convergence
clc
N_grid = 30;
%h = linspace(0,0.5,N_grid);
h = zeros(N_grid,1);
E = zeros(N_grid,1);
N_min = 200;
N_max = 2000;

N_vect = linspace(N_min,N_max,N_grid);
r_min = 0;
r_max = 10;
for j = 1:N_grid
    N = round(N_vect(j));
    r = linspace(r_min,r_max,N+2);
    r = r(2:end-1)';
    h(j) = r(2) - r(1);

    phi=ones(N,1)/N;   % Initial guess for wave function
    E1 = 0; % Set energies to enter while-loop
    E2 = 3;
    tol=1e-7;   % Convergence tolerance
    Z=2;    % Helium
    i = 0;  % Number of iterations
    % Self consistency loop
    while abs(E2-E1) > tol
        E1 = E2;
        V = calc_potHA3(phi,r); % V_sH

        A = hamiltonianHA3(V,r,Z);    
        [F,lambda] = eig(A);

        u = F(:,1); % Eigenfunction to radial Kohm-Sham
        phi = u./(sqrt(4*pi)*r);      % Wave function  
        phi = phi / sqrt(trapz(r,4*pi*r.^2.*phi.^2));  % Normalization
        E2 = 2*lambda(1,1) - trapz(r,4*pi*r.^2.*V.*abs(phi).^2);  % Calculate energy  
        i = i+1;
    end
    E(j) = E2;
end
%%
clc
figure(1)
subplot(1,2,1)
plot(N_vect,E,'o-')
xlabel('N','Interpreter','Latex','fontsize',14)
ylabel('Ground-state energy [E$_H$] ','Interpreter','Latex','fontsize',14)
subplot(1,2,2)
plot(N_vect(end-8:end),E(end-8:end),'o-')
%xlim([1500 2000])
xlabel('N','Interpreter','Latex','fontsize',14)
%ylabel('Ground-state energy [E$_H$] ','Interpreter','Latex','fontsize',14)

figure(2)
subplot(1,2,1)
plot(h,E,'o-')
xlabel('Grid spacing [a$_0$]','Interpreter','Latex','fontsize',12)
ylabel('E [E$_H$] ','Interpreter','Latex','fontsize',12)
subplot(1,2,2)
plot(h(end-8:end),E(end-8:end),'o-')
xlabel('h [a$_0$]','Interpreter','Latex','fontsize',12)
ylabel('E [E$_H$] ','Interpreter','Latex','fontsize',12)
%% r conv
clc, close all
% r
Nr = 25;
r_c = zeros(Nr,1);
E = zeros(Nr,1);
for r_max=1:Nr
r_c(r_max) = r_max;;   
r_min = 0;
h=0.01;
N = r_max/h;
r = linspace(r_min,r_max,N+2);
r = r(2:end-1)';

phi=ones(N)/N;
n=phi.^2;
Z=2;

V = calc_potHA3(phi,r);

E1 = 0;
E2 = 3;
tol=1e-7;
i = 0;
while abs(E2-E1) > tol
    E1 = E2;
    V = calc_potHA3(phi,r);
    
    A = hamiltonianHA3(V,r,Z);    
     
    [F,lambda] = eig(A);
    
    f = F(:,1);
    phi = f./(sqrt(4*pi)*r);
    phi = phi / sqrt(trapz(r,(4*pi*r.^2.*phi.^2)));  
    E2 = 2*lambda(1,1) - trapz(r,4*pi*r.^2.*V.*abs(phi).^2);
   abs(E2-E1);
   i = i+1;
end
i;
E(r_max) = E2;
n = abs(phi).^2;
end
plot(r_c,E,'-o')
%%
plot(r_c,E,'-o')
xlabel('$r_{max}$ [a.u.]','Interpreter','Latex','fontsize',16)
ylabel('$E_0$ [a.u.]','Interpreter','Latex','fontsize',16)
title('Convergence for $E_0$','Interpreter','Latex','fontsize',16)

