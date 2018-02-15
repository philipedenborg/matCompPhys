%% 2
clc, close all
% r
r_min = 0;
r_max = 10;
N = 1000;
r = linspace(r_min,r_max,N+2);
r = r(2:end-1)';


phi=1/sqrt(pi)*exp(-r);

V = calc_potHA3(phi,r);

plot(r,V,'cyan')
hold on
V_H = 1./r - (1+1./r).*exp(-2*r);
plot(r,V_H,'--r') 

%% 3
clc, close all
% r
r_min = 0;
r_max = 10;
N = 1000;
r = linspace(r_min,r_max,N+2);
r = r(2:end-1)';
Z=1;

V=zeros(N);

H=hamiltonianHA3(V,r,Z);

[F,lambda] = eig(H);
phi = F(:,1);
%phi = u / sqrt(trapz(r,phi.^2))
%phi=phi/trapz(r,4*pi*r.^2.*phi.^2);

E=lambda(1,1)
n=abs(phi).^2;
n=n/trapz(r,n);
plot(r,n,'black')

%% 4
clc, close all
% r
r_min = 0;
r_max = 10;
N = 1000;
r = linspace(r_min,r_max,N+2);
r = r(2:end-1)';

phi=ones(N)/N;
n=phi.^2;
Z=2;

V = calc_potHA3(phi,r);

E1 = 0;
E2 = 3;
tol=1e-5;
i = 0;
while abs(E2-E1) > tol
    E1 = E2;
    V = calc_potHA3(phi,r);
    
    A = hamiltonianHA3(V,r,Z);    
     
    [F,lambda] = eig(A);
    
    f = F(:,1);
    phi = f./(sqrt(4*pi)*r);
    phi = phi / sqrt(trapz(r,4*pi*r.^2.*phi.^2));  
    E2 = 2*lambda(1,1) - trapz(r,4*pi*r.^2.*V.*abs(phi).^2);
   abs(E2-E1)
   i = i+1;
end
i
E = E2
n = abs(phi).^2;
n=2*n/trapz(r,4*pi*r.^2.*n);
plot(r,n,'g')
xlim([0 2])
hold on

%% 5
clc, %close all
% r
r_min = 0;
r_max = 10;
N = 1000;
r = linspace(r_min,r_max,N+2);
r = r(2:end-1)';

% Calc Hartree potential
phi=1/sqrt(pi)*exp(-r); % Hydrogen
V_H = calc_potHA3(phi,r); 


phi=ones(N,1)/N;
n=phi.^2;
eps_x = -3/4*(3*n/pi).^(1/3);

E1 = 0;
E2 = 3;
tol=1e-7;
Z=2;
i = 0;
V_H = 1./r - (1+1./r).*exp(-2*r);

while abs(E2-E1) > tol
    E1 = E2;
    V_H = calc_potHA3(phi,r);
    n=2*abs(phi).^2;
    eps_x = -3/4*(3*n/pi).^(1/3);
    deps_x = -1/4*(3*n/pi).^(1/3);
    V_x = eps_x + deps_x;
    V = 2*V_H + V_x;
    
    A = hamiltonianHA3(V,r,Z);    
    [F,lambda] = eig(A);

    u = F(:,1);
    phi = u./(sqrt(4*pi)*r);
    phi = phi / sqrt(trapz(4*pi*r.^2.*phi.^2));
    E2 = 2*lambda(1,1) - 2*trapz(u.^2.*(V_H + V_x - eps_x))   
    i=i+1;
end
i
E = E2

plot(r,phi,'r')
xlim([0 2])

%% 6
clc, %close all
% r
r_min = 0;
r_max = 15;
N = 2000;
r = linspace(r_min,r_max,N+2);
r = r(2:end-1)';

A = 0.0311;
B = -0.048;
C = 0.0020;
D = -0.0116;
gamma = -0.1423;
beta1 = 1.0529;
beta2 = 0.03334;

% Calc Hartree potential
phi=1/sqrt(pi)*exp(-r); % Hydrogen
V_H = calc_potHA3(phi,r); 


phi=ones(N,1)/N;
u = phi;
n=phi.^2;
eps_x = -3/4*(3*n/pi).^(1/3);

E1 = 0;
E2 = 3;
tol=1e-7;
Z=2;
i = 0;
V_H = 1./r - (1+1./r).*exp(-2*r);
eps_c = zeros(N,1);
V_c = zeros(N,1);
while abs(E2-E1) > tol
    E1 = E2;
    V_H = calc_potHA3(phi,r);
    n=2*abs(phi).^2;
    eps_x = -3/4*(3*n/pi).^(1/3);
    deps_x = -1/4*(3*n/pi).^(1/3);
    
    r_s = (3./(4*pi*n)).^(1/3);
    
    for j = 1:length(r_s)
        
        if r_s(j) < 1
            eps_c(j) = A*log(r_s(j)) + B + C*r_s(j)*log(r_s(j)) + D*r_s(j);
            V_c(j) = A*log(r_s(j)) + B - A/3 + 2/3*C*r_s(j)*log(r_s(j)) + (2*D-C)*r_s(j)/3;
        else
            eps_c(j) = gamma/(1 + beta1*sqrt(r_s(j)) + beta2*r_s(j));
            V_c(j) = eps_c(j)*(1+7/6*beta1*sqrt(r_s(j))+beta2*r_s(j))/(1+beta1*sqrt(r_s(j))+beta2*r_s(j));           
        end
    end
    
    V_x = eps_x + deps_x;
    V = 2*V_H + V_x + V_c;
    
    H = hamiltonianHA3(V,r,Z);    
    [F,lambda] = eig(H);

    u = F(:,1);
    u=u/sqrt(trapz(r,u.^2));
    phi = u./(sqrt(4*pi)*r);
    %phi = phi / sqrt(trapz(r,4*pi*r.^2.*phi.^2));
    E2 = 2*lambda(1,1) - 2*trapz(r,u.^2.*(V_H + V_x + V_c - eps_x - eps_c))   
    i=i+1;
end
i
E = E2

plot(r,phi,'r')
xlim([0 2])