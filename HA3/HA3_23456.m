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
phi=
E=lambda(1)
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
n=n/trapz(r,4*pi*r.^2.*n);
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
    phi = phi / sqrt(trapz(r,4*pi*r.^2.*phi.^2));
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
    
    r_s= (4*pi*n/3)^(1/3);
    
    if r_s < 1
        
    end
    
    V_x = eps_x + deps_x;
    V = 2*V_H + V_x;
    
    A = hamiltonianHA3(V,r,Z);    
    [F,lambda] = eig(A);

    u = F(:,1);
    phi = u./(sqrt(4*pi)*r);
    phi = phi / sqrt(trapz(r,4*pi*r.^2.*phi.^2));
    E2 = 2*lambda(1,1) - 2*trapz(u.^2.*(V_H + V_x - eps_x))   
    i=i+1;
end
i
E = E2

plot(r,phi,'r')
xlim([0 2])