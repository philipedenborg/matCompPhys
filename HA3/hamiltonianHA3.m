function H = hamiltonianHA3(V,r,Z)
% Returns hamiltonian 
%
%

N = length(r);
H = zeros(N,N);
h = r(2) - r(1);    % Stepping distance 
for i = 2:N-1
   H(i,i) = 1/h^2 -Z/r(i) + V(i);   % Kinetic part + nuclear part + potential
   H(i,i-1) = -1/(2*h^2);
   H(i,i+1) = -1/(2*h^2);   
end
H(1,1) = 1/h^2 - Z/r(1) + V(1);
H(N,N) = 1/h^2 - Z/r(N) + V(N);
H(N,N-1) = -1/(2*h^2);
H(1,2) = -1/(2*h^2);

end