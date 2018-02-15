function A = hamiltonianHA3(V,r,Z)
N = length(r);
A = zeros(N,N);
h = r(2) - r(1);
for i = 2:N-1
   A(i,i) = 1/h^2 -Z/r(i) + V(i);
   A(i,i-1) = -1/(2*h^2);
   A(i,i+1) = -1/(2*h^2);   
end
A(1,1) = 1/h^2 - Z/r(1) + V(1);
A(N,N) = 1/h^2 - Z/r(N) + V(N);
A(N,N-1) = -1/(2*h^2);
A(1,2) = -1/(2*h^2);

end