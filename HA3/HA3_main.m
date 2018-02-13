clc
xi=@(r)[exp(-0.297104*r) exp(-1.236745*r) exp(-5.749982*r) exp(-38.216677*r)];
alpha=[0.297104 1.236745 5.749982 38.216677];
r = linspace(0,100,1000);
xi1=@(r)exp(-alpha(1)*r);
xi2=@(r)exp(-alpha(2)*r);
xi3=@(r)exp(-alpha(3)*r);
xi4=@(r)exp(-alpha(4)*r);

xi = {@(r)exp(-alpha(1)*r),  @(r)exp(-alpha(2)*r), @(r)exp(-alpha(3)*r), @(r)exp(-alpha(4)*r)};

S = zeros(4,4);
for i=1:4
    for j=1:4
        S(i,j) = 2*pi/(alpha(i)+alpha(j));
    end
end

h = zeros(4,4);

for i=1:4
    for j=1:4
        h_func=@(r) xi{j}(r).*(4*pi*xi{i}(r).*alpha(i).*(alpha(i)*r-1) - 2*xi{i}(r)/r);
        haids = h_func(r);
        h(i,j) = trapz(r,haids);
    end
end

Q_prsq = zeros(4,4,4,4);
for i=1:4
    for j=1:4
        for k=1:4
            for l=1:4
                Q_prsq(i,j,k,l) = 2*pi^(5/2)/((alpha(i) + alpha(j))*(alpha(k) + alpha(l))*sqrt(alpha(i) + alpha(j) + alpha(k) +alpha(l)));
            end
        end
    end
end
                 
%%
clc
f1 = @(x) sin(x)+1;
f2 =  @(x) x*f1(x)+1;
trapz(r,f1)           