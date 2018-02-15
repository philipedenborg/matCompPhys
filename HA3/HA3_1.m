format long
clc

alpha=[0.297104 1.236745 5.749982 38.216677]; %


E_old = -20;
E_new = -12; 
tol = 1e-5; 
count = 0; 

S = zeros(4,4);
for i=1:4
    for j=1:4
        S(i,j) = (pi/(alpha(i) + alpha(j)))^(3/2);
    end
end

h = zeros(4,4);

for i=1:4
    for j=1:4
        h(i,j) = 3*alpha(i)*alpha(j)*pi^(3/2)/(alpha(i)+alpha(j))^(5/2) - 4*pi/(alpha(i)+alpha(j));
    end
end

Q = zeros(4,4,4,4);
for i=1:4
    for j=1:4
        for k=1:4
            for l=1:4
                Q(i,j,k,l) = 2*pi^(5/2)/((alpha(i) + alpha(j))*(alpha(k) + alpha(l))*sqrt(alpha(i) + alpha(j) + alpha(k) +alpha(l)));
            end
        end
    end
end    

C = [1, 1, 1, 1]; %Initial guess for C
C_n = 0; %Sum to normalize C
for i=1:4
    for j=1:4
        C_n = C_n + C(i)*C(j)*S(i,j);
    end
end
C = C/sqrt(C_n); %Normalizing C

F = zeros(4,4);

while abs(E_new - E_old) > tol %Iterate until tolerance is met
    E_old = E_new;
    E_new = 0; 
    
    for i=1:4
        for j=1:4
            F(i,j) = h(i,j); 
            for k=1:4
                for l=1:4
                    F(i,j) = F(i,j) + Q(i,j,k,l)*C(k)*C(l); %Calc F for every iteration
                end
            end
        end
    end

    [W lambda V] = eig(F,S); %Solve the eigenvalue problem
    C = V(:,1); %Set c as eigenvector corresponding to lowest eigenvalue 
        C_n = 0;
    for i=1:4
        for j=1:4
            C_n = C_n + C(i)*C(j)*S(i,j);
        end
    end
    
    C = C/sqrt(C_n); %Normalize eigenvector

    for i=1:4
        for j=1:4
            E_new = E_new + 2*C(i)*C(j)*h(i,j);
            for l=1:4
                for k=1:4
                    E_new = E_new + Q(i,j,k,l)*C(i)*C(j)*C(k)*C(l); %Calc E for this iteration
                end
            end
        end
    end
    
    count = count + 1;

end

E_new

E_old