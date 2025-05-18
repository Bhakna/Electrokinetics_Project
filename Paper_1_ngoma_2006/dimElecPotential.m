function [ phi_star ] = dimElecPotential(K, h1_star, zeta1_star, zeta2_star, N)

%-------------------------------------------------------------------------%

% ELECTRIC POTENTIAL ( NUMERICAL )

% Finite Difference Method

ymesh = linspace(0, h1_star, N);
del_y = ymesh(2) - ymesh(1);

ymesh = ymesh';

C = 2 + ((K * del_y)^2);

A = zeros(N, N);
b = zeros(N, 1);

A(1, 1) = 1;
A(N, N) = 1;

for i = 2:(N-1)
    A(i, i-1) = 1;
    A(i, i) = -C;
    A(i, i+1) = 1;
end

b(1) = zeta1_star;
b(end) = zeta2_star;

phi_star = inv(A) * b;

end

%-------------------------------------------------------------------------%