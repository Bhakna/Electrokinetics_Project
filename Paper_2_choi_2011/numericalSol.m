function [ sol ] = numericalSol()

data = nondimData();

%-------------------------------------------------------------------------%

% GENERATING RESULTS

% 1-D Mesh/Grid
sol.y1_star = linspace(0, data.h1_star, data.N1);
sol.y2_star = linspace(data.h1_star, 1, data.N2);
sol.y_star = [sol.y1_star sol.y2_star];

sol.y1_star = sol.y1_star';
sol.y2_star = sol.y2_star';
sol.y_star = sol.y_star';

% Electric Potential
[ sol.phi_star, sol.phi_1_star, sol.phi_2_star ] = ...
    electricPotential(data.K1, data.K2, ...
    data.zeta1_star, data.zeta2_star, data.zeta_int_star, ...
    data.alpha_e, data.q_int_star, ...
    sol.y1_star, sol.y2_star);

% Velocity Field
[ sol.u_star, sol.u_1_star, sol.u_2_star ] = ...
    velocityField(data.K1, data.K2, ...
    data.zeta1_star, data.zeta2_star, ...
    data.alpha_e, data.q_int_star, sol.phi_star, ...
    sol.y1_star, sol.y2_star);

end

%-------------------------------------------------------------------------%

% ANALYTICAL FUNCTION

function [ phi_star, phi_1_star, phi_2_star ] = electricPotential(K1, K2, zeta1_star, zeta2_star, zeta_int_star, alpha_e, q_int_star, y1_star, y2_star)

del_y1_star = y1_star(2) - y1_star(1);
N1 = length(y1_star);

del_y2_star = y2_star(2) - y2_star(1);
N2 = length(y2_star);

beta_1 = alpha_e * del_y1_star / del_y2_star;
beta_2 = q_int_star * del_y1_star;

% Initializing matrix

A = zeros((N1+N2), (N1+N2));
b = zeros((N1+N2), 1);

% Bulk part

for i = 2:(N1-1)

    A(i, i-1) = 1;
    A(i, i) = -2 - ((K1 * del_y1_star) ^ 2);
    A(i, i+1) = 1;

end

for i = (N1+2):(N1+N2-1)

    A(i, i-1) = 1;
    A(i, i) = -2 - ((K2 * del_y2_star) ^ 2);
    A(i, i+1) = 1;

end

% Boundary part

A(1, 1) = 1;
b(1) = zeta1_star;

A(end, end) = 1;
b(end) = zeta2_star;

A(N1, N1) = -1;
A(N1, N1+1) = 1;
b(N1) = zeta_int_star;

A(N1+1, N1-1) = -1;
A(N1+1, N1) = 1;
A(N1+1, N1+1) = beta_1;
A(N1+1, N1+2) = -1 * beta_1;
b(N1+1) = beta_2;

% Final solution

phi_star = inv(A) * b;

phi_1_star = phi_star(1:N1);

phi_2_star = phi_star((N1+1):end);

end

function [ u_star, u_1_star, u_2_star ] = velocityField(K1, K2, zeta1_star, zeta2_star, alpha_e, q_int_star, phi_star, y1_star, y2_star)

del_y1_star = y1_star(2) - y1_star(1);
N1 = length(y1_star);

del_y2_star = y2_star(2) - y2_star(1);
N2 = length(y2_star);

beta_1 = alpha_e * del_y1_star / del_y2_star;
beta_2 = q_int_star * del_y1_star;

% Initializing matrix

A = zeros((N1+N2), (N1+N2));
b = zeros((N1+N2), 1);

% Bulk part

for i = 2:(N1-1)

    A(i, i-1) = 1;
    A(i, i) = -2;
    A(i, i+1) = 1;
    b(i) = (-1 * phi_star(i) / zeta1_star) * ((K1 * del_y1_star)^2);

end

for i = (N1+2):(N1+N2-1)

    A(i, i-1) = 1;
    A(i, i) = -2;
    A(i, i+1) = 1;
    b(i) = (-1 * phi_star(i) / zeta2_star) * ((K2 * del_y2_star)^2);

end

% Boundary part

A(1, 1) = 1;
b(1) = 0;

A(end, end) = 1;
b(end) = 0;

A(N1, N1) = -1;
A(N1, N1+1) = 1;

A(N1+1, N1-1) = -1;
A(N1+1, N1) = 1;
A(N1+1, N1+1) = beta_1 * zeta2_star / zeta1_star;
A(N1+1, N1+2) = -1 * beta_1 * zeta2_star / zeta1_star;
b(N1+1) = -1 * beta_2 / zeta1_star;

% Final solution

u_star = inv(A) * b;

u_1_star = u_star(1:N1);

u_2_star = u_star((N1+1):end);

end

%-------------------------------------------------------------------------%