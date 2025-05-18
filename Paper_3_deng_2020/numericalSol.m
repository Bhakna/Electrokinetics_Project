function [ sol ] = numericalSol()

data = nondimData();

%-------------------------------------------------------------------------%

% GENERATING RESULTS

% 1-D Mesh/Grid
sol.r_c_star = linspace(0, data.r_o_star, data.N_c);
sol.r_a_star = linspace(data.r_o_star, 1, data.N_a);
sol.r_star = [sol.r_c_star sol.r_a_star];

sol.r_c_star = sol.r_c_star';
sol.r_a_star = sol.r_a_star';
sol.r_star = sol.r_star';

% Electric Potential
[ sol.phi_star, sol.phi_c_star, sol.phi_a_star ] = ...
    electricPotential(data.Kc, data.Ka, ...
    data.zeta_wall_star, data.zeta_int_star, ...
    data.beta, data.q_int_star, ...
    sol.r_c_star, sol.r_a_star);

% Velocity Field
[ sol.u_star, sol.u_c_star, sol.u_a_star ] = ...
    velocityField(data.Kc, data.Ka, ...
    data.zeta_wall_star, data.beta, data.q_int_star, ...
    sol.phi_c_star, sol.phi_a_star, ...
    sol.r_c_star, sol.r_a_star);

end

%-------------------------------------------------------------------------%

% ANALYTICAL FUNCTION

function [ phi_star, phi_c_star, phi_a_star ] = electricPotential(Kc, Ka, zeta_wall_star, zeta_int_star, beta, q_int_star, r_c_star, r_a_star)

del_r_c_star = r_c_star(2) - r_c_star(1);
N_c = length(r_c_star);

del_r_a_star = r_a_star(2) - r_a_star(1);
N_a = length(r_a_star);

% Initializing matrix

A = zeros((N_c + N_a), (N_c + N_a));
b = zeros((N_c + N_a), 1);

% Bulk part

for i = 2:(N_c - 1)

    A(i, i-1) = (2 * r_c_star(i)) - del_r_c_star;
    A(i, i) = -2 * r_c_star(i) * (2 + ((Kc * del_r_c_star) ^ 2));
    A(i, i+1) = (2 * r_c_star(i)) + del_r_c_star;

end

for i = (N_c + 2):(N_c + N_a - 1)

    A(i, i-1) = (2 * r_a_star(i - N_c)) - del_r_a_star;
    A(i, i) = -2 * r_a_star(i - N_c) * (2 + ((Ka * del_r_a_star) ^ 2));
    A(i, i+1) = (2 * r_a_star(i - N_c)) + del_r_a_star;

end

% Boundary part

A(1, 1) = 1;
A(1, 2) = -1;

A(end, end) = 1;
b(end) = zeta_wall_star;

A(N_c, N_c) = 1;
A(N_c, N_c+1) = -1;
b(N_c) = zeta_int_star;

A(N_c+1, N_c-1) = -1 * beta * del_r_a_star;
A(N_c+1, N_c) = beta * del_r_a_star;
A(N_c+1, N_c+1) = del_r_c_star;
A(N_c+1, N_c+2) = -1 * del_r_c_star;
b(N_c+1) = q_int_star * del_r_c_star * del_r_a_star;

% Final solution

phi_star = A\b;

phi_c_star = phi_star(1:N_c);

phi_a_star = phi_star((N_c+1):end);

end

function [ u_star, u_c_star, u_a_star ] = velocityField(Kc, Ka, zeta_wall_star, beta, q_int_star, phi_c_star, phi_a_star, r_c_star, r_a_star)

del_r_c_star = r_c_star(2) - r_c_star(1);
N_c = length(r_c_star);

del_r_a_star = r_a_star(2) - r_a_star(1);
N_a = length(r_a_star);

% Initializing matrix

A = zeros((N_c + N_a), (N_c + N_a));
b = zeros((N_c + N_a), 1);

% Bulk part

for i = 2:(N_c - 1)

    A(i, i-1) = (2 * r_c_star(i)) - del_r_c_star;
    A(i, i) = -4 * r_c_star(i);
    A(i, i+1) = (2 * r_c_star(i)) + del_r_c_star;
    b(i) = -2 * r_c_star(i) * (del_r_c_star^2) * ...
        (Kc^2) * phi_c_star(i) / zeta_wall_star;

end

for i = (N_c + 2):(N_c + N_a - 1)

    A(i, i-1) = (2 * r_a_star(i - N_c)) - del_r_a_star;
    A(i, i) = -4 * r_a_star(i - N_c);
    A(i, i+1) = (2 * r_a_star(i - N_c)) + del_r_a_star;
    b(i) = -2 * r_a_star(i - N_c) * (del_r_a_star^2) * ...
        (Ka^2) * phi_a_star(i - N_c) / zeta_wall_star;

end

% Boundary part

% r = 0
A(1, 1) = 1;
A(1, 2) = -1;

% no slip (r = R)
A(end, end) = 1;
b(end) = 0;

% no slip (r = r_o)
A(N_c, N_c) = 1;
A(N_c, N_c+1) = -1;
b(N_c) = 0;

% stress equilibria
A(N_c+1, N_c-1) = beta * del_r_a_star;
A(N_c+1, N_c) = -1 * beta * del_r_a_star;
A(N_c+1, N_c+1) = -1 * del_r_c_star;
A(N_c+1, N_c+2) = del_r_c_star;
b(N_c+1) = q_int_star * del_r_c_star * del_r_a_star / zeta_wall_star;

% Final solution

u_star = A\b;

u_c_star = u_star(1:N_c);

u_a_star = u_star((N_c+1):end);

end

%-------------------------------------------------------------------------%