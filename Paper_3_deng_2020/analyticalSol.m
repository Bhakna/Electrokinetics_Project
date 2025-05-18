function [ sol ] = analyticalSol()

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
    electricPotential(data.Kc, data.Ka, data.r_o_star, ...
    data.zeta_wall_star, data.zeta_int_star, ...
    data.beta, data.q_int_star, ...
    sol.r_c_star, sol.r_a_star);

% Velocity Field
[ sol.u_star, sol.u_c_star, sol.u_a_star ] = ...
    velocityField(data.zeta_int_star, data.zeta_wall_star, ...
    sol.phi_c_star, sol.phi_a_star);

end

%-------------------------------------------------------------------------%

% ANALYTICAL FUNCTION

function [ phi_star, phi_c_star, phi_a_star ] = electricPotential(Kc, Ka, r_o_star, zeta_wall_star, zeta_int_star, beta, q_int_star, r_c_star, r_a_star)

D = (Kc * beta * besseli(1, (Kc * r_o_star)) * ...
    ( (-1 * besseli(0, Ka*r_o_star) * besselk(0, Ka)) + (besseli(0, Ka) * besselk(0, Ka*r_o_star)) )) + ...
    (Ka * besseli(0, (Kc * r_o_star)) * ...
    ( (besseli(1, Ka*r_o_star) * besselk(0, Ka)) + (besseli(0, Ka) * besselk(1, Ka*r_o_star)) ));

An = (Kc * beta * besseli(1, (Kc * r_o_star)) * ...
    ( (zeta_int_star * besselk(0, Ka)) + (zeta_wall_star * besselk(0, Ka*r_o_star)) )) + ...
    (besseli(0, (Kc * r_o_star)) * ...
    ( (-1 * q_int_star * besselk(0, Ka)) + (Ka * zeta_wall_star * besselk(1, Ka*r_o_star)) ));

A = An / D;

Bn = (-1 * Kc * beta * besseli(1, (Kc * r_o_star)) * ...
    ( (zeta_int_star * besseli(0, Ka)) + (zeta_wall_star * besseli(0, Ka*r_o_star)) )) + ...
    (besseli(0, (Kc * r_o_star)) * ...
    ( (q_int_star * besseli(0, Ka)) + (Ka * zeta_wall_star * besseli(1, Ka*r_o_star)) ));

B = Bn / D;

Cn = (zeta_wall_star / r_o_star) + ...
    (besselk(0, Ka) * ...
    ( (-1 * q_int_star * besseli(0, Ka*r_o_star)) + (Ka * zeta_int_star * besseli(1, Ka*r_o_star)) )) + ...
    (besseli(0, Ka) * ...
    ( (q_int_star * besselk(0, Ka*r_o_star)) + (Ka * zeta_int_star * besselk(1, Ka*r_o_star)) ));

C = Cn / D;

phi_c_star = C .* besseli(0, (Kc .* r_c_star));

phi_a_star = (A .* besseli(0, (Ka .* r_a_star))) + ...
    (B .* besselk(0, (Ka .* r_a_star)));

phi_star = [phi_c_star; phi_a_star];

end

function [ u_star, u_c_star, u_a_star ] = velocityField(zeta_int_star, zeta_wall_star, phi_c_star, phi_a_star)

u_c_star = 1 + (zeta_int_star / zeta_wall_star) - (phi_c_star ./ zeta_wall_star);

u_a_star = 1 - (phi_a_star ./ zeta_wall_star);

u_star = [u_c_star; u_a_star];

end

%-------------------------------------------------------------------------%