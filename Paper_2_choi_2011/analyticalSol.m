function [ sol ] = analyticalSol()

data = nondimData();

N = 15;

%-------------------------------------------------------------------------%

% GENERATING RESULTS

% 1-D Mesh/Grid
sol.y1_star = linspace(0, data.h1_star, N);
sol.y2_star = linspace(data.h1_star, 1, N);
sol.y_star = [sol.y1_star sol.y2_star];

sol.y1_star = sol.y1_star';
sol.y2_star = sol.y2_star';
sol.y_star = sol.y_star';

% Electric Potential
[ sol.phi_star, sol.phi_1_star, sol.phi_2_star ] = ...
    electricPotential(data.K1, data.K2, data.h1_star, data.h2_star, ...
    data.zeta1_star, data.zeta2_star, data.zeta_int_star, ...
    data.alpha_e, data.q_int_star, ...
    sol.y1_star, sol.y2_star);

% Velocity Field
[ sol.u_star, sol.u_1_star, sol.u_2_star ] = ...
    velocityField(data.h1_star, data.zeta1_star, data.zeta2_star, ...
    data.alpha_e, sol.phi_1_star, sol.phi_2_star, ...
    sol.y1_star, sol.y2_star);

end

%-------------------------------------------------------------------------%

% ANALYTICAL FUNCTION

function [ phi_star, phi_1_star, phi_2_star ] = electricPotential(K1, K2, h1_star, h2_star, zeta1_star, zeta2_star, zeta_int_star, alpha_e, q_int_star, y1_star, y2_star)

f1 = K1 * coth(K1 * h1_star) * ( exp(K2 * (h1_star - 2)) - exp(-K2*h1_star) );
f2 = K1 * coth(K1 * h1_star) * ( zeta_int_star - (zeta2_star * exp(-K2*h2_star)) );
f3 = alpha_e * K2 * ( exp(K2 * (h1_star - 2)) + exp(-K2*h1_star) );
f4 = -q_int_star - (alpha_e*K2*zeta2_star*exp(-K2*h2_star));
f5 = K1 * zeta1_star * (-2 / (exp(K1 * h1_star) - exp(-K1 * h1_star)));

b2 = (f5 + f4 - f2) / (f1 - f3);

b1_1 = ( 1 / (2 * sinh(K1*h1_star)) ) * ...
    (zeta_int_star - ...
    (zeta2_star * exp(-K2*h2_star)) + ...
    ( b2 * ( exp(K2 * (h1_star - 2)) - exp(-K2*h1_star) ) ) );

b1_2 = zeta1_star / (1 - exp(-2 * K1 * h1_star));

b1 = b1_1 + b1_2;

a2 = (zeta2_star * exp(-K2)) - (b2 * exp(-2*K2));

a1 = (zeta1_star - b1_2) - b1_1;

phi_1_star = (a1 .* exp(K1 .* y1_star)) + (b1 .* exp(-K1 .* y1_star));

phi_2_star = (a2 .* exp(K2 .* y2_star)) + (b2 .* exp(-K2 .* y2_star));

phi_star = [phi_1_star; phi_2_star];

end

function [ u_star, u_1_star, u_2_star ] = velocityField(h1_star, zeta1_star, zeta2_star, alpha_e, phi_1_star, phi_2_star, y1_star, y2_star)

Cf = alpha_e * zeta2_star / zeta1_star;

C0 = ( (phi_1_star(end) / zeta1_star) - (phi_2_star(1) / zeta2_star) );

C3 = C0 / ( 1 + ( h1_star * (Cf - 1) ) );

C1 = Cf * C3;

C4 = 1 - C3;

C2 = 1;

u_1_star = (-1 .* phi_1_star ./ zeta1_star) + (C1 .* y1_star) + C2;

u_2_star = (-1 .* phi_2_star ./ zeta2_star) + (C3 .* y2_star) + C4;

u_star = [u_1_star; u_2_star];

end

%-------------------------------------------------------------------------%