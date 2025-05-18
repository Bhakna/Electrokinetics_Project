function [ sol ] = analyticalSol(data)

const = readstruct("constants.xml");
N = 30;

%-------------------------------------------------------------------------%

% DIMENSIONLESS FORM

% Height of Interface
h1_star = data.h1 / data.h;

% Dybe-Heckel Parameter
kappa = data.z_nott * const.e * (( (2 * data.n_inf) / (data.eps * const.eps_nott * const.k_b * data.T) )^0.5);
K = kappa * data.h;

% Elementary Thermal Voltage
phi_T = const.k_b * data.T / (data.z_nott * const.e);

zeta1_star = data.zeta1 / phi_T;
zeta2_star = data.zeta2 / phi_T;

% Helmholtz–Smoluchowski electroosmotic velocity
u_h = data.eps * const.eps_nott * data.E_ext * const.k_b * data.T ...
    / (data.mu * data.z_nott * const.e);

% Pressure per Length
Px = data.del_p / data.L;

C0 = (data.h^2) * Px / ( data.mu * data.u_nott );
C1 = u_h / data.u_nott;
C2 = (data.h^2) * Px / ( data.alpha * data.mu * data.u_nott );

%-------------------------------------------------------------------------%

% GENERATING RESULTS

% 1-D Mesh/Grid
sol.y1_star = linspace(0, h1_star, N);
sol.y2_star = linspace(h1_star, 1, N);
sol.y_star = [sol.y1_star sol.y2_star(2:end)];

sol.y1_star = sol.y1_star';
sol.y2_star = sol.y2_star';
sol.y_star = sol.y_star';

% Electric Potential
sol.phi_star = electricPotential(K, h1_star, zeta1_star, zeta2_star, sol.y1_star);

% Fluid Velocity
sol.u_star = velocityProfile(K, C0, C1, C2, data.alpha, h1_star, zeta1_star, zeta2_star, sol.phi_star, N, sol.y_star);
sol.u1_star = sol.u_star(1:N);
sol.u2_star = sol.u_star(N:end);

% Interface Velocity
sol.u_int_star = sol.u_star(N);

% Interface Shear Stress
sol.tau_int_star = interfaceShearStress(K, C0, C1, C2, data.alpha, h1_star, zeta1_star, zeta2_star);

% Volumetric Flowrate
sol.Q_star = volumetricFlowRate(K, C0, C1, C2, data.alpha, h1_star, zeta1_star, zeta2_star);

% Mean Velocity
sol.u_mean = [(sol.Q_star(1) / h1_star), ...
    (sol.Q_star(2) / (1 - h1_star))]';

end

%-------------------------------------------------------------------------%

% ANALYTICAL FUNCTION

function [ output ] = electricPotential(K, h1_star, zeta1_star, zeta2_star, y_star)

var1 = zeta2_star .* sinh( K .* y_star );
var2 = zeta1_star .* sinh( K .* (y_star - h1_star) );
var3 = sinh( K * h1_star );

output = (var1 - var2) ./ var3;

end

function [ output ] = velocityProfile(K, C0, C1, C2, alpha, h1_star, zeta1_star, zeta2_star, phi_star, N, ymesh)

a2 = -C1 * zeta1_star;

den = 1 + (h1_star * (alpha - 1));
psi_h1_star = (K * h1_star / sinh(K * h1_star)) * ...
    ((zeta2_star * cosh(K * h1_star)) - zeta1_star);

b1 = (1/den) * ((C1 * psi_h1_star) + (alpha*C2*(h1_star^2)) + (C2/2) - ...
    ((h1_star^2)*(C0+C2)/2) + (C1*(zeta1_star - zeta2_star)));

a1 = (alpha*b1) + (C0*h1_star) - (alpha*C2*h1_star) - (C1 * psi_h1_star / h1_star);

b2 = (C2/2) - b1;

output = zeros((2*N)-1, 1);

y_star = ymesh(1:N);
output(1:N) = (a1 .* y_star) + a2 + (C1 .* phi_star) - ((C0./2) .* (y_star.^2));

y_star = ymesh(N+1:end);
output(N+1:end) = (b1 .* y_star) + b2 - ((C2./2) .* (y_star.^2));

end

function [ output ] = interfaceShearStress(K, C0, C1, C2, alpha, h1_star, zeta1_star, zeta2_star)

den = 1 + (h1_star * (alpha - 1));
psi_h1_star = (K * h1_star / sinh(K * h1_star)) * ...
    ((zeta2_star * cosh(K * h1_star)) - zeta1_star);

b1 = (1/den) * ((C1 * psi_h1_star) + (alpha*C2*(h1_star^2)) + (C2/2) - ...
    ((h1_star^2)*(C0+C2)/2) + (C1*(zeta1_star - zeta2_star)));

output = alpha * (b1 - (C2 * h1_star));

end

function [ output ] = volumetricFlowRate(K, C0, C1, C2, alpha, h1_star, zeta1_star, zeta2_star)

a2 = -C1 * zeta1_star;

den = 1 + (h1_star * (alpha - 1));
psi_h1_star = (K * h1_star / sinh(K * h1_star)) * ...
    ((zeta2_star * cosh(K * h1_star)) - zeta1_star);

b1 = (1/den) * ((C1 * psi_h1_star) + (alpha*C2*(h1_star^2)) + (C2/2) - ...
    ((h1_star^2)*(C0+C2)/2) + (C1*(zeta1_star - zeta2_star)));

a1 = (alpha*b1) + (C0*h1_star) - (alpha*C2*h1_star) - (C1 * psi_h1_star / h1_star);

b2 = (C2/2) - b1;

var1 = (C1 / (K * sinh(K * h1_star))) * ...
    ((zeta2_star * (cosh(K * h1_star) - 1)) - ...
    (zeta1_star * (1 - cosh(K * h1_star))));

Q1_star = (-C0 * (h1_star^3) / 6) + (a1 * (h1_star^2) / 2) + (a2*h1_star) + var1;

Q2_star = ((-C2/6)*(1-(h1_star^3))) + ((b1/2)*(1-(h1_star^2))) + (b2*(1-h1_star));

output = [Q1_star Q2_star]';

end

%-------------------------------------------------------------------------%