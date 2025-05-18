function [ sol ] = numericalSol(data)

const = readstruct("constants.xml");

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
sol.y1_star = linspace(0, h1_star, data.N);
sol.y2_star = linspace(h1_star, 1, data.N);
sol.y_star = [sol.y1_star sol.y2_star(2:end)];

sol.y1_star = sol.y1_star';
sol.y2_star = sol.y2_star';
sol.y_star = sol.y_star';

% Electric Potential
sol.phi_star = dimElecPotential(K, h1_star, zeta1_star, zeta2_star, data.N);

% Fluid Velocity
sol.u_star = dimVelocityField(K, C0, C1, C2, data.alpha, h1_star, data.N, sol.phi_star, data.BC_interface);
sol.u1_star = sol.u_star(1:data.N);
sol.u2_star = sol.u_star(data.N:end);

% Interface Velocity
sol.u_int_star = sol.u_star(data.N);

% Interface Shear Stress
sol.tau_int_star = (sol.u1_star(data.N) - sol.u1_star(data.N - 1)) / ...
    (sol.y1_star(data.N) - sol.y1_star(data.N - 1));

% Volumetric Flowrate
sol.Q_star = dimVolFlowRate(h1_star, data.N, sol.u_star);

% Mean Velocity
sol.u_mean = [(sol.Q_star(1) / h1_star), ...
    (sol.Q_star(2) / (1 - h1_star))]';

%-------------------------------------------------------------------------%

end