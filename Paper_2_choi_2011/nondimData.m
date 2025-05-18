function [ output ] = nondimData()

data = readstruct("dimData.xml");

%-------------------------------------------------------------------------%

% DIMENSIONLESS FORM

% Grid Info

output.N1 = data.fluid1.N1;

output.N2 = data.fluid2.N2;

% Height

output.h1_star = data.fluid1.h1 / data.system.h;

output.h2_star = data.fluid2.h2 / data.system.h;

% Dybe-Heckel Parameter

output.K1 = data.system.h * data.fluid1.z1 * data.const.e * ...
    (( (2 * data.fluid1.n1_inf) / ...
    (data.fluid1.eps_r1 * data.const.eps_nott * data.const.k_b * data.system.T) )^0.5);

output.K2 = data.system.h * data.fluid2.z2 * data.const.e * ...
    (( (2 * data.fluid2.n2_inf) / ...
    (data.fluid2.eps_r2 * data.const.eps_nott * data.const.k_b * data.system.T) )^0.5);

% Elementary Thermal Voltage

phi_T = data.const.k_b * data.system.T / data.const.e;

% Electrostatic BC variables

output.zeta1_star = data.fluid1.zeta_1 / phi_T;

output.zeta2_star = data.fluid2.zeta_2 / phi_T;

output.zeta_int_star = data.system.zeta_int / phi_T;

output.alpha_e = data.fluid2.eps_r2 / data.fluid1.eps_r1;

output.q_int_star = (data.system.q_int * data.system.h / ...
    (data.fluid1.eps_r1 * data.const.eps_nott)) / phi_T;

output.q_1_star = (data.fluid1.q_1 * data.system.h / ...
    (data.fluid1.eps_r1 * data.const.eps_nott)) / phi_T;

output.q_2_star = (data.fluid2.q_2 * data.system.h / ...
    (data.fluid2.eps_r2 * data.const.eps_nott)) / phi_T;

% Hydrodynamic variables

output.u_h1 = -1 * data.fluid1.eps_r1 * data.const.eps_nott * ...
    data.system.E_ext * data.fluid1.zeta_1 / data.fluid1.mu_1;

output.u_h2 = -1 * data.fluid2.eps_r2 * data.const.eps_nott * ...
    data.system.E_ext * data.fluid2.zeta_2 / data.fluid2.mu_2;

%-------------------------------------------------------------------------%

writestruct(output, "nondimData.xml")

end