function [ output ] = nondimData()

data = readstruct("dimData.xml");

%-------------------------------------------------------------------------%

% DIMENSIONLESS FORM

% Grid Info

output.N_c = data.fluid_c.Nc;

output.N_a = data.fluid_a.Na;

% Interface

output.r_o_star = data.system.r_o / data.system.R;

% Debye Length

output.lambda_D_c = (1 / (data.fluid_c.z_c * data.const.e)) * ...
    sqrt((data.const.eps_o * data.const.k_b * data.system.T * data.fluid_c.eps_c) ...
    / (2 * data.fluid_c.n_c_inf));

output.lambda_D_a = (1 / (data.fluid_a.z_a * data.const.e)) * ...
    sqrt((data.const.eps_o * data.const.k_b * data.system.T * data.fluid_a.eps_a) ...
    / (2 * data.fluid_a.n_a_inf));

% Debye-Huckel Parameter

output.Kc = data.system.R / output.lambda_D_c;

output.Ka = data.system.R / output.lambda_D_a;

% Elementary Thermal Voltage

phi_T = data.const.k_b * data.system.T / data.const.e;

% Electrostatic BC variables

output.zeta_wall_star = data.system.zeta_wall / phi_T;

output.beta = data.fluid_c.eps_c / data.fluid_a.eps_a;

output.zeta_int_star = data.system.zeta_int / phi_T;

output.q_int_star = (data.system.q_int * data.system.R / ...
    (data.fluid_a.eps_a * data.const.eps_o)) / phi_T;

% Hydrodynamic variables

output.u_h = -1 * data.fluid_a.eps_a * data.const.eps_o * ...
    data.system.E_ext * data.system.zeta_wall / data.fluid_a.mu_a;

%-------------------------------------------------------------------------%

writestruct(output, "nondimData.xml")

end