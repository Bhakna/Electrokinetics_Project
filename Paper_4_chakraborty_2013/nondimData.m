function [ output ] = nondimData()

data = readstruct("dimData.xml");

%-------------------------------------------------------------------------%

% DIMENSIONLESS FORM

% Grid Info

output.N_y = data.system.N_y;

output.N_z = data.system.N_z;

% Interface

output.D_h = 4 * data.system.L * data.system.H_int / ...
    (data.system.L + data.system.H_int);

output.h1_star = data.system.H_int / output.D_h;

output.l_star = data.system.L / output.D_h;

% Debye number (lambda_D)

output.lambda_D = sqrt( data.fluid.eps * data.const.eps_o * ...
    data.const.k_b * data.system.T / ...
    (2 * data.fluid.n_inf * ...
    ((data.const.e * data.fluid.z)^2)) );

% Debye number (De)

output.De = output.lambda_D / output.D_h;

% Electrostatic variables

output.Z_b_star = data.system.Z_b / data.system.Z_b;

output.Z_s_star = data.system.Z_s / data.system.Z_b;

output.Z_int_star = data.system.Z_int / data.system.Z_b;

% Hydrodynamic variables

output.u_ref = -1 * data.fluid.eps * data.const.eps_o * ...
    data.system.E_ext * data.system.Z_b / data.fluid.mu;

output.gamma_r = data.fluid.eps * data.const.eps_o * ...
    data.system.E_ext * data.system.Z_b / (data.fluid.mu * output.u_ref);

% Schemes

output.anaSolFunc = data.scheme.analyticalSol;

output.numSolFunc = data.scheme.numericalSol;

output.numSORw = data.scheme.numSORw;

%-------------------------------------------------------------------------%

writestruct(output, "nondimData.xml")

end