%-------------------------------------------------------------------------%

clc;
clear variables;

%-------------------------------------------------------------------------%

data_const = struct;

% universal constants
data_const.eps_o = 8.854e-12;
data_const.e = 1.6022e-19;
data_const.k_b = 1.3807e-23;
data_const.N_A = 6.022e23;

%-------------------------------------------------------------------------%

data_scheme = struct;

data_scheme.analyticalSol = "mathematica"; % or paper
data_scheme.numericalSol = "gaussSeidel";
data_scheme.numSORw = 0.85;

%-------------------------------------------------------------------------%

data_system = struct;

% temperature
data_system.T = 298;                  % in kelvin

% external electric field in x-direction
data_system.E_ext = 5e5;              % in V/m

% model geometry
data_system.L = 15e-6;                % 1e-6 = micro meter
data_system.H = 30e-6;                % 1e-6 = micro meter
data_system.H_int = 20e-6;            % 1e-6 = micro meter

% no of grid points
data_system.N_y = 100;
data_system.N_z = 50;

% zeta-potential at bottom wall
data_system.Z_b = -24e-3;             % 1e-3 = milli volts

% zeta-potential at side wall
data_system.Z_s = -24e-3;             % 1e-3 = milli volts

% zeta-potential at interface
data_system.Z_int = -24e-3;           % 1e-3 = milli volts

%-------------------------------------------------------------------------%

data_fluid = struct;

% Assumption : fluid molecules have 1 cation and 1 anion each
% so, |z+| = |z-| = z
% and, n_inf+ = n_inf- = n_inf

% valency
data_fluid.z = 1;

% molarity of solution (in bulk) in mol/L
data_fluid.M_bulk = 1e-3;

% ionic conc (in bulk) in m^-3
data_fluid.n_inf = data_const.N_A * data_fluid.M_bulk * 1e3;

% dielectric constant
data_fluid.eps = 80;

% viscosity in Pa*s
data_fluid.mu = 1e-3;

%-------------------------------------------------------------------------%

% Un-comment code to create/edit files

data.const = data_const;
data.system = data_system;
data.fluid = data_fluid;
data.scheme = data_scheme;

% writestruct(data, "defaultData.xml")

%-------------------------------------------------------------------------%