%-------------------------------------------------------------------------%

clc;
clear variables;

data_const = struct;

data_const.eps_o = 8.854e-12;
data_const.e = 1.6022e-19;
data_const.k_b = 1.3807e-23;

%-------------------------------------------------------------------------%

data_system = struct;

data_system.T = 298;                  % in kelvin
data_system.R = 35e-6;                % 1e-6 = micro meter
data_system.r_o = 17e-6;              % 1e-6 = micro meter
data_system.E_ext = 5e5;              % in V/m
data_system.zeta_wall = -24e-3;       % in volts
data_system.zeta_int = 0;             % 1e-3 = milli volts
data_system.q_int = 0;                % in C/m^2

%-------------------------------------------------------------------------%

data_fluidc = struct;

data_fluidc.z_c = 1;                   % valence
data_fluidc.n_c_inf = 6.022e20;        % in m^-3
data_fluidc.eps_c = 80;                % dielectric const
data_fluidc.mu_c = 1e-3;               % in Pa*s
data_fluidc.Nc = 200;                  % grid points

%-------------------------------------------------------------------------%

data_fluida = struct;

data_fluida.z_a = 1;                   % valence
data_fluida.n_a_inf = 6.022e20;        % in m^-3
data_fluida.eps_a = 80;                % dielectric const
data_fluida.mu_a = 1e-3;               % in Pa*s
data_fluida.Na = 200;                  % grid points

%-------------------------------------------------------------------------%

% Un-comment code to create/edit files

data.const = data_const;
data.system = data_system;
data.fluid_c = data_fluidc;
data.fluid_a = data_fluida;

% writestruct(data, "defaultData.xml")

%-------------------------------------------------------------------------%