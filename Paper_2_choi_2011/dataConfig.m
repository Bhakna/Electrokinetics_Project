%-------------------------------------------------------------------------%

clc;
clear variables;

data_const = struct;

data_const.eps_nott = 8.854e-12;
data_const.e = 1.6022e-19;
data_const.k_b = 1.3807e-23;
data_const.N_A = 6.022e23;

%-------------------------------------------------------------------------%

data_system = struct;

data_system.T = 296.15;               % in kelvin
data_system.h = 100e-9;               % in meter
data_system.E_ext = 5e5;              % in V/m
data_system.zeta_int = 0;             % in volts
data_system.q_int = -4.8e-3;          % in C/m^2

%-------------------------------------------------------------------------%

% Bottom fluid (water)

data_fluid1 = struct;

data_fluid1.z1 = 1;                   % valence
data_fluid1.n1_inf = data_const.N_A;  % in m^-3
data_fluid1.eps_r1 = 79;              % dielectric const
data_fluid1.h1 = 80e-9;               % in meters

data_fluid1.zeta_1 = 3.8090e-02;          % in volts
data_fluid1.q_1 = -8.2e-3;            % in C/m^2

data_fluid1.mu_1 = 0.93e-3;           % in Pa*s
data_fluid1.N1 = 200;                 % grid points

%-------------------------------------------------------------------------%

% Top fluid (CHB)

data_fluid2 = struct;

data_fluid2.z2 = 1;                   % valence
data_fluid2.n2_inf = data_const.N_A;  % in m^-3
data_fluid2.eps_r2 = 2.0;             % dielectric const

data_fluid2.zeta_2 = -1.4490e-20;      % in volts
data_fluid2.q_2 = -0.21e-3;            % in C/m^2

data_fluid2.mu_2 = 0.566e-3;          % in Pa*s
data_fluid2.N2 = 200;                 % grid points

data_fluid2.h2 = data_system.h - data_fluid1.h1;

%-------------------------------------------------------------------------%

% Un-comment code to create/edit files

data.const = data_const;
data.system = data_system;
data.fluid1 = data_fluid1;
data.fluid2 = data_fluid2;

writestruct(data, "defaultData.xml")
writestruct(data, "dimData.xml")

%-------------------------------------------------------------------------%