% Prof Project under Dr Renganathan
% Paper - 3
% Paper : Deng 2020

% Aayush Bhakna
% CH22B008

% Simulation Results

clc;
clear variables;
close all;

font_size = 24;

%-------------------------------------------------------------------------%

% Changing values

% Changing Interface potential jump to 12e-3
setValue('system', 'zeta_int', 12e-3);

% % Changing Wall Zeta Potential to 12e-3
% setValue('system', 'zeta_wall', -12e-3);

% % Changing Interface Charge density to 1e-3
% setValue('system', 'q_int', 1e-3);

% % Changing Dielectric const of annular fluid to 10
% setValue('fluid_a', 'eps_a', 10);

% % Changing Viscosity of core fluid to 5e-3
% setValue('fluid_c', 'mu_c', 5e-3);

% For more changes, look at the dimData.xml file for variable names
% Look into defaultData.xml for default values

%-------------------------------------------------------------------------%

% Electric Potential Field

sol = analyticalSol();
xData1 = [flip(sol.phi_star); sol.phi_star];
yData1 = [flip(-1 .* sol.r_star); sol.r_star];

sol = numericalSol();
xData2 = [flip(sol.phi_star); sol.phi_star];
yData2 = [flip(-1 .* sol.r_star); sol.r_star];

figure(1)
hold on
plot(xData1, yData1, Color='blue', LineWidth=1.5, LineStyle='-')
plot(xData2, yData2, Color='red', LineWidth=1.5, LineStyle='--')
hold off
legend('Analytical', 'Numerical')
xlabel('Dimensionless Electric Potential')
ylabel('Dimensionless Pipe Cross-section')
fontsize(font_size, "points")


%-------------------------------------------------------------------------%

% Fluid Velocity Field

sol = analyticalSol();
xData1 = [flip(sol.u_star); sol.u_star];
yData1 = [flip(-1 .* sol.r_star); sol.r_star];

sol = numericalSol();
xData2 = [flip(sol.u_star); sol.u_star];
yData2 = [flip(-1 .* sol.r_star); sol.r_star];

figure(2)
hold on
plot(xData1, yData1, Color='blue', LineWidth=1.5, LineStyle='-')
plot(xData2, yData2, Color='red', LineWidth=1.5, LineStyle='--')
hold off
legend('Analytical', 'Numerical')
xlabel('Dimensionless Velocity Field')
ylabel('Dimensionless Pipe Cross-section')
fontsize(font_size, "points")

%-------------------------------------------------------------------------%