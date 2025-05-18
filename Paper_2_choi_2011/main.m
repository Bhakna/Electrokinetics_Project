% Prof Project under Dr Renganathan
% Paper - 2
% Paper : Choi 2011

% Aayush Bhakna
% CH22B008

% Simulation Results

clc;
clear variables;
close all;

data = readstruct("dimData.xml");
phiT = data.const.k_b * data.system.T / data.const.e;

font_size = 24;

%-------------------------------------------------------------------------%

% Water + CHB
setValue('fluid1', 'eps_r1', 79)
setValue('fluid1', 'mu_1', 0.93e-3)
setValue('fluid2', 'eps_r2', 7.9)
setValue('fluid2', 'mu_2', 2.269e-3)

% 1 nM concentration
setValue('fluid1', 'n1_inf', 01*6.022e23)
setValue('fluid2', 'n2_inf', 01*6.022e23)
setValue('system', 'q_int', -4.8e-3)

sol = numModified();
xData3 = phiT .* sol.phi_star();
column1 = [xData3(1), xData3(200), xData3(end)]';

setValue('fluid1', 'zeta_1', xData3(1))
setValue('fluid2', 'zeta_2', xData3(end))

% plotting

sol = analyticalSol();
xData1 = phiT .* sol.phi_star;
yData1 = sol.y_star;

sol = numericalSol();
xData2 = phiT .* sol.phi_star;
yData2 = sol.y_star;

figure(1)
hold on
scatter(xData1, yData1, 40, "red", "filled", "o", ...
    DisplayName='Analytical (1 mM)')
plot(xData2, yData2, Color='red', ...
    LineWidth=1.5, LineStyle='-', DisplayName='Numerical (1 mM)')
hold off
legend(Location="best")
grid on
xlabel('Electric Potential  ( \Phi )')
ylabel('Channel Width Fraction ( y^* )')
title('Water + CHB')
fontsize(font_size, "points")

u_h1 = -1 * data.fluid1.eps_r1 * data.const.eps_nott * ...
    data.system.E_ext * data.fluid1.zeta_1 / data.fluid1.mu_1;

u_h2 = -1 * data.fluid2.eps_r2 * data.const.eps_nott * ...
    data.system.E_ext * data.fluid2.zeta_2 / data.fluid2.mu_2;

sol = analyticalSol();
xData1 = [(u_h1 .* sol.u_1_star); (u_h2 .* sol.u_2_star)];

sol = numericalSol();
xData2 = [(u_h1 .* sol.u_1_star); (u_h2 .* sol.u_2_star)];

figure(2)
hold on
scatter(xData1, yData1, 40, "red", "filled", "o", ...
    DisplayName='Analytical (1 mM)')
plot(xData2, yData2, Color='red', ...
    LineWidth=1.5, LineStyle='-', DisplayName='Numerical (1 mM)')
hold off
grid on
legend(Location="best")
xlabel('Velocity Field ( u )')
ylabel('Channel Width Fraction ( y^* )')
title('Water + CHB')
fontsize(font_size, "points")

% 10 nM concentration
setValue('fluid1', 'n1_inf', 10*6.022e23)
setValue('fluid2', 'n2_inf', 10*6.022e23)
setValue('system', 'q_int', -5.9e-3)

sol = numModified();
xData3 = phiT .* sol.phi_star();
column2 = [xData3(1), xData3(200), xData3(end)]';

setValue('fluid1', 'zeta_1', xData3(1))
setValue('fluid2', 'zeta_2', xData3(end))

% plotting

sol = analyticalSol();
xData1 = phiT .* sol.phi_star;
yData1 = sol.y_star;

sol = numericalSol();
xData2 = phiT .* sol.phi_star;
yData2 = sol.y_star;

figure(1)
hold on
scatter(xData1, yData1, 30, "blue", "filled", "o", ...
    DisplayName='Analytical (10 mM)')
plot(xData2, yData2, Color='blue', ...
    LineWidth=1.5, LineStyle='-', DisplayName='Numerical (10 mM)')
hold off
legend(Location="best")
grid on
xlabel('Electric Potential  ( \Phi )')
ylabel('Channel Width Fraction ( y^* )')
title('Water + CHB')
fontsize(font_size, "points")

saveas(gcf, "choi/CHB_phi.png")

u_h1 = -1 * data.fluid1.eps_r1 * data.const.eps_nott * ...
    data.system.E_ext * data.fluid1.zeta_1 / data.fluid1.mu_1;

u_h2 = -1 * data.fluid2.eps_r2 * data.const.eps_nott * ...
    data.system.E_ext * data.fluid2.zeta_2 / data.fluid2.mu_2;

sol = analyticalSol();
xData1 = [(u_h1 .* sol.u_1_star); (u_h2 .* sol.u_2_star)];

sol = numericalSol();
xData2 = [(u_h1 .* sol.u_1_star); (u_h2 .* sol.u_2_star)];

figure(2)
hold on
scatter(xData1, yData1, 30, "blue", "filled", "o", ...
    DisplayName='Analytical (10 mM)')
plot(xData2, yData2, Color='blue', ...
    LineWidth=1.5, LineStyle='-', DisplayName='Numerical (10 mM)')
hold off
grid on
legend(Location="best")
xlabel('Velocity Field ( u )')
ylabel('Channel Width Fraction ( y^* )')
title('Water + CHB')
fontsize(font_size, "points")

saveas(gcf, "choi/CHB_u.png")

%-------------------------------------------------------------------------%

% Water + Octane
setValue('fluid1', 'eps_r1', 79)
setValue('fluid1', 'mu_1', 0.93e-3)
setValue('fluid2', 'eps_r2', 2.0)
setValue('fluid2', 'mu_2', 0.566e-3)

% 1 nM concentration
setValue('fluid1', 'n1_inf', 01*6.022e23)
setValue('fluid2', 'n2_inf', 01*6.022e23)
setValue('system', 'q_int', -4.8e-3)

sol = numModified();
xData3 = phiT .* sol.phi_star();
column3 = [xData3(1), xData3(200), xData3(end)]';

setValue('fluid1', 'zeta_1', xData3(1))
setValue('fluid2', 'zeta_2', xData3(end))

% plotting

sol = analyticalSol();
xData1 = phiT .* sol.phi_star;
yData1 = sol.y_star;

sol = numericalSol();
xData2 = phiT .* sol.phi_star;
yData2 = sol.y_star;

figure(3)
hold on
scatter(xData1, yData1, 40, "red", "filled", "o", ...
    DisplayName='Analytical (1 mM)')
plot(xData2, yData2, Color='red', ...
    LineWidth=1.5, LineStyle='-', DisplayName='Numerical (1 mM)')
hold off
legend(Location="best")
grid on
xlabel('Electric Potential  ( \Phi )')
ylabel('Channel Width Fraction ( y^* )')
title('Water + Octane')
fontsize(font_size, "points")

u_h1 = -1 * data.fluid1.eps_r1 * data.const.eps_nott * ...
    data.system.E_ext * data.fluid1.zeta_1 / data.fluid1.mu_1;

u_h2 = -1 * data.fluid2.eps_r2 * data.const.eps_nott * ...
    data.system.E_ext * data.fluid2.zeta_2 / data.fluid2.mu_2;

sol = analyticalSol();
xData1 = [(u_h1 .* sol.u_1_star); (u_h2 .* sol.u_2_star)];

sol = numericalSol();
xData2 = [(u_h1 .* sol.u_1_star); (u_h2 .* sol.u_2_star)];

figure(4)
hold on
scatter(xData1, yData1, 40, "red", "filled", "o", ...
    DisplayName='Analytical (1 mM)')
plot(xData2, yData2, Color='red', ...
    LineWidth=1.5, LineStyle='-', DisplayName='Numerical (1 mM)')
hold off
grid on
legend(Location="best")
xlabel('Velocity Field ( u )')
ylabel('Channel Width Fraction ( y^* )')
title('Water + Octane')
fontsize(font_size, "points")

% 10 nM concentration
setValue('fluid1', 'n1_inf', 10*6.022e23)
setValue('fluid2', 'n2_inf', 10*6.022e23)
setValue('system', 'q_int', -5.9e-3)

sol = numModified();
xData3 = phiT .* sol.phi_star();
column4 = [xData3(1), xData3(200), xData3(end)]';

setValue('fluid1', 'zeta_1', xData3(1))
setValue('fluid2', 'zeta_2', xData3(end))

% plotting

sol = analyticalSol();
xData1 = phiT .* sol.phi_star;
yData1 = sol.y_star;

sol = numericalSol();
xData2 = phiT .* sol.phi_star;
yData2 = sol.y_star;

figure(3)
hold on
scatter(xData1, yData1, 30, "blue", "filled", "o", ...
    DisplayName='Analytical (10 mM)')
plot(xData2, yData2, Color='blue', ...
    LineWidth=1.5, LineStyle='-', DisplayName='Numerical (10 mM)')
hold off
legend(Location="best")
grid on
xlabel('Electric Potential  ( \Phi )')
ylabel('Channel Width Fraction ( y^* )')
title('Water + Octane')
fontsize(font_size, "points")

saveas(gcf, "choi/Octane_phi.png")

u_h1 = -1 * data.fluid1.eps_r1 * data.const.eps_nott * ...
    data.system.E_ext * data.fluid1.zeta_1 / data.fluid1.mu_1;

u_h2 = -1 * data.fluid2.eps_r2 * data.const.eps_nott * ...
    data.system.E_ext * data.fluid2.zeta_2 / data.fluid2.mu_2;

sol = analyticalSol();
xData1 = [(u_h1 .* sol.u_1_star); (u_h2 .* sol.u_2_star)];

sol = numericalSol();
xData2 = [(u_h1 .* sol.u_1_star); (u_h2 .* sol.u_2_star)];

figure(4)
hold on
scatter(xData1, yData1, 30, "blue", "filled", "o", ...
    DisplayName='Analytical (10 mM)')
plot(xData2, yData2, Color='blue', ...
    LineWidth=1.5, LineStyle='-', DisplayName='Numerical (10 mM)')
hold off
grid on
legend(Location="best")
xlabel('Velocity Field ( u )')
ylabel('Channel Width Fraction ( y^* )')
title('Water + Octane')
fontsize(font_size, "points")

saveas(gcf, "choi/Octane_u.png")

%-------------------------------------------------------------------------%

% Table

output = table;
output.Variable = ["Phi Water/Wall"; "Phi Water/Oil"; "Phi Oil/Wall"];
output.CHB_1 = column1;
output.CHB_10 = column2;
output.Oct_1 = column3;
output.Oct_10 = column4;

disp(output)

%-------------------------------------------------------------------------%