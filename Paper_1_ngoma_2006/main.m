% Prof Project under Dr Renganathan
% Paper - 1
% Paper : Ngoma 2006

% Aayush Bhakna
% CH22B008

% Simulation Results

clear variables;
clc;
close all;

%-------------------------------------------------------------------------%

% Reading stored data
data = readstruct("defaultData.xml");
const = readstruct("constants.xml");
data.BC_interface = "cheng"; % or "ngoma"
font_size = 24;

data.T = 296.15;
phiT = const.k_b * data.T / (data.z_nott * const.e);

%-------------------------------------------------------------------------%

% Water + CHB
data.alpha = 2.269/0.93;

% n_inf = 1 mM
data.n_inf = 6.022e23;
data.zeta1 = 0.11511;
data.zeta2 = -0.0511;
data.u_nott = - data.eps * const.eps_nott * data.E_ext * data.zeta1 / (data.mu);

u_1nM = data.u_nott;

sol = anaSolMod(data);

xData1 = phiT .* sol.phi_star;
yData1 = sol.y1_star;

zData1 = u_1nM .* sol.u_star;
wData1 = sol.y_star;

sol = numericalSol(data);

xData2 = phiT .* sol.phi_star;
yData2 = sol.y1_star;

zData2 = u_1nM .* sol.u_star;
wData2 = sol.y_star;

% n_inf = 10 mM
data.n_inf = 10 * 6.022e23;
data.zeta1 = 0.03809;
data.zeta2 = -0.020753;
data.u_nott = - data.eps * const.eps_nott * data.E_ext * data.zeta1 / (data.mu);

u_10nM = data.u_nott;

sol = anaSolMod(data);

xData3 = phiT .* sol.phi_star;
yData3 = sol.y1_star;

zData3 = u_10nM .* sol.u_star;
wData3 = sol.y_star;

sol = numericalSol(data);

xData4 = phiT .* sol.phi_star;
yData4 = sol.y1_star;

zData4 = u_10nM .* sol.u_star;
wData4 = sol.y_star;

figure(1)
hold on
scatter(xData1, yData1, 40, "blue", "filled", "o", ...
    DisplayName="Analytical (1 mM)")
plot(xData2, yData2, ...
    LineStyle="-", LineWidth=1.5, Color="blue", ...
    DisplayName="Numerical (1 mM)")
scatter(xData3, yData3, 40, "red", "filled", "o", ...
    DisplayName="Analytical (10 mM)")
plot(xData4, yData4, ...
    LineStyle="-", LineWidth=1.5, Color="red", ...
    DisplayName="Numerical (10 mM)")
hold off
grid on
legend(Location="best")
xlabel('Electric Potential ( \Phi )')
ylabel('Channel Width Fraction ( y^* )')
title('Water + CHB')
fontsize(font_size, "points")

% saveas(gcf, "ngoma/CHB_phi.png")

figure(2)
hold on
scatter(zData1, wData1, 40, "blue", "filled", "o", ...
    DisplayName="Analytical (1 mM)")
plot(zData2, wData2, ...
    LineStyle="-", LineWidth=1.5, Color="blue", ...
    DisplayName="Numerical (1 mM)")
scatter(zData3, wData3, 40, "red", "filled", "o", ...
    DisplayName="Analytical (10 mM)")
plot(zData4, wData4, ...
    LineStyle="-", LineWidth=1.5, Color="red", ...
    DisplayName="Numerical (10 mM)")
hold off
grid on
legend(Location="best")
xlabel('Fluid Velocity ( u )')
ylabel('Channel Width Fraction ( y^* )')
title('Water + CHB')
fontsize(font_size, "points")

% saveas(gcf, "ngoma/CHB_u.png")

%-------------------------------------------------------------------------%

% Water + Octane
data.alpha = 0.566/0.93;

% n_inf = 1 mM
data.n_inf = 6.022e23;
data.zeta1 = 0.11511;
data.zeta2 = -0.05819;
data.u_nott = data.eps * const.eps_nott * data.E_ext * const.k_b * data.T ...
    / (data.mu * data.z_nott * const.e);

u_1nM = data.u_nott;

sol = anaSolMod(data);

xData1 = phiT .* sol.phi_star;
yData1 = sol.y1_star;

zData1 = u_1nM .* sol.u_star;
wData1 = sol.y_star;

sol = numericalSol(data);

xData2 = phiT .* sol.phi_star;
yData2 = sol.y1_star;

zData2 = u_1nM .* sol.u_star;
wData2 = sol.y_star;

% n_inf = 10 mM
data.n_inf = 10 * 6.022e23;
data.zeta1 = 0.03809;
data.zeta2 = -0.023765;
data.u_nott = data.eps * const.eps_nott * data.E_ext * const.k_b * data.T ...
    / (data.mu * data.z_nott * const.e);

u_10nM = data.u_nott;

sol = anaSolMod(data);

xData3 = phiT .* sol.phi_star;
yData3 = sol.y1_star;

zData3 = u_10nM .* sol.u_star;
wData3 = sol.y_star;

sol = numericalSol(data);

xData4 = phiT .* sol.phi_star;
yData4 = sol.y1_star;

zData4 = u_10nM .* sol.u_star;
wData4 = sol.y_star;

figure(3)
hold on
scatter(xData1, yData1, 40, "blue", "filled", "o", ...
    DisplayName="Analytical (1 mM)")
plot(xData2, yData2, ...
    LineStyle="-", LineWidth=1.5, Color="blue", ...
    DisplayName="Numerical (1 mM)")
scatter(xData3, yData3, 40, "red", "filled", "o", ...
    DisplayName="Analytical (10 mM)")
plot(xData4, yData4, ...
    LineStyle="-", LineWidth=1.5, Color="red", ...
    DisplayName="Numerical (10 mM)")
hold off
grid on
legend(Location="best")
xlabel('Electric Potential ( \Phi )')
ylabel('Channel Width Fraction ( y^* )')
title('Water + Octane')
fontsize(font_size, "points")

% saveas(gcf, "ngoma/Octane_phi.png")

figure(4)
hold on
scatter(zData1, wData1, 40, "blue", "filled", "o", ...
    DisplayName="Analytical (1 mM)")
plot(zData2, wData2, ...
    LineStyle="-", LineWidth=1.5, Color="blue", ...
    DisplayName="Numerical (1 mM)")
scatter(zData3, wData3, 40, "red", "filled", "o", ...
    DisplayName="Analytical (10 mM)")
plot(zData4, wData4, ...
    LineStyle="-", LineWidth=1.5, Color="red", ...
    DisplayName="Numerical (10 mM)")
hold off
grid on
legend(Location="best")
xlabel('Fluid Velocity ( u )')
ylabel('Channel Width Fraction ( y^* )')
title('Water + Octane')
fontsize(font_size, "points")

% saveas(gcf, "ngoma/Octane_u.png")

%-------------------------------------------------------------------------%