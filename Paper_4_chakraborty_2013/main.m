% Prof Project under Dr Renganathan
% Paper - 4
% Paper : Mayur & Chakraborty 2013

% Aayush Bhakna
% CH22B008

% Simulation Results

clc;
clear variables;
close all;

addpath analyticalFunc
addpath numericalFunc

setValue("default", "", "", "");
font_size = 20;

%-------------------------------------------------------------------------%

% Changing values

% Fig 4

N = 25;
N = (2*N) + 1;

% Changing N_y = 200
setValue("dim", "system", "N_y", N);

% Changing N_z = 200
setValue("dim", "system", "N_z", N);

% Setting schemes
setValue("dim", "scheme", "analyticalSol", "mathematica");
setValue("dim", "scheme", "numericalSol", "gaussJacobi");
setValue("dim", "scheme", "numSORw", 1.0);

% Changing temp
setValue("dim", "system", "T", 296.15);

% Changing geometry
setValue("dim", "system", "L", 100e-9);
setValue("dim", "system", "H", 100e-9);
setValue("dim", "system", "H_int", 80e-9);

% c = 10 mM
setValue("dim", "fluid", "M_bulk", 10e-3)

% water
setValue("dim", "fluid", "eps", 79)
setValue("dim", "fluid", "mu", 0.93e-3)

% Changing zeta-values (10 mM)
Z_b = 0.03809;
setValue("dim", "system", "Z_b", Z_b);
setValue("dim", "system", "Z_s", Z_b);
setValue("dim", "system", "Z_int", -0.020753);

% For more changes, look at the dimData.xml file for variable names
% Look into defaultData.xml for default values

%-------------------------------------------------------------------------%

% Electric Potential Field

data = readstruct("nondimData.xml");

sol = analyticalSol();
colorData1 = Z_b .* sol.phi_star';
yData1 = sol.Y_mesh';
zData1 = sol.Z_mesh';
Y1 = yData1 ./ data.h1_star;
Z1 = zData1 ./ data.l_star;

del_phi_1 = sol.phi_star(2:end-1, end) - sol.phi_star(2:end-1, end-1);

sol = numericalSol();
colorData2 = Z_b .* sol.phi_star';
yData2 = sol.Y_mesh';
zData2 = sol.Z_mesh';
Y2 = yData2 ./ data.h1_star;
Z2 = zData2 ./ data.l_star;

del_phi_2 = sol.phi_star(2:end-1, end) - sol.phi_star(2:end-1, end-1);

figure(Name="Electric Potential")

subplot(2, 2, 1)
contourf(Z1, Y1, colorData1, 'ShowText','on')
colormap("default")
colorbar
xlabel('z^* / l^*')
ylabel('y^* / h_1^*')
title('Analytical Solution')
fontsize(font_size, "points")

subplot(2, 2, 2)
contourf(Z2, Y2, colorData2, 'ShowText','on')
colormap("default")
colorbar
xlabel('z^* / l^*')
ylabel('y^* / h_1^*')
title('Numerical Solution')
fontsize(font_size, "points")

subplot(2, 2, 3)
hold on
y_index = round(data.N_y / 2);
plot(Z1(y_index, :), colorData1(y_index, :), ...
    LineWidth=1.75, LineStyle="-", DisplayName="Analytical")
plot(Z2(y_index, :), colorData2(y_index, :), ...
    LineWidth=1.75, LineStyle="--", DisplayName="Numerical")
hold off
grid on
legend;
xlabel('z^* / l^*')
ylabel('\Phi')
title('Solution at ( y / h_1 = 0.5 ) ')
fontsize(font_size, "points")

subplot(2, 2, 4)
hold on
z_index = round(data.N_z / 2);
plot(Y1(:, z_index), colorData1(:, z_index), ...
    LineWidth=1.75, LineStyle="-", DisplayName="Analytical")
plot(Y2(:, z_index), colorData2(:, z_index), ...
    LineWidth=1.75, LineStyle="--", DisplayName="Numerical")
hold off
grid on
legend;
xlabel('y^* / h_1^*')
ylabel('\Phi')
title('Solution at ( z = 0 ) ')
fontsize(font_size, "points")

%-------------------------------------------------------------------------%

% Fluid Velocity Field

sol = analyticalSol();
colorData1 = data.u_ref .* sol.u_star';

del_u_1 = sol.u_star(2:end-1, end) - sol.u_star(2:end-1, end-1);

sol = numericalSol();
colorData2 = data.u_ref .* sol.u_star';

del_u_2 = sol.u_star(2:end-1, end) - sol.u_star(2:end-1, end-1);

figure(Name="Fluid Velocity")

subplot(2, 2, 1)
contourf(Z1, Y1, colorData1, 10, 'ShowText','on')
colormap("default")
colorbar
xlabel('z^* / l^*')
ylabel('y^* / h_1^*')
title('Analytical Solution')
fontsize(font_size, "points")

subplot(2, 2, 2)
contourf(Z2, Y2, colorData2, 10, 'ShowText','on')
colormap("default")
colorbar
xlabel('z^* / l^*')
ylabel('y^* / h_1^*')
title('Numerical Solution')
fontsize(font_size, "points")

subplot(2, 2, 3)
hold on
y_index = round(data.N_y / 2);
plot(Z1(y_index, :), colorData1(y_index, :), ...
    LineWidth=1.75, LineStyle="-", DisplayName="Analytical")
plot(Z2(y_index, :), colorData2(y_index, :), ...
    LineWidth=1.75, LineStyle="--", DisplayName="Numerical")
hold off
grid on
legend;
xlabel('z^* / l^*')
ylabel('u')
title('Solution at ( y / h_1 = 0.5 ) ')
fontsize(font_size, "points")

subplot(2, 2, 4)
hold on
z_index = round(data.N_z / 2);
plot(Y1(:, z_index), colorData1(:, z_index), ...
    LineWidth=1.75, LineStyle="-", DisplayName="Analytical")
plot(Y2(:, z_index), colorData2(:, z_index), ...
    LineWidth=1.75, LineStyle="--", DisplayName="Numerical")
hold off
grid on
legend;
xlabel('y^* / h_1^*')
ylabel('u')
title('Solution at ( z = 0 ) ')
fontsize(font_size, "points")

err1 = del_u_1 - (data.gamma_r .* del_phi_1);
RMSE1 = sqrt(mean(err1.^2, "all"));

err2 = del_u_2 - (data.gamma_r .* del_phi_2);
RMSE2 = sqrt(mean(err2.^2, "all"));


% figure(Name="3D Fluid Velocity")
% 
% subplot(1, 2, 1)
% surf(Z1, Y1, colorData1)
% colormap("default")
% colorbar
% xlabel('z^* / l^*')
% ylabel('y^* / h_1^*')
% title('Analytical Solution')
% fontsize(font_size, "points")
% 
% subplot(1, 2, 2)
% surf(Z2, Y2, colorData2)
% colormap("default")
% colorbar
% xlabel('z^* / l^*')
% ylabel('y^* / h_1^*')
% title('Numerical Solution')
% fontsize(font_size, "points")

figure(Name="3D Fluid Velocity")

subplot(1, 2, 1)
surf(Z1(2:end-1, 2:end-1), Y1(2:end-1, 2:end-1), ...
    colorData1(2:end-1, 2:end-1))
colormap("default")
colorbar
xlabel('z^* / l^*')
ylabel('y^* / h_1^*')
title('Analytical Solution')
fontsize(font_size, "points")

subplot(1, 2, 2)
surf(Z2(2:end-1, 2:end-1), Y2(2:end-1, 2:end-1), ...
    colorData2(2:end-1, 2:end-1))
colormap("default")
colorbar
xlabel('z^* / l^*')
ylabel('y^* / h_1^*')
title('Numerical Solution')
fontsize(font_size, "points")

%-------------------------------------------------------------------------%