function [ sol ] = numericalSol()

data = readstruct("nondimData.xml");

%-------------------------------------------------------------------------%

% GENERATING RESULTS

% 1-D Mesh/Grid
sol.y_star = linspace(0, data.h1_star, data.N_y);
sol.z_star = linspace(-data.l_star, data.l_star, data.N_z);
[sol.Y_mesh, sol.Z_mesh] = meshgrid(sol.y_star, sol.z_star);

% Electric Potential
[ sol.phi_star ] = ...
    electricPotential(data.De, ...
    data.Z_b_star, data.Z_s_star, data.Z_int_star, ...
    sol.y_star, sol.z_star);

% Velocity Field

if data.numSolFunc == "gaussJacobi"

[ sol.u_star ] = ...
    gaussJacobi(data.De, data.gamma_r, ...
    sol.phi_star, sol.y_star, sol.z_star, ...
    "NONE");

elseif data.numSolFunc == "gaussSeidel"

[ sol.u_star ] = ...
    gaussSeidel(data.De, data.gamma_r, ...
    sol.phi_star, sol.y_star, sol.z_star, ...
    "NONE");

elseif data.numSolFunc == "SOR"

[ sol.u_star ] = ...
    SOR(data.De, data.gamma_r, ...
    sol.phi_star, sol.y_star, sol.z_star, ...
    "NONE", data.numSORw);

elseif data.numSolFunc == "iterativeGaussSeidel"

[ sol.u_int ] = ...
    iterativeGaussSeidel(data.De, data.gamma_r, ...
    sol.phi_star, sol.y_star, sol.z_star);
[ sol.u_star ] = ...
    gaussSeidel(data.De, data.gamma_r, ...
    sol.phi_star, sol.y_star, sol.z_star, ...
    sol.u_int);

elseif data.numSolFunc == "directAxb"

[ sol.u_star ] = ...
    directAxb(data.De, data.gamma_r, ...
    sol.phi_star, sol.y_star, sol.z_star);

elseif data.numSolFunc == "iterativeAxb"

[ sol.u_star ] = ...
    iterativeAxb(data.De, data.gamma_r, ...
    sol.phi_star, sol.y_star, sol.z_star);

else
    error("Please enter correct numerical scheme...")
end

end

%-------------------------------------------------------------------------%

% ELECTRIC POTENTIAL

function [ phi_star ] = ...
    electricPotential(De, Z_b_star, Z_s_star, Z_int_star, y_star, z_star)

% grid points
N_y = length(y_star);
N_z = length(z_star);

% step size
del_y = y_star(2) - y_star(1);
del_z = z_star(2) - z_star(1);

% initializing mesh
% phi_j_k -> phi( k, j ) 
phi_star = zeros(N_z, N_y);

% boundary conditions

% bottom wall
phi_star(:, 1) = Z_b_star;

% side walls
phi_star(1, :) = Z_s_star;
phi_star(end, :) = Z_s_star;

% interface
phi_star(:, end) = Z_int_star;

% while loop
phi_old = phi_star;
iter = 0;
RMSE = 1;
maxIter = 9999;
tol = 1e-9;
while ((iter < maxIter) & (RMSE > tol))

    iter = iter + 1;

    phi_star = phi_next(phi_old, N_y, N_z, del_y, del_z, De);

    err = (phi_star - phi_old);
    RMSE = sqrt(mean(abs(err), "all"));

    phi_old = phi_star;

end

end

function [ phi_new ] = phi_next(phi_old, N_y, N_z, del_y, del_z, De)

phi_new = phi_old;

dn = (2 * (del_y^2)) + (2 * (del_z^2)) + ((del_z * del_y / De)^2);
cj = (del_z^2) / dn;
ck = (del_y^2) / dn;

for j = 2:N_y-1
    for k = 2:N_z-1

        phi_new(k, j) = (cj * phi_old(k, j+1)) + ...
            (cj * phi_new(k, j-1)) + ...
            (ck * phi_old(k+1, j)) + ...
            (ck * phi_new(k-1, j));

    end
end

end

%-------------------------------------------------------------------------%