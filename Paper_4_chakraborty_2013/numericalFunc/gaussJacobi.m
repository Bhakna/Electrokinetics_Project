function [ u_star ] = ...
    gaussJacobi(De, gamma_r, phi_star, y_star, z_star, u_int)

% grid points
N_y = length(y_star);
N_z = length(z_star); % rows
N_z = ((N_z - 1)/2) + 1;

% step size
del_y = y_star(2) - y_star(1);
del_z = z_star(2) - z_star(1);

% initializing mesh
% u_j_k -> u( k, j ) 
u_star = zeros(N_z, N_y);

% boundary conditions

% bottom wall
u_star(:, 1) = 0;

% right side wall
u_star(end, :) = 0;

% interface and left side
if u_int == "NONE"
    interfaceNeumann = "YES";
else
    u_star(1, :) = 0;
    u_star(:, end) = u_int;
    interfaceNeumann = "NO";
end

% while loop
u_old = u_star;
iter = 0;
RMSE = 1;
maxIter = 9999;
tol = 1e-9;
while ((iter < maxIter) & (RMSE > tol))

    iter = iter + 1;

    u_star = u_next(u_old, N_y, N_z, del_y, del_z, De, gamma_r, phi_star, interfaceNeumann);

    err = (u_star - u_old);
    RMSE = sqrt(mean(abs(err), "all"));

    u_old = u_star;

end

u_star_new = zeros(length(z_star), N_y); % actual domain
for k = 1:(N_z-1)
    u_star_new(k, :) = u_star(N_z+1-k, :);
end
u_star_new(N_z:end, :) = u_star;

u_star = u_star_new;

end

function [ u_new ] = u_next(u_old, N_y, N_z, del_y, del_z, ...
    De, gamma_r, phi_star, interfaceNeumann)

u_new = u_old;

dn = (2 * (del_y^2)) + (2 * (del_z^2));
cj = (del_z^2) / dn;
ck = (del_y^2) / dn;
cphi = (-1 * gamma_r * ((del_y * del_z / De)^2)) / dn;

% bulk
for j = 2:N_y-1
    for k = 2:N_z-1

        u_new(k, j) = (cj * u_old(k, j+1)) + ...
            (cj * u_old(k, j-1)) + ...
            (ck * u_old(k+1, j)) + ...
            (ck * u_old(k-1, j)) + ...
            (cphi * phi_star(k, j));

    end
end

if interfaceNeumann == "YES"

% interface
j = N_y;
for k = 2:N_z-1
    
    del_phi = phi_star(k, j) - phi_star(k, j-1);
    u_new(k, j) = u_new(k, j-1) + (gamma_r * del_phi);

end

% left side
k = 1;
for j = 2:N_y-1

    u_new(k, j) = u_new(k+1, j);

end

end

end