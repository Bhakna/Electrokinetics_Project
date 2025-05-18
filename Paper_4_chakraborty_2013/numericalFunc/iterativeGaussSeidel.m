function [ u_int_val ] = iterativeGaussSeidel(De, gamma_r, phi_star, y_star, z_star)

u_int_guess = 0;

u_int_val = fminunc(@(u_int) errFunc(De, gamma_r, phi_star, y_star, z_star, u_int), ...
    u_int_guess);

end

function [ RMSE ] = errFunc(De, gamma_r, phi_star, y_star, z_star, u_int)

% grid points
N_y = length(y_star);
N_z = length(z_star);

u_star = gaussSeidel(De, gamma_r, phi_star, y_star, z_star, u_int);

err = zeros(N_z);
j = N_y;
for k = 2:N_z-1

    del_phi = phi_star(k, j) - phi_star(k, j-1);
    del_u = u_star(k, j) - u_star(k, j-1);
    err(k) = del_u - (gamma_r * del_phi);

end

RMSE = sqrt(mean(err.^2, "all"));

end

