function [ u_star ] = ...
    velocityPaper(gamma_r, h1_star, l_star, ...
    Z_b_star, Z_s_star, ...
    phi_star, Y_mesh, Z_mesh)

old_u_m = zeros(size(Y_mesh));

maxIter = 999;
tol = 1e-6;

iter = 0;
RMSE = 1;

while ((iter < maxIter) & (RMSE > tol))

    iter = iter + 1;

    u_m = old_u_m + u_mth(iter, gamma_r, h1_star, l_star, ...
    Z_b_star, Y_mesh, Z_mesh);

    err = (u_m - old_u_m);
    RMSE = sqrt(mean(abs(err), "all"));

    old_u_m = u_m;

end

fprintf('\n max m = %d \n', iter)
fprintf('\n RMSE m = %0.4e \n', RMSE)

old_u_n = zeros(size(Y_mesh));

iter = 0;
RMSE = 1;

while ((iter < maxIter) & (RMSE > tol))

    iter = iter + 1;

    u_n = old_u_n + u_nth(iter, gamma_r, h1_star, l_star, ...
    Z_s_star, Y_mesh, Z_mesh);

    err = (u_n - old_u_n);
    RMSE = sqrt(mean(abs(err), "all"));

    old_u_n = u_n;

end

fprintf('\n max n = %d \n', iter)
fprintf('\n RMSE n = %0.4e \n', RMSE)

u_star = (u_m - u_n) + (gamma_r .* phi_star);

end

function [ output ] = ...
    u_mth(m, gamma_r, h1_star, l_star, ...
    Z_b_star, Y_mesh, Z_mesh)

% lambda_m -> f_m
f_m = pi * ((2*m) - 1) / (2 * l_star);

c1 = 2 * ((-1)^m) * gamma_r;
c2 = (Z_b_star .* cosh(f_m .* (h1_star - Y_mesh))) ...
    ./ ( l_star * f_m * cosh(f_m * h1_star) );
c3 = cos(f_m .* Z_mesh);

output = c1 .* c2 .* c3;

end

function [ output ] = ...
    u_nth(n, gamma_r, h1_star, l_star, ...
    Z_s_star, Y_mesh, Z_mesh)

% lambda_n -> f_n
f_n = pi * n / h1_star; % Chakraborty

c1 = 2 * (1 + ((-1)^(n-1))) * gamma_r; % Chakraborty
c2 = (Z_s_star .* cosh(f_n .* Z_mesh)) ...
    ./ ( h1_star * f_n * cosh(f_n * l_star) );
c3 = sin(f_n .* Y_mesh);

output = c1 .* c2 .* c3;

end