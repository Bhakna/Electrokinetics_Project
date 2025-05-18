function [ phi_star ] = ...
    electricPotential(De, h1_star, l_star, ...
    Z_b_star, Z_s_star, Z_int_star, ...
    Y_mesh, Z_mesh)

old_phi_m = zeros(size(Y_mesh));

iter = 0;

maxIter = 500;
tol = 1e-9;

while (iter < maxIter)

    iter = iter + 1;

    phi_m = old_phi_m + phi_mth(iter, De, h1_star, l_star, ...
    Z_b_star, Z_int_star, ...
    Y_mesh, Z_mesh);

    err = (phi_m - old_phi_m);
    RMSE = sqrt(mean(abs(err), "all"));

    if (RMSE <= tol) || (isnan(RMSE))
        phi_m = old_phi_m;
        break;
    end

    old_phi_m = phi_m;

end

fprintf('Phi m = %d \n', iter)

old_phi_n = zeros(size(Y_mesh));

iter = -1;

while (iter < maxIter)

    iter = iter + 2;

    phi_n = old_phi_n + phi_nth(iter, De, h1_star, l_star, ...
    Z_s_star, Y_mesh, Z_mesh);

    err = (phi_n - old_phi_n);
    RMSE = sqrt(mean(abs(err), "all"));

    if (RMSE <= tol) || (isnan(RMSE))
        phi_n = old_phi_n;
        break;
    end

    old_phi_n = phi_n;

end

fprintf('Phi n = %d \n', iter)

phi_star = phi_m + phi_n;

end

function [ output ] = ...
    phi_mth(m, De, h1_star, l_star, ...
    Z_b_star, Z_int_star, ...
    Y_mesh, Z_mesh)

% lambda_m -> f_m
f_m = pi * ((2*m) - 1) / (2 * l_star);

c0 = sqrt((1/(De^2)) + ((f_m^2)));
c1 = 2 * ((-1)^(m-1));
c2 = ( (Z_int_star .* sinh(c0 .* Y_mesh)) + ...
    (Z_b_star .* sinh(c0 .* (h1_star - Y_mesh))) ) ...
    ./ ( l_star * f_m * sinh(c0 * h1_star) );
c3 = cos(f_m .* Z_mesh);

output = c1 .* c2 .* c3;

end

function [ output ] = ...
    phi_nth(n, De, h1_star, l_star, ...
    Z_s_star, Y_mesh, Z_mesh)

% lambda_n -> f_n
f_n = pi * n / h1_star;

c0 = sqrt((1/(De^2)) + ((f_n^2)));
c1 = 2 * (1 + ((-1)^(n-1)));
c2 = (Z_s_star .* cosh(c0 .* Z_mesh)) ...
    ./ ( h1_star * f_n * cosh(c0 * l_star) );
c3 = sin(f_n .* Y_mesh);

output = c1 .* c2 .* c3;

end