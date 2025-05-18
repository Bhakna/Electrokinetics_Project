function [ u_star ] = iterativeAxb(De, gamma_r, phi_star, y_star, z_star)

% grid points
N_y = length(y_star);
N_z = length(z_star); % rows
N_yz = N_y * N_z;

% step size
del_y = y_star(2) - y_star(1);
del_z = z_star(2) - z_star(1);

% x -> reshape(y_z, y*z, 1)

% calculating coeffs in A
A = zeros(N_yz, N_yz);
b = zeros(N_yz, 1);

% bulk
for j = 2:N_y-1
    for k = 2:N_z-1
    
        A(mapping(k, j, N_z), mapping(k, j, N_z)) = ...
            -2 * ( (1/del_y^2) + (1/del_z^2) );

        A(mapping(k, j, N_z), mapping(k+1, j, N_z)) = 1/del_z^2;
        A(mapping(k, j, N_z), mapping(k-1, j, N_z)) = 1/del_z^2;

        A(mapping(k, j, N_z), mapping(k, j+1, N_z)) = 1/del_y^2;
        A(mapping(k, j, N_z), mapping(k, j-1, N_z)) = 1/del_y^2;

        b(mapping(k, j, N_z)) = (gamma_r/(De^2)) * phi_star(k, j);

    end
end

% bottom plate
j = 1;
for k = 2:N_z-1    
    A(mapping(k, j, N_z), mapping(k, j, N_z)) = 1;
end

% left plate
k = 1;
for j = 1:N_y    
    A(mapping(k, j, N_z), mapping(k, j, N_z)) = 1;
end

% right plate
k = N_z;
for j = 1:N_y    
    A(mapping(k, j, N_z), mapping(k, j, N_z)) = 1;
end

% interface
j = N_y;
for k = 2:N_z-1
    
    A(mapping(k, j, N_z), mapping(k, j, N_z)) = 1;
    A(mapping(k, j, N_z), mapping(k, j-1, N_z)) = -1;

    b(mapping(k, j, N_z)) = gamma_r * (phi_star(k, j) - phi_star(k, j-1));

end

% solution
x = lsqminnorm(A, b);

% reverse mapping
u_star = reshape(x, N_z, N_y);

end

function [ output ] = mapping(i, j, rows)

output = (rows * (j-1)) + i;

end