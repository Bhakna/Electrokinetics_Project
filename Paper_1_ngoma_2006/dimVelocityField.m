function [ u_star ] = dimVelocityField(K, C0, C1, C2, alpha, h1_star, N, phi_star, BC_interface)

%-------------------------------------------------------------------------%

% CREATING MESH FOR SYSTEM

% Fluid 1
ymesh1 = linspace(0, h1_star, N);
del_y1 = ymesh1(2) - ymesh1(1);

% Fluid 2
ymesh2 = linspace(h1_star, 1, N);
del_y2 = ymesh2(2) - ymesh2(1);

% Combining the two regions
ymesh = [ymesh1 ymesh2(2:end)];
C3 = alpha / ( del_y2 / del_y1 );

% Modification due to Cheng et al. 2017
C4 = (phi_star(end) - phi_star(end-1)) / del_y1;
C5 = -1 * C1 * C4 * del_y1;

% Note : ymesh2(1) was ignored as its already represented by ymesh1(end)

%-------------------------------------------------------------------------%

% FINITE DIFFERENCE METHOD

Nn = length( ymesh );
A = zeros( Nn, Nn );
b = zeros( Nn, 1 );

% Fluid 1 : Navier Stokes
for i = 2:(N-1)
    A(i, i-1) = 1;
    A(i, i) = -2;
    A(i, i+1) = 1;
    b(i) = (C1 * phi_star(i) * ((K * del_y1)^2)) - (C0 * (del_y1^2));
end

% Fluid 1 : Boundary Cond
A(1, 1) = 1;
b(1) = 0;

% Fluid 2 : Navier Stokes
for i = (N+1):(2*(N-1))
    A(i, i-1) = 1;
    A(i, i) = -2;
    A(i, i+1) = 1;
    b(i) =  -1 * C2 * (del_y2^2);
end

% Fluid 2 : Boundary Cond
A(end, end) = 1;
b(end) = 0;

% Interface : Boundary Cond
A(N, N-1) = 1;
A(N, N) = -1 - C3;
A(N, N+1) = C3;

if BC_interface == "ngoma"
    b(N) = 0;
elseif BC_interface == "cheng"
    b(N) = C5;
else
   error('Please specify type of Boundary Condn : "ngoma" or "cheng"')
end

% Output
u_star = inv(A) * b;

end

%-------------------------------------------------------------------------%