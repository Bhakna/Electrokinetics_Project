function [ sol ] = analyticalSol()

data = readstruct("nondimData.xml");

%-------------------------------------------------------------------------%

% GENERATING RESULTS

% 1-D Mesh/Grid
sol.y_star = linspace(0, data.h1_star, data.N_y);
sol.z_star = linspace(-data.l_star, data.l_star, data.N_z);
[sol.Y_mesh, sol.Z_mesh] = meshgrid(sol.y_star, sol.z_star);

% Electric Potential
[ sol.phi_star ] = ...
    electricPotential(data.De, data.h1_star, data.l_star, ...
    data.Z_b_star, data.Z_s_star, data.Z_int_star, ...
    sol.Y_mesh, sol.Z_mesh);

% Velocity Field
if data.anaSolFunc == "mathematica"

[ sol.u_star ] = ...
    velocityMath(data.gamma_r, ...
    data.h1_star, data.l_star, ...
    data.Z_b_star, data.Z_s_star, ...
    sol.phi_star, sol.Y_mesh, sol.Z_mesh);

elseif (data.anaSolFunc == "paper") || (data.anaSolFunc == "chakraborty")

[ sol.u_star ] = ...
    velocityPaper(data.gamma_r, ...
    data.h1_star, data.l_star, ...
    data.Z_b_star, data.Z_s_star, ...
    sol.phi_star, sol.Y_mesh, sol.Z_mesh);

else
    error("Enter correct analytical sol type : 'mathematica' or 'paper'")
end

end

%-------------------------------------------------------------------------%