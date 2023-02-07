function [disp_er] = Test1D_check(d, Material, BC, Mesh)
%TEST1D_CHECK Calculates the error between FEA and analytical solutions
%   [disp_er] = Test1D_check(d, Material, BC, Mesh)
%   calculates the error between FEA and analytical displacements along the
%   free edge of a one dimensional problem
%
%   ----------------------------------------------------------
%   Input
%   ----------------------------------------------------------
%   d:                  Displacement vectors
%   BC:                 Boundary Conditions
%   Material:           Material Properties
%   Mesh:               Mesh data structure
%
%   ----------------------------------------------------------
%   Output
%   ----------------------------------------------------------
%   disp_er:             Error related to Displacement

% The analytical displacement along the free edge of a 1D bar is given by
% d = F*L/(E*A) + b*L^2/(2*E*A)
% F, b and A variables are the applied force along the free edge, distributed body force and
% sectional area of the bar, respectively.

global E traction b

% Calculate analytical displacement
% prescribed traction [N]
Force = traction;
% distributed body force [N/m]
Body_force = b;
L = max(Mesh.x(:,1));
rightnode = BC.traction_force_node;
d_exact = Force*L/E/Material.t(Mesh.x(rightnode,1)) + Body_force*L^2/(2*E*Material.t(Mesh.x(rightnode,1)));

% Calculate the error
disp_er = abs(d(rightnode) - d_exact)/abs(d_exact);

end