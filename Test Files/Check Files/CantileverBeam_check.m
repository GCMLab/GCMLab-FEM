function [disp_er, time_er] = CantileverBeam_check(d, Material, BC, Mesh)
%CANTILEVERBEAM_CHECK Calculates the error between FEA and analytical solutions
%   [disp_er] = CantileverBeam_check(d, Material, BC, Mesh)
%   calculates the error between FEA and analytical displacements along the
%   free edge of a cantilever beam
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

% The analytical displacement along the free edge of a cantilever beam is given by
% d_y = F*Lx^3/(3*E*I)
% F and I variables are the applied force along the free edge and
% moment of inertia of the beam section, respectively.
% F = t*Ly*Material.t
% I = (1/12)*Material.t*Ly^3

% Calculate analytical displacement
sigma = BC.traction;
Lx = max(Mesh.x(:,1));
Ly = max(Mesh.x(:,2));
toprightnode = BC.traction_force_node(find(Mesh.x(BC.traction_force_node,2) == max(Mesh.x(:,2))));
I = (1/12)*Material.t([Mesh.x(toprightnode,1),Mesh.x(toprightnode,2)])*Ly^3;
d_exact = sigma*Material.t([Mesh.x(toprightnode,1),Mesh.x(toprightnode,2)])*Ly*(Lx)^3/3/Material.Prop(1).E0/I;

% Calculate the error
disp_er = abs(d(toprightnode*2,3) - d_exact)/abs(d_exact);

% Time error
tt = linspace(0,2*pi,9);
d_time = d_exact*sin(tt);

time_er = sqrt(sum((d_time - d(toprightnode*2,:)).^2))/sqrt(sum(d_time.^2));


end