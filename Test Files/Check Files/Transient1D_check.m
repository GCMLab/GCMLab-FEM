function [L2d] = Transient1D_check(d, Material, ~, Control,BC)
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
%% NOTE: 
% 1. the external force is applied to all nodes in the x-direction along the free-end of the beam. 
%    This is accounted by F = 21*Fn. 
% 2. Index 863 corresponds to the displacement in x along the free-end of the beam at the midpoint in the y-direction. 

% Import Material Parameters
C = Material.Prop(1).C;
K = Material.Prop(1).E0;

% Time-step Parameters
t = Control.StartTime+Control.TimeStep:Control.TimeStep:Control.EndTime;

% Force at Free End
F = BC.Fn;

% Compute Exact Solution
dsoln = @(t) ((F./(C.^2+K.^2)).*(-C.*cos(t)+K.*sin(t)+C.*cos(t).*exp(-K.*t./C)))*10;
dexact = dsoln(t);

% Index Displacements at Mid-Height at the Free-end of Bar
dnum = d(863,2:end);

% Calculate the error
L2d = norm(dnum-dexact);

end
