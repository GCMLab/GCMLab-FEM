function [TH_er] = Diffusion_check(d, Mesh)
%DIFFUSION_CHECK Calculates the error between FEA and analytical solutions
%   [TH_er] = DIFFUSION_CHECK(d, Mesh)
%   calculates the error between FEA and analytical displacements at the
%   center point of the domain
%
%   ----------------------------------------------------------
%   Input
%   ----------------------------------------------------------
%   d:                  Displacement vectors
%   Mesh:               Mesh data structure
%
%   ----------------------------------------------------------
%   Output
%   ----------------------------------------------------------
%   TH_er:             Error related to Temperature

% The analytical solution to the reference problem is 
% T(x,y) = sin(pi*x) * sin(pi*y)

% Calculate analytical temperature at reference node ni
ni = 14;
xp = Mesh.x(ni,1);
yp = Mesh.x(ni,2);
T_exact = sin(pi*xp)*sin(pi*yp);
T_FEA = d(ni);

% Calculate the error
TH_er = abs(T_FEA - T_exact)/abs(T_exact);


end