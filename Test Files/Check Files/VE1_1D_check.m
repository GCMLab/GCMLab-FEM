function [L2d] = VE1_1D_check(d, Material, ~, Control,BC)
%TEST1D_CHECK Calculates the error between FEA and Kelvin-Voigt anyaltical
%   solutions for a 1D bar. 
%   [disp_er] = Test1D_check(d, Material, BC, Mesh)
%   Calculates the error between FEA and analytical displacements along the
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
% 3. Anyalitical Solution for Creep Loading ε(t) = σ0/E*(1-exp(-Et/C))
% 4. Analytical Solution for Recovery ε(t) =
%    σ0/E*exp(-Et/C)*(exp(Et1/C)-1), where t1 is the time at which the load is
%    released.

% Import Material Parameters
C = Material.Prop(1).C;
K = Material.Prop(1).E0;
t1 = 0.5; % Time at Load Removed

% Time-step Parameters
tcreep = 0:Control.TimeStep:t1;
trel = t1+Control.TimeStep:Control.TimeStep:Control.EndTime;
t = Control.StartTime:Control.TimeStep:Control.EndTime;
% Force at Free End
sigma_0 = BC.Fn_total/0.1; % Applied Traction at Free-End [Pa]

% Compute Exact Solution
dcreep = @(t) sigma_0/K*(1-exp(-K*t/C)); % Displacement Under Creep for 1D KV Model
drel = @(t) sigma_0/K*exp(-K*t/C)*(exp(K*tcreep(end)/C)-1);
dexact = [dcreep(tcreep),drel(trel)];


% Index Displacements at Mid-Height at the Free-end of Bar
dnum = d(863,:);

% Calculate the error
L2d = norm(dnum-dexact)/norm(dexact);

end
