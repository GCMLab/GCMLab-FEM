function [K, R, Fint] = getK_TH1(~, ~, ~, ~, Fext, Fextnm1, Klin, ~, d_m, dt, ~, C, alpha)
%GETK_TH1 Conductivity matrix for diffusion case
%   K = GETK_TH1(Mesh, Quad, Material) returns the stiffness
%   matrix K for the iterative solver where the problem uses a linear elastic material
%
%   [K, R] = GETK_TH1(Mesh, Quad, Material) also returns the residual vector R for the 
%   iterative solver where the problem uses a linear elastic material
%
%   [K, R, Fint] = GETK_TH1(Mesh, Quad, Material) also returns the internal force vector for the 
%   iterative solver where the problem uses a linear elastic material
%
%   This file generates the stiffness matrices for both TH1 (Steady-state)
%   and TH2 (transient) diffusion cases. In the steady-state case, C is 
%   automatically set to zero and alpha is automatically set to 1.
%
%   --------------------------------------------------------------------
%   Accepted Inputs (in order)
%   --------------------------------------------------------------------
%   getK_TH1(~, ~, ~, Fext, Fextnm1, Klin, ~, d_m, dt, ~, C, alpha)
%
%   Fext:       External flux vector at timestep n
%   Fextnm1:    External flux vector at timestep n-1
%   Klin:       Linear conductivity matrix
%   d:          unconverged degree of freedom vector at current timestep n and iteration
%   dnm1:       converged degree of freedom vector at timestep n-1
%   dt:         timestep size between timesteps n-1 and n
%   alpha:      time intagration parameter

dnm1 = d_m.dnm1;
d = d_m.d;

% stiffness matrix in transient case
K = alpha*Klin + C./dt;

% internal forces
Fint = C*(d-dnm1)./dt + (1-alpha)*Klin*dnm1 + alpha*Klin*d;

% residual
R = alpha*Fext + (1-alpha)*Fextnm1 - Fint;

end
