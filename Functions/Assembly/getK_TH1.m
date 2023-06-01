function [K, R, Fint] = getK_TH1(~, ~, ~, Fext, Fextnm1, Klin, ~, d, dnm1, ~, dt, ~, C, alpha)
%GETK_TH1 Conductivity matrix for transient diffusion case
%   [K, R, Fint] = GETK_ELASTIC(Mesh, Quad, Material) returns the stiffness
%   matrix K, the residual vector R, and the internal force vector for the 
%   iterative solver where the problem uses a linear elastic material
%   
%   Template file for other tangent matrix files 
%   --------------------------------------------------------------------
%   Accepted Inputs (in order)
%   --------------------------------------------------------------------
%   getK_TH1(~, ~, ~, Fext, Fextnm1, Klin, ~, d, dnm1, ~, dt, ~, C, alpha)
%
%   Fext:       External flux vector at timestep n
%   Fextnm1:    External flux vector at timestep n-1
%   Klin:       Linear conductivity matrix
%   d:          unconverged degree of freedom vector at current timestep n and iteration
%   dnm1:       converged degree of freedom vector at timestep n-1
%   dt:         timestep size between timesteps n-1 and n
%   alpha:      time intagration parameter


K = 1/dt*C + alpha*Klin ;
Fint = 1/dt*C*(d - dnm1) + (1-alpha)*Klin*dnm1 + alpha*Klin*d;
R = alpha*Fext + (1-alpha)*Fextnm1 - Fint;

end