function [K, R, Fint] = getK_THLE1(~, ~, ~, ~, Fext, ~, Klin, ~, d_m, ~, ~,~,~)
%GETK_THLE1 Stiffness matrix for iterative thermoelastic case
%   [K, R, Fint] = GETK_THLE1(Mesh, Quad, Material) returns the stiffness
%   matrix K, the residual vector R, and the internal force vector for the 
%   iterative solver where the problem uses a linear thermoelastic material
%   
%   Template file for other tangent matrix files 
%   --------------------------------------------------------------------
%   Accepted Inputs (in order)
%   --------------------------------------------------------------------
%   getK_LET1(~, ~, ~, ~, Fext, Fextnm1, Klin, ~, d_m, dt, ~,~,~)
%   Mesh:       ~
%   Quad:       ~
%   Material:   ~
%   Fintnm1:    ~
%   Fext:       External force vector at timestep n
%   Fextnm1:    External force vector at timestep n-1
%   Klin:       Linear elastic stiffness matrix
%   M:          ~
%   d_m:        Structure array with the following fields
%               d:          unconverged degree of freedom vector at current timestep n and iteration
%               dnm1:       converged degree of freedom vector at timestep n-1
%               dnm2:       converged degree of freedom vector at timestep n-2
%               dnm3:       converged degree of freedom vector at timestep n-3
%   dt:         timestep size between timesteps n-1 and n
%   dtnm1:      ~
%   C:          ~
%   alpha:      ~

d = d_m.d;

K = Klin;
Fint = K*d;
R = Fext - Fint;


end