function [K, R, Fint] = getK_THLE2(~, ~, ~, ~, Fext, Fextnm1, Klin, ~, d_m, dt, ~,C,alpha)
%GETK_THLE2 Stiffness matrix for iterative transient thermoelastic case
%   [K, R, Fint] = GETK_THLE2(Mesh, Quad, Material) returns the stiffness
%   matrix K, the residual vector R, and the internal force vector for the 
%   iterative solver where the problem uses a linear transient 
%   thermoelastic material
%   
%   Template file for other tangent matrix files 
%   --------------------------------------------------------------------
%   Accepted Inputs (in order)
%   --------------------------------------------------------------------
%   Accepted Inputs (in order)
%   --------------------------------------------------------------------
%   getK_LET1(~, ~, ~, ~, Fext, Fextnm1, Klin, ~, d_m, dt, ~,C,alpha)
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
%   C:          Material damping matrix
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