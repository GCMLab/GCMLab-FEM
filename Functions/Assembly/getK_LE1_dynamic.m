function [K, R, Fint] = getK_LE1_dynamic(~, ~, ~, Fext, ~, Klin, M, d_m, dt, ~ ,C,alpha)
%GETK_HHT Relevant matrices for HHT dynamic solver
%   [K, R, Fint] = GETK_HHT(Mesh, Quad, Material) returns the coefficient matrix
%   , the residual vector R, and the internal force vector for the 
%   iterative solver where the problem uses a linear elastic material
%   
%   Template file for other tangent matrix files 
%   --------------------------------------------------------------------
%   Accepted Inputs (in order)
%   --------------------------------------------------------------------
%   getK_elastic(Mesh, Quad, Material, Klin, M, d, dnm1, dnm2, stress, strain, dt, dtnm1)
%   Mesh:       Structure array with the following fields, may be updated
%               with new fields
%               .ne:    Total number of elements in the mesh
%               .nne:   Vector of number of nodes per element (size nsd x 1)
%               .nsd:   Number of spatial dimensions
%               .conn:  Array of element connectivity (size ne x nne)
%               .x:     Array of nodal spatial locations for
%                       undeformed mesh (size nn x nsd)
%               .DOF:   Array of DOF indices (size nn x nsd)
%               .nDOFe: Number of DOFs per element
%               .nDOF:  Total number of DOFs
%  
%   Quad:       Structure array with the following fields, may be updated
%               with new fields
%               .W:      Vector of quadrature weights (size nq x 1)      
%               .nq:     Number of quadrature points 
%               .Nq:     Cell array (size nq x 1) with shape functions  
%                        evaluated at each quadrature point
%               .dNdxiq: Cell array (size nq x 1) with derivative of shape 
%                        functions w.r.t. parent coordinates evaluated at 
%                        each quadrature point
% 
%   Material:   Structure array with the following fields, may be updated
%               with new fields
%               .t:         Material thickness
%
%   Fext:       External force vector at timestep n
%   Fextnm1:    External force vector at timestep n-1
%   Klin:       Linear elastic stiffness matrix
%   M:          Mass matrix
%   d_m:        Structure array with the following fields
%               d:          unconverged degree of freedom vector at current timestep n and iteration
%               dnm1:       converged degree of freedom vector at timestep n-1
%               dnm2:       converged degree of freedom vector at timestep n-2
%               dnm3:       converged degree of freedom vector at timestep n-3
%   dt:         timestep size between timesteps n-1 and n
%   dtnm1:      timestep size between timesteps n-2 and n-1
%   alpha:       intagration parameter
%   C:          Linear damping coefficient matrix

% Compute constants
gam = 1/2-alpha;
bet = (1-alpha)^2/4;

d = d_m.d;
dnm1 = d_m.dnm1;
dnm2 = d_m.dnm2;
dnm3 = d_m.dnm3;
dnm4 = d_m.dnm4;

% 2nd order accurate backwards difference approximation
vnm1 = 1/2/dt* (3*dnm1 - 4*dnm2 + dnm3);
anm1 = 1/dt^2 * (2*dnm1  - 5*dnm2 + 4*dnm3 - dnm4);

d_temp = dnm1+ dt*vnm1 + dt^2/2*(1-2*bet)*anm1;
v_temp = vnm1 + dt*(1-gam)*anm1;

% Internal forces
a = (d - d_temp)./dt^2/bet;
Fint = M*a +(1+alpha)*C*(v_temp-gam*d_temp/(dt*bet)+gam*d/(dt*bet))...
    +(1+alpha)*Klin*d - alpha*(C*vnm1 + Klin*dnm1);

% Jacobian
K = M/dt^2/bet + (1+alpha)*gam/dt/bet*C + (1+alpha)*Klin;

% Residual
R = Fext + alpha*(C*vnm1 + K*dnm1) + M*d_temp/dt^2/bet -...
    (1+alpha)*C*(v_temp-gam*d_temp/dt/bet);

end
