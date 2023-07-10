function [K, R, Fint] = getK_LE1(~, ~, ~, ~, Fext, ~, Klin, ~, d_m, ~, ~,~,~)
%GETK_LE1 Stiffness matrix for iterative elastic case
%   [K, R, Fint] = GETK_LE1(Mesh, Quad, Material) returns the stiffness
%   matrix K, the residual vector R, and the internal force vector for the 
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
%   Klin:       Linear elastic stiffness matrix
%   M:          Mass matrix
%   d_m:        Structure array with the following fields
%               d:          unconverged degree of freedom vector at current timestep n and iteration
%               dnm1:       converged degree of freedom vector at timestep n-1
%               dnm2:       converged degree of freedom vector at timestep n-2
%               dnm3:       converged degree of freedom vector at timestep n-3
%   dt:         timestep size between timesteps n-1 and n
%   dtnm1:      timestep size between timesteps n-2 and n-1

d = d_m.d;

K = Klin;
Fint = K*d;
R = Fext - Fint;


end