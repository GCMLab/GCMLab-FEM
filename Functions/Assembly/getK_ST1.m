function [K, R, Fint] = getK_ST1(Mesh, Quad, Material, Fext, ~, ~, d_m, ~, ~, ~, ~)
%GETK_NLELASTIC Stiffness matrix for iterative non linear elastic case
%   [K, R, Fint] = GETK_ST1(Mesh, Quad, Material) returns the 
%   stiffness matrix K, the residual vector R, and the internal force 
%   vector for the iterative solver where the problem uses a non linear 
%   elastic material
%   
%   Template file for other tangent matrix files 
%   --------------------------------------------------------------------
%   Accepted Inputs (in order)
%   --------------------------------------------------------------------
%   getK_NLelastic(Mesh, Quad, Material, Klin, M, d, dnm1, dnm2, stress, strain, dt, dtnm1)
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

d = d_m.d; %Get data from d_m structure

% initialize D matrix file pointer
[~,DMatrix_functn] = fileparts(Material.ConstitutiveLawFile);

% initialize stiffness matrix
vec_size = Mesh.ne*(Mesh.nne * Mesh.nsd)^2; % vector size (solid dofs)
row = zeros(vec_size, 1);                   % vector of row indices
col = zeros(vec_size, 1);                   % vector of column indices
Kvec = zeros(vec_size, 1);                  % vectorized stiffness matrix

count = 1;                                  % DOF counter

% for each element, compute element stiffness matrix and add to global
for e = 1:Mesh.ne

    %% Element variables
        % nodal ids of the element's nodes
        enodes = Mesh.conn(e,:);    
        % global coordinates of the element's nodes
        xI = Mesh.x(enodes,:);  
        % DOFs of element nodes
        dofE = Mesh.DOF(enodes,:);
        dofE = reshape(dofE',Mesh.nDOFe,[]);
        % number of degrees of freedom
        ndofE = numel(dofE);
        % displacements of the element's nodes
        de = d(dofE);
        % element material type
        nMat = Mesh.MatList(e); 

    %% Shape functions and derivatives in parent coordinates
        W = Quad.W;
        nq = Quad.nq;

    %% Assemble stiffness matrix
    
        % length of element. used to check that quadrature points and weights 
        % are correct. A = Sum Wi*Ji
        A = 0;

        % initialize element stiffness matrix
        Ke = zeros(ndofE, ndofE); 

        % loop through all quadrature points
        for q = 1:nq     

            % Shape functions and derivatives in parent coordinates
            N = Quad.Nq{q};
            dNdxi = Quad.dNdxiq{q};

            % quadrature point in physical coordinates
            Xi = xI'*N;

            % Jacobian of the transformation between parent and global 
            % coordinates
            Je = dNdxi'*xI;
           
            % determinant of the Jacobian
            dJe = det(Je);
            if dJe < 0
               error('Element %d has a negative Jacobian.', e)
            end

            % derivative of shape function in physical coordinates 
            % (tensor form)
            dNdxi = dNdxi';
            B = Je\dNdxi;

            % convert B matrix to Voigt form
            Bv = getBv(B', Mesh.nsd);
            
            % for 2D, volume integral includes the thickness
            switch Mesh.nsd 
                case 1
                    L = Material.t(Xi);
                case 2
                    L = Material.t(Xi);
                case 3
                    L = 1;                
            end

            % compute strains in quadrature point
            strain_e = Bv'*de;

            % Compute constitutive matrix
            D = feval(DMatrix_functn, nMat, Material, Mesh, strain_e);

            % Calculate local stiffness matrix
            Ke = Ke + W(q)*Bv*D*Bv'*L*dJe;

            % quadrature debug tool
            A = A + W(q)*dJe; 
        end
   
    %% Forming the vectorized stiffness matrix
        count = count + ndofE^2;
        Ke = reshape(Ke, [ndofE^2, 1]);
        rowmatrix = dofE*ones(1,ndofE);        
        rowe = reshape(rowmatrix,[ndofE^2, 1]);
        cole = reshape(rowmatrix',[ndofE^2, 1]);
        
        Kvec(count-ndofE^2:count-1) = Ke;
        row(count-ndofE^2:count-1) = rowe;
        col(count-ndofE^2:count-1) = cole;
end

% sparse stiffness matrix
K = sparse(row, col, Kvec, Mesh.nDOF, Mesh.nDOF);

% internal forces
Fint = K*d;

% residual
R = Fext - Fint;

end