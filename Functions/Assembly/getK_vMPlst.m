function [K_hat, R_hat, Fint, stress, strain, strain_e, acc_plastic] = getK_vMPlst(Mesh, Quad, Material, ~, Fext, ~, ~, ~, d, ~, ~, Dd, ~, ~, strain_e, acc_plastic)

%GETK_TR2 Stiffness matrix for iterative non linear elastic transient case
%   K_hat = getK_TR2(Mesh, Quad, Material) returns the stiffness matrix K 
%   for the iterative solver where the problem uses a nonlinear transient elastic material
%   
%   [K_hat, R_hat] = getK_TR2(Mesh, Quad, Material) also returns the residual vector R for the 
%   iterative solver where the problem uses a nonlinear elastic transient material
%
%   [K_hat, R_hat, Fint] = getK_TR2(Mesh, Quad, Material) also returns 
%   the internal force vector for the iterative solver where the problem uses 
%   a nonlinear elastic transient material
%
%   Template file for other tangent matrix files 
%   --------------------------------------------------------------------
%   Accepted Inputs (in order)
%   --------------------------------------------------------------------
%   getK_TR2(Mesh, Quad, Material, Klin, M, d, dnm1, dnm2, stress, strain, dt, dtnm1)
%
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
%   d:          unconverged degree of freedom vector at current timestep n and iteration
%   dnm1:       converged degree of freedom vector at timestep n-1
%   dnm2:       converged degree of freedom vector at timestep n-2
%   dt:         timestep size between timesteps n-1 and n
%   dtnm1:      timestep size between timesteps n-2 and n-1
%   alpha:       intagration parameter

% initialize D matrix file pointer
[~,DMatrix_functn] = fileparts(Material.ConstitutiveLawFile);

% initialize stiffness matrix
vec_size = Mesh.ne*(Mesh.nne * Mesh.nsd)^2; % vector size (solid dofs)
row = zeros(vec_size, 1);                   % vector of row indices
col = zeros(vec_size, 1);                   % vector of column indices
Kvec = zeros(vec_size, 1);                  % vectorized stiffness matrix

count = 1;                                  % DOF counter

% initialize internal force vector
Fint = zeros(Mesh.nDOF, 1);  

stress = zeros(Mesh.ne, Quad.nq, 3);
strain = zeros(Mesh.ne, Quad.nq, 3);

% initialize d vector
dnm1 = d.dnm1;
d = d.d;


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
        Dde = Dd(dofE);
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
        
        % initialize element internal force vector
        fe = zeros(ndofE, 1); 

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
            Dstrain_e = Bv'*Dde;
            

            % Compute constitutive matrix
            [D] = feval(DMatrix_functn, nMat, Material, Mesh);
            
            ee(1,1) = strain_e(e,q,1);
            ee(2,1) = strain_e(e,q,2);
            ee(3,1) = strain_e(e,q,3);
            
            strain_e(e,q,:) = Dstrain_e + ee;
            strain(e,q,:) = Bv'*de;
            stress(e,q,:) = D*ee;
            

            % Calculate local stiffness matrix
            Ke = Ke + W(q)*Bv*D*Bv'*L*dJe;
            % Ke = Ke + W(q)*Bv*D*Bv'*L*dJe;
            
            % Calculate local internal force vector
            fe = fe + W(q)*Bv*D*(Dstrain_e + ee)*L*dJe;

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
        
        %% Assemble element internal forces
        Fint(dofE) = Fint(dofE) + fe;
end

% sparse stiffness matrix
K = sparse(row, col, Kvec, Mesh.nDOF, Mesh.nDOF);

% stiffness matrix in transient case
K_hat = K;

% internal forces
% Fint = C*(d-dnm1)./dt + K*(1-alpha)*dnm1 + alpha*K*d;

% residual
R_hat = Fext - Fint;


end