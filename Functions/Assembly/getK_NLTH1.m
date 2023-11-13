function [K_hat, R_hat, Fint] = getK_NLTH1(Mesh, Quad, Material, Fintnm1, Fext, Fextnm1, ~, ~, d_m, dt, ~, C, beta)
%GETK_NLTH1 computes the Hessian matrix (H) and the residual vector for a  
%   nonlinear heat diffusion problem, where the nonlinearity arises as
%   k(T) = k_0 + αT^n. 
%
%   [K_hat] = getK_NLTH1(Mesh, Quad, Material, ~, Fext, Fextnm1, ~, ~,
%   d_m, dt, ~, C, beta) returns the stiffness matrix (Hessian matrix) used in the 
%   computation of the Newton-Raphson solver. 

%   [K_hat, R_hat] = getK_NLTH1(Mesh, Quad, Material, ~, Fext, Fextnm1, ~, ~,
%   d_m, dt, ~, C, beta) returns the stiffness matrix (Hessian matrix),and the residual 
%   vector used in the computation of the Newton-Raphson solver. 
%
%   [K_hat, R_hat, Fint] = getK_NLTH1(Mesh, Quad, Material, ~, Fext, Fextnm1, ~, ~,
%   d_m, dt, ~, C, beta) returns the stiffness matrix (Hessian matrix), the residual 
%   vector, and the internal force vector used in the 
%   computation of the Newton-Raphson solver. 
%   
%   --------------------------------------------------------------------
%   Accepted Inputs (in order)
%   --------------------------------------------------------------------
%   getK_NLTH1(Mesh, Quad, Material, ~, Fext, Fextnm1, ~, ~, d_m, dt, ~, C, beta)
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
%   K:          Hessian matrix for nonlinear heat diffusion
%   K_hat:      Equivalent stiffness matrix for nonlinear heat diffusion
%   R_hat:      Equivalent residual force for nonlinear heat diffusion
%   Fint:       Internal force for nonlinear heat diffusion at timestep n
%   Fintnm1:    Internal force for nonlinear heat diffusion at timestep n-1
%   C:          Heat capacity matrix
%   d_m:        Structure array with the following fields
%               d:          unconverged degree of freedom vector at current
%                           timestep n and iteration
%               de:         unconverged degree of freedom vector at current timestep
%               denm1:      unconverged degree of freedom vector at timestep n-1 for element e 
%               denm2:      unconverged degree of freedom vector at timestep n-2 for element e 
%               denm3:      unconverged degree of freedom vector at timestep n-3 for element e
%               dG:         temperature at Gaussian Point at current timestep
%               dGnm1:      temperature at Gaussian Point at timestep n-1 for element e
%               dGnm2:      temperature at Gaussian Point at timestep n-1 for element e
%               dGnm3:      temperature at Gaussian Point at timestep n-1 for element e
%               dnm1:       converged degree of freedom vector at timestep n-1
%               dnm2:       converged degree of freedom vector at timestep n-2
%               dnm3:       converged degree of freedom vector at timestep n-3
%   dt:         timestep size between timesteps n-1 and n
%   dtnm1:      timestep size between timesteps n-2 and n-1
%
% Acknowledgements: Jonathan Zingaro, Saeed Hatefi Ardakani



% initialize stiffness matrix
vec_size = Mesh.ne*(Mesh.nne)^2; % vector size (solid dofs)
row = zeros(vec_size, 1);                   % vector of row indices
col = zeros(vec_size, 1);                   % vector of column indices
Kvec = zeros(vec_size, 1);                  % vectorized stiffness matrix

count = 1;                                  % DOF counter

% initialize internal force vector
Fint = zeros(Mesh.nDOF, 1);  

% initialize d vector
d = d_m.d; % Get data from d_m structure for timestep n
dnm1 = d_m.dnm1; % Get temperature from d_m structure for timestep n-1

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
        % temperature of the element's nodes
        de = d(dofE);
        % element material type
        nMat = Mesh.MatList(e); 
        % Constutitive Law Parameters
        n = Material.Prop(nMat).n;

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

            % for 2D, volume integral includes the thickness
            switch Mesh.nsd 
                case 1
                    L = Material.t(Xi);
                case 2
                    L = Material.t(Xi);         
            end

            % Compute Temperature at Gaussian Point
            dG = N'*de;

            % Constutitive Law for Nonlinear Diffusion
            [D, alpha] = getD_NLTH1(nMat, Material, Mesh, dG);

            % Calculate Hessian Matrix
            Ke = Ke + W(q)*B'*D*B*dJe*L + W(q)*B'*n*alpha*(dG.^(n-1))*B*de*N'*dJe*L;

            % Calculate local internal force vector
            fe = fe + W(q)*B'*D*B*de*L*dJe;

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

% Sparse Stiffness Matrix
K = sparse(row, col, Kvec, Mesh.nDOF, Mesh.nDOF);

% Stiffness matrix in transient case
K_hat = beta*K + C./dt;

% Residual
R_hat = beta*Fext + (1-beta)*Fextnm1 - (C*(d-dnm1)./dt + beta*Fint + (1-beta)*Fintnm1);

end
