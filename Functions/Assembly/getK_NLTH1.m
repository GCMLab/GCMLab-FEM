function [H, R, Fint] = getK_NLTH1(Mesh, Quad, Material, ~, Fext, Fextnm1, ~, ~, d_m, dt, ~, C, beta)
%GETK_NLTH1 computes the Hessian matrix (H) and the residual vector for a  
%   nonlinear heat diffusion problem, where the nonlinearity arises as
%   k(T) = k_0 + αT^n. 
%   [K, R, Fint] = getK_NLTH1(Mesh, Quad, Material, ~, Fext, Fextnm1, ~, ~,
%   d_m, dt, ~, C, beta) returns the Hessian matrix (H), the residual 
%   vector (R), and the internal force vector (Fint) used in the 
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
%   H:          Hessian matrix for nonlinear heat diffusion
%   K:          Stiffness matrix for nonlinear heat diffusion at timestep n
%   Knm1:       Stiffness matrix for nonlinear heat diffusion at timestep n-1
%   M:          Mass matrix
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
% Acknowledgements: Jonathan Zingaro

d = d_m.d; % Get data from d_m structure for timestep n
dnm1 = d_m.dnm1; % Get temperature from d_m structure for timestep n-1

% initialize stiffness matrix
vec_size = Mesh.ne*(Mesh.nne)^2; % vector size (solid dofs)
row = zeros(vec_size, 1);                   % vector of row indices
col = zeros(vec_size, 1);                   % vector of column indices
Hvec = zeros(vec_size, 1);                  % vectorized stiffness matrix
count = 1;                                  % DOF counter

% initialize stiffness matrices for K and Knm1
Kvec = zeros(vec_size, 1);                  % vectorized stiffness matrix
Knm1vec = zeros(vec_size, 1);               % vectorized stiffness matrix

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
        denm1 = dnm1(dofE);
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
        He = zeros(ndofE, ndofE); 
        Ke = zeros(ndofE, ndofE); 
        Knm1e = zeros(ndofE, ndofE); 
      
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

            % Tangent Matrix 
            % Compute Temperature at Gaussian Point
            dG = N'*de;

            % Constutitive Law for Nonlinear Diffusion -> TH3
            [D,alpha] = getD_NLTH1(nMat, Material, Mesh, dG);

            % Calculate Linear Portion of Hessian Matrix
            He = He + W(q)*B'*D*B*dJe*L;

            % Contruct Remaining Nonlinear Portion of Hessian Matrix
            for i = 1:Mesh.nne
                for j = 1:Mesh.nne
                    He(i,j) = He(i,j) + W(q)*B(:,i)'*n*alpha*(dG.^(n-1))*B(:,j)*dG*dJe*L;
                end
            end

            % Compute Elemental Stiffness Matrix at timestep n
            Ke = Ke + W(q)*B'*D*B*dJe*L;

            % Compute Elemental Stiffness Matrix at timestep n-1
            dGnm1 = N'*denm1; % Temperature at Gaussian Point at timestep n-1
            Dnm1 = getD_NLTH1(nMat, Material, Mesh, dGnm1);
            Knm1e = Knm1e + W(q)*B'*Dnm1*B*dJe*L;
          
            % quadrature debug tool
            A = A + W(q)*dJe; 
        end
   
    %% Forming the vectorized stiffness matrix
        count = count + ndofE^2;
        He = reshape(He, [ndofE^2, 1]);
        Ke = reshape(Ke, [ndofE^2, 1]);
        Knm1e = reshape(Knm1e, [ndofE^2, 1]);
        rowmatrix = dofE*ones(1,ndofE);        
        rowe = reshape(rowmatrix,[ndofE^2, 1]);
        cole = reshape(rowmatrix',[ndofE^2, 1]);
        
        Hvec(count-ndofE^2:count-1) = He;
        Kvec(count-ndofE^2:count-1) = Ke;
        Knm1vec(count-ndofE^2:count-1) = Knm1e;
        row(count-ndofE^2:count-1) = rowe;
        col(count-ndofE^2:count-1) = cole;
       
end

% Sparse Stiffness Matrix
H = sparse(row, col, Hvec, Mesh.nDOF, Mesh.nDOF);
K = sparse(row, col, Kvec, Mesh.nDOF, Mesh.nDOF);
Knm1 = sparse(row, col, Knm1vec, Mesh.nDOF, Mesh.nDOF);

% Final  Hessian Matrix  
H = beta*H + C./dt;

% Residual
Fint = C*(d-dnm1)./dt + (1-beta)*Knm1*dnm1 + beta*K*d;
R = beta*Fext + (1-beta)*Fextnm1 - Fint;

end
