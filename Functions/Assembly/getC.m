function C = getC(Mesh, Quad, Material)
%GETC Damping/Capacity matrix 
%   C = GETC(Mesh, Quad, Material) is the global damping matrix for 
%   parameters defined in the structure arrays Mesh, Quad, and Material. 
%   The sparse matrix has size Mesh.nDOF x Mesh.nDOF.
%
%   This function also assembles the capacity matrix for diffusion
%   problems. In this case, the damping coefficient is replaced with the
%   heat capacity: (specific heat * density) = [J/K kg]*[kg/m^3] = [J/K m^3]
%   
%   --------------------------------------------------------------------
%   Input
%   --------------------------------------------------------------------
%   Mesh:       Structure array with the following fields,
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
%   Quad:       Structure array with the following fields,
%               .W:      Vector of quadrature weights (size nq x 1)      
%               .nq:     Number of quadrature points 
%               .Nq:     Cell array (size nq x 1) with shape functions  
%                        evaluated at each quadrature point
%               .dNdxiq: Cell array (size nq x 1) with derivative of shape 
%                        functions w.r.t. parent coordinates evaluated at 
%                        each quadrature point
% 
%   Material:   Structure array with the following fields,
%               .t:         Material thickness
%               .C:         Damping Coefficient / heat capacity

% Acknowledgements: Jonathan Zingaro

% initialize damping matrix
vec_size = Mesh.ne*(Mesh.nne * Mesh.nDOFn)^2; % vector size (solid dofs)
row = zeros(vec_size, 1);                   % vector of row indices
col = zeros(vec_size, 1);                   % vector of column indices
Cvec = zeros(vec_size, 1);                  % vectorized damping matrix

count = 1;                                  % DOF counter

% for each element, compute element damping matrix and add to global
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
        
    %% Retrieve Material Properties of Element
        nMat = Mesh.MatList(e); % element material type        
        ce = Material.Prop(nMat).C; % viscosity of element material properties
    %% Shape Functions and Derivatives in Parent Coordinates
        W = Quad.W;
        nq = Quad.nq;

    %% Assemble Damping matrix
    
        % length of element. used to check that quadrature points and weights 
        % are correct. A = Sum Wi*Ji
        A = 0;

        % initialize element damping matrix
        Ce = zeros(ndofE, ndofE); 

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

            % Convert N to Matrix to Voigt Form
            Nv = getNv(N, Mesh.nDOFn);

            % for 2D, volume integral includes the thickness
            switch Mesh.nsd 
                case 1
                    L = Material.t(Xi);
                case 2
                    L = Material.t(Xi);
                case 3
                    L = 1;                
            end

            % Calculate local damping matrix
            Ce = Ce + W(q)*Nv*ce*Nv'*L*dJe;

            % quadrature debug tool
            A = A + W(q)*dJe; 
        end
   
    %% Forming the vectorized damping matrix
        count = count + ndofE^2;
        Ce = reshape(Ce, [ndofE^2, 1]);
        rowmatrix = dofE*ones(1,ndofE);        
        rowe = reshape(rowmatrix,[ndofE^2, 1]);
        cole = reshape(rowmatrix',[ndofE^2, 1]);
        
        Cvec(count-ndofE^2:count-1) = Ce;
        row(count-ndofE^2:count-1) = rowe;
        col(count-ndofE^2:count-1) = cole;
end

% sparse stiffness matrix
C = sparse(row, col, Cvec, Mesh.nDOF, Mesh.nDOF);

end