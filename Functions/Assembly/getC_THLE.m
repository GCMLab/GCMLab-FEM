function C = getC_THLE(Mesh, Quad, Material)
%GETC_THLE Damping matrix for thermoelasticity problem
%   C = GETC_THLE(Mesh, Quad, Material) is the global thermoelasticity 
%   damping matrix for parameters defined in the structure arrays Mesh, 
%   Quad, and Material. The sparse matrix has size Mesh.nDOF x Mesh.nDOF.
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

% Acknowledgements: Jonathan Zingaro, Bruna Campos

% initialize damping matrix
Cttvec_size = Mesh.ne*(Mesh.nne)^2;                     % vector size
Cutvec_size = Mesh.ne*Mesh.nDOFe_mech*Mesh.nDOFe_sca;   % vector size

rowt = zeros(Cttvec_size, 1);                   % vector of row indices
colt = zeros(Cttvec_size, 1);                   % vector of column indices

rowut = zeros(Cutvec_size, 1);                   % vector of row indices
colut = zeros(Cutvec_size, 1);                   % vector of column indices

Cttvec = zeros(Cttvec_size, 1);                  % vectorized damping matrix
Cutvec = zeros(Cutvec_size, 1);                  % vectorized damping matrix

count_t = 1;                                  % DOF counter
count_ut = 1;                                  % DOF counter

% for each element, compute element damping matrix and add to global
for e = 1:Mesh.ne

    %% Element variables
        % nodal ids of the element's nodes
        enodes = Mesh.conn(e,:);    
        % global coordinates of the element's nodes
        xI = Mesh.x(enodes,:);  
        % number of degrees of freedom mechanical field
        ndofE_mech = Mesh.nne*Mesh.nsd;
        % number of degrees of freedom scalar field
        ndofE_sca = Mesh.nne;
        % displacement DOFs (mechanical problem)
        dofE_mech = Mesh.DOF_mech(enodes,:);
        dofE_mech = reshape(dofE_mech', Mesh.nDOFe_mech,[]);
        % temperature DOFs (scalar problem)
        dofE_sca = Mesh.DOF_sca(enodes,:);
        dofE_sca = reshape(dofE_sca', Mesh.nDOFe_sca,[]);
        
    %% Retrieve Material Properties of Element
        nMat = Mesh.MatList(e); % element material type        
        beta = Material.Prop(nMat).beta; % coupled matrix coefficient
        c = Material.Prop(nMat).C; % heat capacity
        T0 = Material.Prop(nMat).T0; % reference temperature
        
    %% Shape Functions and Derivatives in Parent Coordinates
        W = Quad.W;
        nq = Quad.nq;

    %% Assemble Damping matrix
    
        % length of element. used to check that quadrature points and weights 
        % are correct. A = Sum Wi*Ji
        A = 0;

        % initialize element damping matrix
        Ctte = zeros(ndofE_sca, ndofE_sca);
        Cute = zeros(ndofE_mech, ndofE_sca);

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
                    m = 1; % mapping vector
                case 2
                    L = Material.t(Xi);
                    m = [1;1;0]; % mapping vector
                case 3
                    L = 1;                
            end

            % Calculate local damping matrix
            Ctte = Ctte + W(q)*N*c*N'*L*dJe; 
            Cute = Cute + W(q)*beta*Bv*m*N'*dJe;

            % quadrature debug tool
            A = A + W(q)*dJe; 
        end
   
    %% Forming the vectorized damping matrix
    % Heat capacity matrix
        count_t = count_t + ndofE_sca^2;
        Ctte = reshape(Ctte, [ndofE_sca^2, 1]);
        rowmatrix_t = dofE_sca*ones(1,ndofE_sca);        
        rowe_t = reshape(rowmatrix_t,[ndofE_sca^2, 1]);
        cole_t = reshape(rowmatrix_t',[ndofE_sca^2, 1]);
        
        Cttvec(count_t-ndofE_sca^2:count_t-1) = Ctte;
        rowt(count_t-ndofE_sca^2:count_t-1) = rowe_t;
        colt(count_t-ndofE_sca^2:count_t-1) = cole_t;
        
    % Coupling matrix
        count_ut = count_ut + ndofE_mech*ndofE_sca;
        Cute = reshape(Cute, [ndofE_mech*ndofE_sca, 1]);       
        rowe_ut = reshape(dofE_mech*ones(1,ndofE_sca),[ndofE_mech*ndofE_sca, 1]);
        cole_ut = reshape(ones(ndofE_mech,1)*dofE_sca',[ndofE_mech*ndofE_sca, 1]);
        
        Cutvec(count_ut-ndofE_mech*ndofE_sca:count_ut-1) = Cute;
        rowut(count_ut-ndofE_mech*ndofE_sca:count_ut-1) = rowe_ut;
        colut(count_ut-ndofE_mech*ndofE_sca:count_ut-1) = cole_ut;
end

% sparse heat capacity matrix
Ctt = sparse(rowt, colt, Cttvec, Mesh.nDOF, Mesh.nDOF);

% sparse coupling matrix 
Cut = sparse(rowut, colut, Cutvec, Mesh.nDOF, Mesh.nDOF);
Ctu = T0*Cut';

% sparse stiffness matrix
C = Ctt + Ctu;

end