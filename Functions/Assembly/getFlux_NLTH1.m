function [strain, stress, gradT, flux] = getFlux_NLTH1(d, Mesh, Material, calc_type, Quad, ~)
%GETFLUX_TH1 Evaluate fluxes from a diffusion problem
%   [gradT, flux] = GETFLUX_NLTH1(d, Mesh, Material) returns one matrix of 
%   temperature gradients and one matrix offluxes. The matrices are of size 
%   dim x ne, in which dim = 1 for 1D elements, 3 for 2D elements, and 6 
%   for 3D elements.
%   
%   Supported input for calc_type:
%   'none' -  does not compute the fluxes
%   'nodal' -  returns matrices of nodal-averaged flux (size dim x nn)
%   'center' - returns the matrices of fluxes at the center of each element
%               (size dim x ne).
%   'L2projection' - returns matrices of L2-projected fluxes at the nodes (size dim x nn).
% 
%   --------------------------------------------------------------------
%   Input
%   --------------------------------------------------------------------
%   d:          Vector of displacements at each DOF (size ndof x 1 in 
%               which ndof is the number of degrees of freedom)
%   dC:         Displacement at centroid of element
% 
%   Mesh:       Structure array with the following fields,
%               .nsd:   Number of spatial dimensions
%               .ne:    Total number of elements in the mesh
%               .nn:    Total number of nodes 
%               .nne:   Vector of number of nodes per element (size nsd x 1)
%               .type:  the toplogical class of finite element; it is in 
%                       the general form 'topology-#of nodes' ie a three 
%                       node triangle is T3 a four node quadralateral is 
%                       Q4 a 4 node tetrahedra is H4 a 27 node brick is 
%                       B27 etc. Presently defined are L2, Q4, and Q9. 
%               .nDOFe: Number of DOFs per element
%               .conn:  Array of element connectivity (size ne x nne)
%               .x:     Array of nodal spatial locations for
%                       undeformed mesh (size nn x nsd)
%               .DOF:   Array of DOF indices (size nn x nsd)
% 
%   Material:   Structure array with the following fields,
%               .E:     Modulus of elasticity
%               .nu:    Poisson's ratio
%               .Dtype: 2D approximation ('PlaneStrain' or 'PlainStress')
%
% 	Quad: 	Structure array with following fields,
%       	.nq:     Number of quadrature points	
%       	.W:      Vector of quadrature weights (size nq x 1)      
%       	.Q:      Vector of quadrature points	(size nq x nsd)
%       	.Nq:     Cell array (size nq x 1) with shape functions  
%       	         evaluated at each quadrature point
%       	.dNdxiq: Cell array (size nq x 1) with derivative of shape 
%       	         functions w.r.t. parent coordinates evaluated at 
%       	         each quadrature point
%       	.Nv:     Cell array (size nq x 1) with shape functions 
%       	         evaluated at each quadrature point in Voigt form

% initialize D matrix file pointer
[~,DMatrix_functn] = fileparts(Material.ConstitutiveLawFile);

% void parameters (used in elasticity problems)
strain = [];
stress = [];

if nargin < 4
    calc_type = 'center';
end

% Specify dimension of the flux matrix
dim = Mesh.nsd;

if strcmp(calc_type, 'none')
    flux = zeros(dim, Mesh.ne);
    gradT = zeros(dim, Mesh.ne);
else
    % Specify type of strain/stress matrix/cell
    switch calc_type
        case 'nodal'
            flux = zeros(dim, Mesh.nn);
            gradT = zeros(dim, Mesh.nn);
            count = zeros(dim, Mesh.nn);
        case 'center'
            flux = zeros(dim, Mesh.ne);
            gradT = zeros(dim, Mesh.ne);
        case 'L2projection'
            vec_size = Mesh.ne*Mesh.nne;
            A = zeros(vec_size,1); % matrix of integrals of shape functions - vectorized
            row = zeros(vec_size,1); % vector of row indices
            col = zeros(vec_size,1); % vector of column indices
            count = 1;

            be_gradT = zeros(Mesh.nn, dim); % column vectors of gradients
            be_flux = zeros(Mesh.nn, dim); % column vectors of fluxes 
    end

    % Loop through all elements
    for e = 1:Mesh.ne

            %% Element variables
                % nodal ids of the element's nodes
                enodes = Mesh.conn(e,:);    
                % global coordinates of the element's nodes
                xI = Mesh.x(enodes,:);      
                % DOFs of element nodes
                dofE = Mesh.DOF(enodes,:);
                dofE = reshape(dofE',Mesh.nDOFe,[]);
                % local displacement vector
                de = d(dofE);
            %% Constitutive matrix
                nMat = Mesh.MatList(e); % element material type
        
        switch calc_type
            case 'center'
                xi = zeros(1,Mesh.nsd);
                [N, dNdxi] = lagrange_basis(Mesh.type, xi, Mesh.nsd);
                Je = dNdxi'*xI;
                B = Je\(dNdxi');
                dC = N'*de;
                gradT(:, e) = B*de;
                [D,~] = feval(DMatrix_functn, nMat, Material, Mesh, dC);
                flux(:, e) = D*gradT(:, e);
            case 'nodal'
                % initialize grad/flux element
                gradT_e = zeros(dim, Mesh.nne);
                flux_e = zeros(dim, Mesh.nne);   

                % loop through all nodes and calculate nodal strains/stresses
                for n = 1:Mesh.nne

                    % node point in parent coordinate
                    xi = getXI(n, Mesh.type);  

                    % Shape function and derivatives in parent coordinates
                    [~, dNdxi] = lagrange_basis(Mesh.type, xi, Mesh.nsd);

                    % Jacobian of the transformation between parent and global 
                    % coordinates
                    Je = dNdxi'*xI;

                    % derivative of shape function in physical coordinates 
                    % (tensor form)
                    dNdxi = dNdxi';
                    B = Je\dNdxi;

                    % Evaluate Constutive Law
                    [D,~] = feval(DMatrix_functn, nMat, Material, Mesh, d(n));

                    % Gradient and Flux at the Nodes
                    gradT_e(:,n) = B*de;
                    flux_e(:,n) = D*gradT_e(:,n);
                end
                % Add to global strains
                    flux(:,enodes) = flux(:,enodes) + flux_e;
                    gradT(:,enodes) = gradT(:,enodes) + gradT_e;
                    count(:,enodes) = count(:,enodes) + 1;
            case 'L2projection'
                % initialize elemental matrices and vectors
                A_e = zeros(Mesh.nne, Mesh.nne);
                gradT_e = zeros(Mesh.nne,dim);    % column vector of gradients
                flux_e = zeros(Mesh.nne,dim);    % column vector of fluxes

                % loop through all quadrature points
                for q = 1:Quad.nq     

                    % Shape functions and derivatives in parent coordinates
                    N = Quad.Nq{q}';
                    dNdxi = Quad.dNdxiq{q};

                    % Jacobian of the transformation between parent and global 
                    % coordinates
                    Je = dNdxi'*xI;

                    % determinant of the Jacobian
                    dJe = det(Je);

                    % derivative of shape function in physical coordinates 
                    % (tensor form)
                    dNdxi = dNdxi';
                    B = Je\dNdxi;

                    % Compute Temperature at Gaussian Point
                    dG = N*de;

                    % Evaluate Constutive Law at Gaussian Point
                    [D,~] = feval(DMatrix_functn, nMat, Material, Mesh, dG);

                    % calculate stress and strain at quadrature point
                    gradT_q = B*de;
                    flux_q = D*gradT_q;

                    % Element level integral of L2-projections
                    A_e = A_e + N'*N*Quad.W(q)*dJe;
                    
                    for i = 1:dim
                       gradT_e(:,i) = gradT_e(:,i) +  N'*Quad.W(q)*dJe*gradT_q(i);
                       flux_e(:,i) = flux_e(:,i) +  N'*Quad.W(q)*dJe*flux_q(i);
                    end
                end

                % Add to global matrices and vectors
                % Form the vectorized A matrix
                nAe = Mesh.nne^2; % number of entries in A matrix
                count = count + nAe;
                A_e = reshape(A_e, [nAe,1]);
                rowmatrix = enodes'*ones(1,Mesh.nne);
                row_e = reshape(rowmatrix, [nAe,1]);
                col_e = reshape(rowmatrix',[nAe,1]);

                A(count-nAe:count-1) = A_e;
                row(count-nAe:count-1) = row_e;
                col(count-nAe:count-1) = col_e;

                % Add element vectors to global vectors
                be_gradT(enodes,:) = be_gradT(enodes,:) + gradT_e;
                be_flux(enodes,:) = be_flux(enodes,:) + flux_e;
         end
    end

    switch calc_type
        case 'nodal' 
           % For nodal strains, divide by count to get the average
           gradT = gradT./count;
           flux = flux./count;
        case 'L2projection'
            % For L2 projection, solve the system of equations and assemble the
            % nodal matrix of stresses and strains
            A = sparse(row, col, A, Mesh.nn, Mesh.nn);
            eL2_gradT = zeros(size(be_gradT));
            eL2_flux = zeros(size(be_flux));
            
            for i = 1:dim
               eL2_gradT(:,i) = A\be_gradT(:,i);
               eL2_flux(:,i) = A\be_flux(:,i);
            end
            
            gradT = eL2_gradT';
            flux = eL2_flux';
    end
end

end