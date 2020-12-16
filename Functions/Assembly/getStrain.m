function [strain, stress] = getStrain(d, Mesh, Material, calc_type)
%GETSTRAIN Evaluate stress and strain
%   strain = GETSTRAIN(d, Mesh, Material) is a cell array of  
%   nodal strains in each element of the mesh. The cell array is of size 
%   dim x ne, in which dim = 1 for 1D elements, 3 for 2D elements, and 6 
%   for 3D elements.
%
%   [strain, stress] = GETSTRAIN(d, Mesh, Material) also returns a cell 
%   array of nodal stresses in each element of the mesh (size dim x ne).
% 
%   [strain, stress] = GETSTRAIN(d, Mesh, Material, 'nodal') returns 
%   matrices of nodal-averaged strains (size dim x nn).
% 
%   [strain, stress] = GETSTRAIN(d, Mesh, Material, 'center') returns 
%   matrices of strains computed at the center of each element 
%   (size dim x ne).
% 
%   --------------------------------------------------------------------
%   Input
%   --------------------------------------------------------------------
%   d:          Vector of displacements at each DOF (size ndof x 1 in 
%               which ndof is the number of degrees of freedom)
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

if nargin < 4
    calc_type = 'none';
end

% Specify dimension of the strain/stress matrix
switch Mesh.nsd
    case 1
        dim = 1;
    case 2
        dim = 3;
    case 3
        dim = 6;
end

% TODO: Add least squares stress/strain projection

% Specify type of strain/stress matrix/cell
switch calc_type
    case 'none'
        strain = cell(dim,Mesh.ne);
        stress = cell(dim,Mesh.ne);
    case 'nodal'
        strain = zeros(dim, Mesh.nn);
        stress = zeros(dim, Mesh.nn);
        count = zeros(dim, Mesh.nn);
    case 'center'
        strain = zeros(dim, Mesh.ne);
        stress = zeros(dim, Mesh.ne);
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
  
    if strcmp(calc_type, 'center')
        xi = zeros(1,Mesh.nsd);
        [N, dNdxi] = lagrange_basis(Mesh.type, xi, Mesh.nsd);
        Xi = xI'*N;
        Je = dNdxi'*xI;
        B = Je\(dNdxi');
        Bv = getBv(B', Mesh.nsd);
        D = getD(Material.E(Xi), Material.nu(Xi), Mesh.nsd, Material.Dtype);
        strain(:, e) = Bv'*de;
        stress(:, e) = D*strain(:, e);
    else
        % initialize strain element
        strain_e = zeros(dim, Mesh.nne);   
        stress_e = zeros(dim, Mesh.nne);
        
        % loop through all nodes and calculate nodal strains/stresses
        for n = 1:Mesh.nne

            % node point in parent coordinate
            xi = getXI(n, Mesh.type);  

            % Shape function and derivatives in parent coordinates
            [N, dNdxi] = lagrange_basis(Mesh.type, xi, Mesh.nsd);

            % quadrature point in physical coordinates
            Xi = xI'*N;

            % Jacobian of the transformation between parent and global 
            % coordinates
            Je = dNdxi'*xI;
           
            % derivative of shape function in physical coordinates 
            % (tensor form)
            dNdxi = dNdxi';
            B = Je\dNdxi;

            % convert B matrix to Voigt form
            Bv = getBv(B', Mesh.nsd);

            D = getD(Material.E(Xi), Material.nu(Xi), Mesh.nsd, Material.Dtype);    
            
            % strain_e = [strainx_n1  strainx_n2...;
                         %strainy_n1  strainy_n2...;
                         %strainxy_n1 strainxy_n2...];
            strain_e(:,n) = Bv'*de;
            stress_e(:,n) = D*strain_e(:,n);
            
        end
            
        %% Add to global strains
        switch calc_type
            case 'none'
                for s = 1:size(strain_e,1)
                    strain{s,e} = strain_e(s,:);
                    stress{s,e} = stress_e(s,:);
                end
            case 'nodal'
                strain(:,enodes) = strain(:,enodes) + strain_e;
                stress(:,enodes) = stress(:,enodes) + stress_e;
                count(:,enodes) = count(:,enodes) + 1;
        end
        
    end
end

%% For nodal strains, divide by count to get the average
    if strcmp(calc_type,'nodal')
       strain = strain./count;
       stress = stress./count;
    end

end