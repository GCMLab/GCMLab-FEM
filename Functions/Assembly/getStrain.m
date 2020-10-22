function [strain, stress] = getStrain(d, Mesh, Control, Material)
% Acknowledgements: Chris Ladubec, Endrina Rivas

% Specify dimension of the strain/stress matrix
switch Mesh.nsd
    case 1
        dim = 1;
    case 2
        dim = 3;
    case 3
        dim = 6;
end

% Specify type of strain/stress matrix/cell
%   if Control.contour = 'none'
%       strain/stress are cell arrays
%       strain = {strainx in e1} {strainx in e2} ... {strainx in enn}
%                {strainy in e1} {strainy in e2} ... {strainy in enn}
%   if Control.contour = 'nodal'
%           strain = [strainx_n1  strainx_n2  ... strainx_nn;
%                     strainy_n1  strainy_n2  ... strainy_nn;
%                     strainxy_n1 strainxy_n2 ... strainxy_nn];
switch Control.contour
    case 'none'
        strain = cell(dim,Mesh.ne);
        stress = cell(dim,Mesh.ne);
    case 'nodal'
        strain = zeros(dim, Mesh.nn);
        stress = zeros(dim, Mesh.nn);
        count = zeros(dim, Mesh.nn);
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
  
    %% Calculate element strains
        % initialize strain element
        strain_e = zeros(dim, Mesh.nne);   
        stress_e = zeros(dim, Mesh.nne);
        
        % loop through all nodes and calculate nodal strains/stresses
        for n = 1:Mesh.nne

            % node point in parent coordinate
            xi = getXI(n,Mesh.type);  

            % Shape function and derivatives in parent coordinates
            [N, dNdxi] = lagrange_basis(Mesh.type, xi);

            % quadrature point in physical coordinates
            Xi = xI'*N;

            % Jacobian of the transformation between parent and global 
            % coordinates
            Je = dNdxi'*xI;
           
            % derivative of shape function in physical coordinates 
            % (tensor form)
            B = dNdxi/Je;

            % convert B matrix to Voigt form
            Bv = getBv(B, Mesh.nsd);

            D = getD(Xi, Mesh.nsd, Material);
            
            % strain_e = [strainx_n1  strainx_n2...;
                         %strainy_n1  strainy_n2...;
                         %strainxy_n1 strainxy_n2...];
            strain_e(:,n) = Bv'*de;
            stress_e(:,n) = D*strain_e(:,n);
            
        end
            
    %% Add to global strains
        switch Control.contour
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

%% For nodal strains, divide by count to get the average
    if strcmp(Control.contour,'nodal')
       strain = strain./count;
       stress = stress./count;
    end

end