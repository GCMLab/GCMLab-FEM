function write2vtk_static(config_name, vtk_dir, Mesh, Control, ...
                        fixedDOF, d, strain, stress, Fint, Fext)
%WRITE2VTK_STATIC Exports results to VTK file 
%   WRITE2VTK_STATIC(Mesh, Control, fixedDOF, d, strain, stress, Fint, Fext) 
%   produces vtk files that can be opened in Paraview to visualize
%   the simulation results. The mesh, boundary conditions, and 
%   simulation parameters are described in the structure arrays
%   Mesh, fixedDOF, and Control, respectively. The simulation results are 
%   the vectors d, Fint, and Fext of size ndof x 1, containing the 
%   displacement, internal forces, and external forces, respectively. 
% 
%   --------------------------------------------------------------------
%   Input
%   --------------------------------------------------------------------
%   Mesh:       Structure array with the following fields, 
%               .x:     Array of nodal spatial locations for
%                       undeformed mesh (size nn x nsd)
%               .conn:  Array of element connectivity (size ne x nne)
%               .xdofs: DOFs in the x-direction
%               .ydofs: DOFs in the y-direction
%               .zdofs: DOFs in the z-direction
%               .nn:    Total number of nodes 
%               .nsd:   Number of spatial dimensions
%               .nDOF:  Total number of DOFs
%               .ne:    Total number of elements in the mesh
%   Control:    Structure array with the following fields,
%               .config_name:   Name of the configuration file
%               .vtk_dir:       Directory where VTK files are stored
%               .MagCoef:       Displacement magnification coefficient
%                               for visualization
%               .stress_calc:   Calculation of values for discontinous variables
%                               ('none', 'nodal', 'center')
%   fixedDOF:   Row vector containing fixed degrees of freedom
%   d:          Column vector of displacements (size ndof x 1 where ndof 
%               is the number of degrees of freedom)
%   strain:     if Control.stress_calc = 'nodal': Matrix of nodal strains 
%                           (size dim x nn in which dim is 1 for 1D 
%                           elements, 3 for 2D elements, and 6 for 3D 
%                           elements and nn is the number of nodes)
%               if Control.stress_calc = 'center': Matrix of element 
%                           strains (calculated at the center of each 
%                           element); size dim x ne in which ne is the 
%                           number of elements
%   stress:     Matrix of stresses (same size as strain matrix)
%   Fint:       Column vector of internal forces (size ndof x 1)
%   Fext:       Column vector of external forces (size ndof x 1)

%% Define variables
    description = config_name;
    filename1 = fullfile(vtk_dir, [config_name '.vtk.0']); 
    filename2 = fullfile(vtk_dir, [config_name '.vtk.1']);

    deformedshape = Mesh.x + ...
                Control.MagCoef*[d(Mesh.xdofs),d(Mesh.ydofs),d(Mesh.zdofs)];

    R = Fext - Fint;

%% Data
    nodedata(1).name = 'nID';
    nodedata(1).data = (1:Mesh.nn)';
    nodedata(1).type = 'int';

    elementdata(1).name = 'eID';
    elementdata(1).data = (1:Mesh.ne)';
    elementdata(1).type = 'int';

    nodedata(end+1).name = 'dof';
    nodedata(end).data = Mesh.DOF;
    nodedata(end).type = 'int';

    fixed = zeros(Mesh.nDOF,1);
    fixed(fixedDOF) = 1;
    nodedata(end+1).name = 'fixed';
    nodedata(end).data = fixed(Mesh.DOF);
    nodedata(end).type = 'int';

    switch Mesh.nsd
        case 1

            nodedata(end+1).name = 'U';
            nodedata(end).data = d(Mesh.DOF);
            nodedata(end).type = 'float';

            nodedata(end+1).name = 'Fext';
            nodedata(end).data = Fext(Mesh.DOF);
            nodedata(end).type = 'float';

            nodedata(end+1).name = 'Fint';
            nodedata(end).data = Fext(Mesh.DOF);
            nodedata(end).type = 'float';

            nodedata(end+1).name = 'R';
            nodedata(end).data = R(Mesh.DOF);
            nodedata(end).type = 'float';

            if strcmp(Control.stress_calc,'nodal')
                nodedata(end+1).name = 'exx';
                nodedata(end).data = strain';
                nodedata(end).type = 'float';
                nodedata(end+1).name = 'sxx';
                nodedata(end).data = stress';
                nodedata(end).type = 'float';                
            elseif strcmp(Control.stress_calc, 'center')
                elementdata(end+1).name = 'exx';
                elementdata(end).data = strain';
                elementdata(end).type = 'float';
            end
            
        case 2
            nodedata(end+1).name = 'U';
            nodedata(end).data = [d(Mesh.DOF), zeros(size(Mesh.xdofs'))];
            nodedata(end).type = 'float';

            nodedata(end+1).name = 'Fext';
            nodedata(end).data = [Fext(Mesh.DOF), zeros(size(Mesh.xdofs'))];
            nodedata(end).type = 'float';
            
            nodedata(end+1).name = 'Fint';
            nodedata(end).data = [Fint(Mesh.DOF), zeros(size(Mesh.xdofs'))];
            nodedata(end).type = 'float';

            nodedata(end+1).name = 'R';
            nodedata(end).data = [R(Mesh.DOF), zeros(size(Mesh.xdofs'))];
            nodedata(end).type = 'float';
            
            if strcmp(Control.stress_calc,'nodal')
                nodedata(end+1).name = 'exx';
                nodedata(end).data = strain(1,:)';
                nodedata(end).type = 'float';

                nodedata(end+1).name = 'eyy';
                nodedata(end).data = strain(2,:)';
                nodedata(end).type = 'float';

                nodedata(end+1).name = 'exy';
                nodedata(end).data = strain(3,:)';
                nodedata(end).type = 'float';

                nodedata(end+1).name = 'Sxx';
                nodedata(end).data = stress(1,:)';
                nodedata(end).type = 'float';

                nodedata(end+1).name = 'Syy';
                nodedata(end).data = stress(2,:)';
                nodedata(end).type = 'float';

                nodedata(end+1).name = 'Sxy';
                nodedata(end).data = stress(3,:)';
                nodedata(end).type = 'float';
            
            elseif strcmp(Control.stress_calc, 'center')
                elementdata(end+1).name = 'exx';
                elementdata(end).data = strain(1,:)';
                elementdata(end).type = 'float';

                elementdata(end+1).name = 'eyy';
                elementdata(end).data = strain(2,:)';
                elementdata(end).type = 'float';

                elementdata(end+1).name = 'exy';
                elementdata(end).data = strain(3,:)';
                elementdata(end).type = 'float';

                elementdata(end+1).name = 'Sxx';
                elementdata(end).data = stress(1,:)';
                elementdata(end).type = 'float';

                elementdata(end+1).name = 'Syy';
                elementdata(end).data = stress(2,:)';
                elementdata(end).type = 'float';

                elementdata(end+1).name = 'Sxy';
                elementdata(end).data = stress(3,:)';
                elementdata(end).type = 'float';
            end
    
        case 3

            nodedata(end+1).name = 'U';
            nodedata(end).data = d(Mesh.DOF);
            nodedata(end).type = 'float';

            nodedata(end+1).name = 'Fext';
            nodedata(end).data = Fext(Mesh.DOF);
            nodedata(end).type = 'float';

            nodedata(end+1).name = 'Fint';
            nodedata(end).data = Fext(Mesh.DOF);
            nodedata(end).type = 'float';

            nodedata(end+1).name = 'R';
            nodedata(end).data = R(Mesh.DOF);
            nodedata(end).type = 'float';

            if strcmp(Control.stress_calc,'nodal')

                nodedata(end+1).name = 'exx';
                nodedata(end).data = strain(1,:)';
                nodedata(end).type = 'float';

                nodedata(end+1).name = 'eyy';
                nodedata(end).data = strain(2,:)';
                nodedata(end).type = 'float';

                nodedata(end+1).name = 'ezz';
                nodedata(end).data = strain(3,:)';
                nodedata(end).type = 'float';

                nodedata(end+1).name = 'eyz';
                nodedata(end).data = strain(4,:)';
                nodedata(end).type = 'float';

                nodedata(end+1).name = 'exz';
                nodedata(end).data = strain(5,:)';
                nodedata(end).type = 'float';

                nodedata(end+1).name = 'exy';
                nodedata(end).data = strain(6,:)';
                nodedata(end).type = 'float';

                nodedata(end+1).name = 'Sxx';
                nodedata(end).data = stress(1,:)';
                nodedata(end).type = 'float';

                nodedata(end+1).name = 'Syy';
                nodedata(end).data = stress(2,:)';
                nodedata(end).type = 'float';

                nodedata(end+1).name = 'Szz';
                nodedata(end).data = stress(3,:)';
                nodedata(end).type = 'float';

                nodedata(end+1).name = 'Syz';
                nodedata(end).data = stress(4,:)';
                nodedata(end).type = 'float';
                
                nodedata(end+1).name = 'Sxz';
                nodedata(end).data = stress(5,:)';
                nodedata(end).type = 'float';
 
                nodedata(end+1).name = 'Sxy';
                nodedata(end).data = stress(6,:)';
                nodedata(end).type = 'float';
 
            elseif strcmp(Control.stress_calc, 'center')
                elementdata(end+1).name = 'exx';
                elementdata(end).data = strain(1,:)';
                elementdata(end).type = 'float';

                elementdata(end+1).name = 'eyy';
                elementdata(end).data = strain(2,:)';
                elementdata(end).type = 'float';

                elementdata(end+1).name = 'ezz';
                elementdata(end).data = strain(3,:)';
                elementdata(end).type = 'float';

                elementdata(end+1).name = 'eyz';
                elementdata(end).data = strain(4,:)';
                elementdata(end).type = 'float';

                elementdata(end+1).name = 'exz';
                elementdata(end).data = strain(5,:)';
                elementdata(end).type = 'float';

                elementdata(end+1).name = 'exy';
                elementdata(end).data = strain(6,:)';
                elementdata(end).type = 'float';

                elementdata(end+1).name = 'Sxx';
                elementdata(end).data = stress(1,:)';
                elementdata(end).type = 'float';

                elementdata(end+1).name = 'Syy';
                elementdata(end).data = stress(2,:)';
                elementdata(end).type = 'float';

                elementdata(end+1).name = 'Szz';
                elementdata(end).data = stress(3,:)';
                elementdata(end).type = 'float';

                elementdata(end+1).name = 'Syz';
                elementdata(end).data = stress(4,:)';
                elementdata(end).type = 'float';
                
                elementdata(end+1).name = 'Sxz';
                elementdata(end).data = stress(5,:)';
                elementdata(end).type = 'float';
 
                elementdata(end+1).name = 'Sxy';
                elementdata(end).data = stress(6,:)';
                elementdata(end).type = 'float';
            end

    end

%% Write to vtk
    WriteMesh2VTK(filename1, description, Mesh.x, ...
                    Mesh.conn, nodedata, elementdata);
    WriteMesh2VTK(filename2, description, deformedshape, ...
            Mesh.conn, nodedata, elementdata);