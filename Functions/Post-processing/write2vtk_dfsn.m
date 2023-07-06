function write2vtk_dfsn(config_name, vtk_dir, Mesh, Control, ...
                        fixedDOF, d, flux, ~, Fint, Fext, timestep)
%WRITE2VTK_DFSN Exports results to VTK file 
%   WRITE2VTK_DFSN(Mesh, Control, fixedDOF, d, flux, ~, Fint, Fext) 
%   produces vtk files that can be opened in Paraview to visualize
%   the simulation results. The mesh, boundary conditions, and 
%   simulation parameters are described in the structure arrays
%   Mesh, fixedDOF, and Control, respectively. The simulation results is
%   the vector d, and the fluxes. 
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
%               .stress_calc:   Calculation of values for discontinous variables
%                               ('none', 'nodal', 'center')
%   fixedDOF:   Row vector containing fixed degrees of freedom
%   d:          Column vector of displacements (size ndof x 1 where ndof 
%               is the number of degrees of freedom)
%   flux:       Fluxes either nodal or at the center of the elements.
%   Fint:       Column vector of internal forces (size ndof x 1)
%   Fext:       Column vector of external forces (size ndof x 1)

%% Define variables
    description = config_name;
    filename = fullfile(vtk_dir, [config_name '.vtk.',num2str(timestep)]); 

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

            nodedata(end+1).name = 'T';
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

            if strcmp(Control.stress_calc,'nodal') || strcmp(Control.stress_calc,'L2projection') 
                nodedata(end+1).name = 'qx';
                nodedata(end).data = flux';
                nodedata(end).type = 'float';

            elseif strcmp(Control.stress_calc, 'center')
                elementdata(end+1).name = 'qx';
                elementdata(end).data = flux';
                elementdata(end).type = 'float';
            end
            
        case 2
            nodedata(end+1).name = 'T';
            nodedata(end).data = [d(Mesh.DOF)];
            nodedata(end).type = 'float';

            nodedata(end+1).name = 'Fext';
            nodedata(end).data = [Fext(Mesh.DOF)];
            nodedata(end).type = 'float';
            
            nodedata(end+1).name = 'Fint';
            nodedata(end).data = [Fint(Mesh.DOF)];
            nodedata(end).type = 'float';

            nodedata(end+1).name = 'R';
            nodedata(end).data = [R(Mesh.DOF)];
            nodedata(end).type = 'float';
            

            if strcmp(Control.stress_calc,'nodal') || strcmp(Control.stress_calc,'L2projection') 
                nodedata(end+1).name = 'q';
                nodedata(end).data = [flux(1,:)',flux(2,:)', zeros(size(Mesh.xdofs'))];
                nodedata(end).type = 'float';
            

            elseif strcmp(Control.stress_calc, 'center')
                elementdata(end+1).name = 'q';
                elementdata(end).data = [flux(1,:)', flux(2,:)', zeros(Mesh.ne,1)];
                elementdata(end).type = 'float';

            end
    
        case 3

            nodedata(end+1).name = 'T';
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

            if strcmp(Control.stress_calc,'nodal') || strcmp(Control.stress_calc,'L2projection') 

                nodedata(end+1).name = 'q';
                nodedata(end).data = [flux(1,:)',flux(2,:)',flux(3,:)'];
                nodedata(end).type = 'float';
                
            elseif strcmp(Control.stress_calc, 'center')

                elementdata(end+1).name = 'q';
                elementdata(end).data = [flux(1,:)',flux(2,:)',flux(3,:)'];
                elementdata(end).type = 'float';


            end

    end

%% Write to vtk
    WriteMesh2VTK(filename, description, Mesh.x, ...
                    Mesh.conn, nodedata, elementdata);
