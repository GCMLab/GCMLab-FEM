function write2vtk_THLE1(config_name, vtk_dir, Mesh, Control, ...
                        fixedDOF, d, strain, stress, gradT, flux, Fint, Fext, timestep, BC)
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
    
    fixedDOF_mech = BC.fix_disp_dof;
    fixedDOF_sca = BC.fix_temp_dof;

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

    fixed_mech = zeros(Mesh.nDOF,1);
    fixed_mech(fixedDOF_mech) = 1;
    nodedata(end+1).name = 'fixed_mech';
    nodedata(end).data = fixed_mech(Mesh.DOF_mech);
    nodedata(end).type = 'int';
    
    fixed_sca = zeros(Mesh.nDOF,1);
    fixed_sca(fixedDOF_sca) = 1;
    nodedata(end+1).name = 'fixed_sca';
    nodedata(end).data = fixed_sca(Mesh.DOF_sca);
    nodedata(end).type = 'int';

    switch Mesh.nsd
        case 1

            nodedata(end+1).name = 'U';
            nodedata(end).data = d(Mesh.DOF_mech);
            nodedata(end).type = 'float';
            
            nodedata(end+1).name = 'T';
            nodedata(end).data = d(Mesh.DOF_sca);
            nodedata(end).type = 'float';
            
            nodedata(end+1).name = 'Fext_mech';
            nodedata(end).data = Fext(Mesh.DOF_mech);
            nodedata(end).type = 'float';
            
            nodedata(end+1).name = 'Fext_sca';
            nodedata(end).data = Fext(Mesh.DOF_sca);
            nodedata(end).type = 'float';

            nodedata(end+1).name = 'Fint_mech';
            nodedata(end).data = Fext(Mesh.DOF_mech);
            nodedata(end).type = 'float';

            nodedata(end+1).name = 'Fint_sca';
            nodedata(end).data = Fext(Mesh.DOF_sca);
            nodedata(end).type = 'float';
            
            nodedata(end+1).name = 'R_mech';
            nodedata(end).data = R(Mesh.DOF_mech);
            nodedata(end).type = 'float';
            
            nodedata(end+1).name = 'R_sca';
            nodedata(end).data = R(Mesh.DOF_sca);
            nodedata(end).type = 'float';

            if strcmp(Control.stress_calc,'nodal') || strcmp(Control.stress_calc,'L2projection') 
                nodedata(end+1).name = 'exx';
                nodedata(end).data = strain';
                nodedata(end).type = 'float';
                
                nodedata(end+1).name = 'sxx';
                nodedata(end).data = stress';
                nodedata(end).type = 'float';   
                
                nodedata(end+1).name = 'gradTx';
                nodedata(end).data = gradT';
                nodedata(end).type = 'float';

                nodedata(end+1).name = 'qx';
                nodedata(end).data = flux';
                nodedata(end).type = 'float';

            elseif strcmp(Control.stress_calc, 'center')
                elementdata(end+1).name = 'exx';
                elementdata(end).data = strain';
                elementdata(end).type = 'float';
                
                elementdata(end+1).name = 'sxx';
                elementdata(end).data = stress';
                elementdata(end).type = 'float';
                
                elementdata(end+1).name = 'gradTx';
                elementdata(end).data = gradT';
                elementdata(end).type = 'float';
                
                elementdata(end+1).name = 'qx';
                elementdata(end).data = flux';
                elementdata(end).type = 'float';
            end
            
        case 2
            nodedata(end+1).name = 'U';
            nodedata(end).data = [d(Mesh.DOF_mech), zeros(size(Mesh.xdofs_u'))];
            nodedata(end).type = 'float';
            
            nodedata(end+1).name = 'T';
            nodedata(end).data = d(Mesh.DOF_sca);
            nodedata(end).type = 'float';

            nodedata(end+1).name = 'Fext_mech';
            nodedata(end).data = [Fext(Mesh.DOF_mech), zeros(size(Mesh.xdofs_u'))];
            nodedata(end).type = 'float';
            
            nodedata(end+1).name = 'Fext_sca';
            nodedata(end).data = Fext(Mesh.DOF_sca);
            nodedata(end).type = 'float';
            
            nodedata(end+1).name = 'Fint_mech';
            nodedata(end).data = [Fint(Mesh.DOF_mech), zeros(size(Mesh.xdofs_u'))];
            nodedata(end).type = 'float';
            
            nodedata(end+1).name = 'Fint_sca';
            nodedata(end).data = Fint(Mesh.DOF_sca);
            nodedata(end).type = 'float';

            nodedata(end+1).name = 'R_mech';
            nodedata(end).data = [R(Mesh.DOF_mech), zeros(size(Mesh.xdofs_u'))];
            nodedata(end).type = 'float';
            
            nodedata(end+1).name = 'R_sca';
            nodedata(end).data = R(Mesh.DOF_sca);
            nodedata(end).type = 'float';
            

            if strcmp(Control.stress_calc,'nodal') || strcmp(Control.stress_calc,'L2projection') 
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
                
                nodedata(end+1).name = 'gradT';
                nodedata(end).data = [gradT(1,:)',gradT(2,:)', zeros(size(Mesh.dofs_t'))];
                nodedata(end).type = 'float';
                
                nodedata(end+1).name = 'q';
                nodedata(end).data = [flux(1,:)',flux(2,:)', zeros(size(Mesh.dofs_t'))];
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
                
                elementdata(end+1).name = 'q';
                elementdata(end).data = [flux(1,:)', flux(2,:)', zeros(Mesh.ne,1)];
                elementdata(end).type = 'float';

            end
    
        case 3

            nodedata(end+1).name = 'U';
            nodedata(end).data = d(Mesh.DOF_mech);
            nodedata(end).type = 'float';
            
            nodedata(end+1).name = 'T';
            nodedata(end).data = d(Mesh.DOF_sca);
            nodedata(end).type = 'float';

            nodedata(end+1).name = 'Fext_mech';
            nodedata(end).data = Fext(Mesh.DOF_mech);
            nodedata(end).type = 'float';
            
            nodedata(end+1).name = 'Fext_sca';
            nodedata(end).data = Fext(Mesh.DOF_sca);
            nodedata(end).type = 'float';

            nodedata(end+1).name = 'Fint_mech';
            nodedata(end).data = Fext(Mesh.DOF_mech);
            nodedata(end).type = 'float';
            
            nodedata(end+1).name = 'Fint_sca';
            nodedata(end).data = Fext(Mesh.DOF_sca);
            nodedata(end).type = 'float';

            nodedata(end+1).name = 'R_mech';
            nodedata(end).data = R(Mesh.DOF_mech);
            nodedata(end).type = 'float';
            
            nodedata(end+1).name = 'R_sca';
            nodedata(end).data = R(Mesh.DOF_sca);
            nodedata(end).type = 'float';

            if strcmp(Control.stress_calc,'nodal') || strcmp(Control.stress_calc,'L2projection') 

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
                
                nodedata(end+1).name = 'q';
                nodedata(end).data = [flux(1,:)',flux(2,:)',flux(3,:)'];
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
                
                elementdata(end+1).name = 'q';
                elementdata(end).data = [flux(1,:)',flux(2,:)',flux(3,:)'];
                elementdata(end).type = 'float';
            end
    end

%% Write to vtk
    WriteMesh2VTK(filename, description, Mesh.x, ...
                    Mesh.conn, nodedata, elementdata);
