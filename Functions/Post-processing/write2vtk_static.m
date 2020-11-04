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
%               .contour:       Nodal averaging for discontinous variables
%                               ('none', 'nodal')
%   fixedDOF:   Row vector containing fixed degrees of freedom
%   d:          Column vector of displacements (size ndof x 1 where ndof 
%               is the number of degrees of freedom)
%   strain:     Matrix of nodal strains (size dim x nn in which dim is 
%               1 for 1D elements, 3 for 2D elements, and 6 for 3D elements
%               and nn is the number of nodes)
%   stress:     Matrix of nodal stresses (size dim x nn in which dim is 
%               1 for 1D elements, 3 for 2D elements, and 6 for 3D elements
%               and nn is the number of nodes)
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
    scalardata(1).name = 'nID';
    scalardata(1).data = (1:Mesh.nn)';
    scalardata(1).type = 'int';

    switch Mesh.nsd
        case 1
            scalardata(end+1).name = 'UX';
            scalardata(end).data = d(Mesh.xdofs);
            scalardata(end).type = 'float';

            if ~strcmp(Control.contour,'none')
                scalardata(end+1).name = 'exx';
                scalardata(end).data = strain(1,:)';
                scalardata(end).type = 'float';
                scalardata(end+1).name = 'sxx';
                scalardata(end).data = stress(1,:)';
                scalardata(end).type = 'float';
                
            end
            
        case 2
            scalardata(end+1).name = 'dofX';
            scalardata(end).data = Mesh.xdofs';
            scalardata(end).type = 'float';
            
            scalardata(end+1).name = 'dofY';
            scalardata(end).data = Mesh.ydofs';
            scalardata(end).type = 'float';

            fixed = zeros(Mesh.nDOF,1);
            fixed(fixedDOF) = 1;
            scalardata(end+1).name = 'fixedX';
            scalardata(end).data = fixed(Mesh.xdofs);
            scalardata(end).type = 'float';

            scalardata(end+1).name = 'fixedY';
            scalardata(end).data = fixed(Mesh.ydofs);
            scalardata(end).type = 'float';

            scalardata(end+1).name = 'UX';
            scalardata(end).data = d(Mesh.xdofs);
            scalardata(end).type = 'float';

            scalardata(end+1).name = 'UY';
            scalardata(end).data = d(Mesh.ydofs);
            scalardata(end).type = 'float';

            scalardata(end+1).name = 'Fext_x';
            scalardata(end).data = Fext(Mesh.xdofs);
            scalardata(end).type = 'float';

            scalardata(end+1).name = 'Fext_y';
            scalardata(end).data = Fext(Mesh.ydofs);
            scalardata(end).type = 'float';
            
            Fmag = sqrt(Fext(Mesh.xdofs).^2 + Fext(Mesh.ydofs).^2);

            scalardata(end+1).name = 'Fext';
            scalardata(end).data = Fmag;
            scalardata(end).type = 'float';

            scalardata(end+1).name = 'Fint_x';
            scalardata(end).data = Fint(Mesh.xdofs);
            scalardata(end).type = 'float';

            scalardata(end+1).name = 'Fint_y';
            scalardata(end).data = Fint(Mesh.ydofs);
            scalardata(end).type = 'float';

            Fintmag = sqrt(Fint(Mesh.xdofs).^2 + Fint(Mesh.ydofs).^2);
            
            scalardata(end+1).name = 'Fint';
            scalardata(end).data = Fintmag;
            scalardata(end).type = 'float';
            
            scalardata(end+1).name = 'Rx';
            scalardata(end).data = R(Mesh.xdofs);
            scalardata(end).type = 'float';

            scalardata(end+1).name = 'Ry';
            scalardata(end).data = R(Mesh.ydofs);
            scalardata(end).type = 'float';
            
            Rmag = sqrt(R(Mesh.xdofs).^2 + R(Mesh.ydofs).^2);
            
            scalardata(end+1).name = 'R';
            scalardata(end).data = Rmag;
            scalardata(end).type = 'float';
            
            if isfield(Control,'contour') && ~strcmp(Control.contour,'none')

                scalardata(end+1).name = 'exx';
                scalardata(end).data = strain(1,:)';
                scalardata(end).type = 'float';

                scalardata(end+1).name = 'eyy';
                scalardata(end).data = strain(2,:)';
                scalardata(end).type = 'float';

                scalardata(end+1).name = 'exy';
                scalardata(end).data = strain(3,:)';
                scalardata(end).type = 'float';

                scalardata(end+1).name = 'Sxx';
                scalardata(end).data = stress(1,:)';
                scalardata(end).type = 'float';

                scalardata(end+1).name = 'Syy';
                scalardata(end).data = stress(2,:)';
                scalardata(end).type = 'float';

                scalardata(end+1).name = 'Sxy';
                scalardata(end).data = stress(3,:)';
                scalardata(end).type = 'float';
            
            end
    
        case 3

            scalardata(end+1).name = 'dofX';
            scalardata(end).data = Mesh.xdofs';
            scalardata(end).type = 'float';
            
            scalardata(end+1).name = 'dofY';
            scalardata(end).data = Mesh.ydofs';
            scalardata(end).type = 'float';

            scalardata(end+1).name = 'dofZ';
            scalardata(end).data = Mesh.zdofs';
            scalardata(end).type = 'float';

            fixed = zeros(Mesh.nDOF,1);
            fixed(fixedDOF) = 1;
            scalardata(end+1).name = 'fixedX';
            scalardata(end+1).name = 'fixedX';
            scalardata(end).data = fixed(Mesh.xdofs);
            scalardata(end).type = 'float';

            scalardata(end+1).name = 'fixedY';
            scalardata(end).data = fixed(Mesh.ydofs);
            scalardata(end).type = 'float';

            scalardata(end+1).name = 'fixedZ';
            scalardata(end).data = fixed(Mesh.zdofs);
            scalardata(end).type = 'float';

            scalardata(end+1).name = 'UX';
            scalardata(end).data = d(Mesh.xdofs);
            scalardata(end).type = 'float';

            scalardata(end+1).name = 'UY';
            scalardata(end).data = d(Mesh.ydofs);
            scalardata(end).type = 'float';

            scalardata(end+1).name = 'UZ';
            scalardata(end).data = d(Mesh.zdofs);
            scalardata(end).type = 'float';

            scalardata(end+1).name = 'Fext_x';
            scalardata(end).data = Fext(Mesh.xdofs);
            scalardata(end).type = 'float';

            scalardata(end+1).name = 'Fext_y';
            scalardata(end).data = Fext(Mesh.ydofs);
            scalardata(end).type = 'float';
 
            scalardata(end+1).name = 'Fext_z';
            scalardata(end).data = Fext(Mesh.zdofs);
            scalardata(end).type = 'float';
            
            Fmag = sqrt(Fext(Mesh.xdofs).^2 + Fext(Mesh.ydofs).^2 + Fext(Mesh.zdofs).^2);

            scalardata(end+1).name = 'Fext';
            scalardata(end).data = Fmag;
            scalardata(end).type = 'float';

            scalardata(end+1).name = 'Fint_x';
            scalardata(end).data = Fint(Mesh.xdofs);
            scalardata(end).type = 'float';

            scalardata(end+1).name = 'Fint_y';
            scalardata(end).data = Fint(Mesh.ydofs);
            scalardata(end).type = 'float';

            scalardata(end+1).name = 'Fint_z';
            scalardata(end).data = Fint(Mesh.zdofs);
            scalardata(end).type = 'float';

            Fintmag = sqrt(Fint(Mesh.xdofs).^2 + Fint(Mesh.ydofs).^2 + Fint(Mesh.zdofs).^2);
            
            scalardata(end+1).name = 'Fint';
            scalardata(end).data = Fintmag;
            scalardata(end).type = 'float';
            
            scalardata(end+1).name = 'Rx';
            scalardata(end).data = R(Mesh.xdofs);
            scalardata(end).type = 'float';

            scalardata(end+1).name = 'Ry';
            scalardata(end).data = R(Mesh.ydofs);
            scalardata(end).type = 'float';

            scalardata(end+1).name = 'Rz';
            scalardata(end).data = R(Mesh.zdofs);
            scalardata(end).type = 'float';
            
            Rmag = sqrt(R(Mesh.xdofs).^2 + R(Mesh.ydofs).^2 + R(Mesh.zdofs).^2);
            
            scalardata(end+1).name = 'R';
            scalardata(end).data = Rmag;
            scalardata(end).type = 'float';

            if isfield(Control,'contour') && ~strcmp(Control.contour,'none')

                scalardata(end+1).name = 'exx';
                scalardata(end).data = strain(1,:)';
                scalardata(end).type = 'float';

                scalardata(end+1).name = 'eyy';
                scalardata(end).data = strain(2,:)';
                scalardata(end).type = 'float';

                scalardata(end+1).name = 'ezz';
                scalardata(end).data = strain(3,:)';
                scalardata(end).type = 'float';

                scalardata(end+1).name = 'eyz';
                scalardata(end).data = strain(4,:)';
                scalardata(end).type = 'float';

                scalardata(end+1).name = 'exz';
                scalardata(end).data = strain(5,:)';
                scalardata(end).type = 'float';

                scalardata(end+1).name = 'exy';
                scalardata(end).data = strain(6,:)';
                scalardata(end).type = 'float';

                scalardata(end+1).name = 'Sxx';
                scalardata(end).data = stress(1,:)';
                scalardata(end).type = 'float';

                scalardata(end+1).name = 'Syy';
                scalardata(end).data = stress(2,:)';
                scalardata(end).type = 'float';

                scalardata(end+1).name = 'Szz';
                scalardata(end).data = stress(3,:)';
                scalardata(end).type = 'float';

                scalardata(end+1).name = 'Syz';
                scalardata(end).data = stress(4,:)';
                scalardata(end).type = 'float';
                
                scalardata(end+1).name = 'Sxz';
                scalardata(end).data = stress(5,:)';
                scalardata(end).type = 'float';
 
                scalardata(end+1).name = 'Sxy';
                scalardata(end).data = stress(6,:)';
                scalardata(end).type = 'float';
 
            end

    end

    celldata(1).name = 'eID';
    celldata(1).data = (1:Mesh.ne)';
    celldata(1).type = 'int';

%% Write to vtk
    WriteMesh2VTK(filename1, description, Mesh.x, ...
                    Mesh.conn, scalardata, celldata);
    WriteMesh2VTK(filename2, description, deformedshape, ...
            Mesh.conn, scalardata, celldata);