function [disp_err, disp_err_stored] = ManufacturedSolution_Dynamic_point_check(d1, d2, d3, Mesh1, Mesh2, Mesh3, Control)
%MANUFACTUREDSOLUTION_CHECK Calculates the convergence rates
%   [m_L2, m_e] = ManufacturedSolution_check(d1, d2, d3, s1, s2, s3, Mesh1,
%   Mesh2, Mesh3) calculates the rates of convergence of the L2 error norm
%   and the energy norm for using a manufactured solution
%
%   Evaluates the solution of the central node of the mesh
%
%   ----------------------------------------------------------
%   Input
%   ----------------------------------------------------------
%   d1, d2, d3:             Displacement vectors from 3 runs with increasing mesh
%                           refinement
%   Mesh1, Mesh2, Mesh3:    Mesh data structure from 3 runs with increasing mesh
%                           refinement
%   1 - coarsest mesh, 2 - medium mesh, 3 - finest mesh

% Manufactured solution
% u := (x1, x2, t) -> -1/1000*sin(1/2*pi*x1)*sin(1/2*pi*x2)*sin(2*pi*t)
% v := (x1, x2, t) -> 1/1000*cos(1/2*pi*x1)*cos(1/2*pi*x2)*cos(2*pi*t)

global n_steps plot_on

load('ManufacturedSolution_Dynamic_check_data.mat') % Load stored data from a previous run

% x and y are the points of analysis
% n_stesp are the number of load stesp

n = 1;  %This variable ensures that it is the center node

time = Control.StartTime:Control.TimeStep:Control.EndTime;
d_diff = [0,0,0];
d_diff_exact = d_diff;
d_data = zeros(2, length( d_coarse_center),3);

for sim = 1:3
    switch sim
        case 1
            Mesh = Mesh1;
            node = (Mesh.nn-1)/2+n;
            dof = Mesh.DOF(node,:);
            d = d1(dof,:);
            d_diff(sim) = norm([d(1,:) - d_coarse_center(1,:), d(2,:) - d_coarse_center(2,:) ]) / ...
                            norm([d_coarse_center(1,:), d_coarse_center(2,:)]);
            d_data(:,:,sim) = d;
        case 2
            Mesh = Mesh2;
            node = (Mesh.nn-1)/2+n;
            dof = Mesh.DOF(node,:);
            d = d2(dof,:);
            d_diff(sim) = norm([d(1,:) - d_fine_center(1,:), d(2,:) - d_fine_center(2,:) ]) / ...
                            norm([d_fine_center(1,:), d_fine_center(2,:)]);
            d_data(:,:,sim) = d;
        case 3
            Mesh = Mesh3;
            node = (Mesh.nn-1)/2+n;
            dof = Mesh.DOF(node,:);
            d = d3(dof,:);
            d_diff(sim) = norm([d(1,:) - d_finer_center(1,:), d(2,:) - d_finer_center(2,:) ]) / ...
                            norm([d_finer_center(1,:), d_finer_center(2,:)]);
            d_data(:,:,sim) = d;
    end

    if plot_on
        figure(1), hold on 
        plot(time, d(1,:), 'r')
        DisplayName = sprintf('Numerical solution, Mesh %d', sim);
        title(DisplayName)

        figure(2), hold on 
        plot(time, d(2,:), 'r' )
        DisplayName = sprintf('Numerical solution, Mesh %d', sim);
        title(DisplayName)
    end
    
end
    
    x = Mesh.x(node,1);
    y = Mesh.x(node,2);
    
    % Calculate exact solutions at point
    d_exact = zeros(2,n_steps+1);
    d_exact(1,:) =  -sin(pi.*x./2).*sin(pi.*y./2).*sin(2.*pi.*time)/1000;
    d_exact(2,:) = cos(pi.*x./2).*cos(pi.*y./2).*cos(2.*pi.*time)/1000;
    
for sim = 1:3
    d_diff_exact(sim) = norm([d_data(1,:,sim) - d_exact(1,:), d_data(2,:,sim) - d_exact(2,:)]) / ...
                    norm([d_exact(1,:), d_exact(2,:)]);
end
     
    disp_err = max(d_diff_exact);
    disp_err_stored = max(d_diff);

        
    if plot_on
        figure(1), hold on 
        plot(time, d_exact(1,:), 'k')
        title('Exact displacement of central node in the X direction of point (0.5,0.5)')
        axis on 
        grid on 
        legend('Location', 'northwest')
        xlabel('time [s]')
        ylabel('x-displacement [m]')

        figure(2), hold on 
        plot(time, d_exact(2,:), 'k')
        title('Exact displacement of central node in the Y direction of point (0.5,0.5)')
        axis on 
        grid on 
        legend('Location', 'northwest')
        xlabel('time [s]')
        ylabel('y-displacement [m]')
    end

end


