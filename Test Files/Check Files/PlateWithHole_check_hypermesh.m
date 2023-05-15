function [error_nodal, error_L2] = PlateWithHole_check_hypermesh(Mesh, stress_nodal, stress_L2)
%PLATEWITHHOLE_CHECK Calculates the error in nodal and L2 projected stresses
%   [error_nodal, error_L2] = PlateWithHole_check(Mesh, stress_nodal, stress_L2)
%   calculates the error in stress along the line x = 0 

%   Adapted from OpenFOAM: https://www.openfoam.com/documentation/tutorial-guide/tutorialse9.php
%
%   ----------------------------------------------------------
%   Input
%   ----------------------------------------------------------
%   Mesh :                  Mesh data structure
%   stress_nodal:           List of nodal stresses determined with nodal
%                           averaging
%   stress_L2:              List of nodal stresses determined with L2-
%                           projection

% The stress along the line x = 0 is given by 
% sigma_xx = { t * ( 1 + R^2/(2y^2) + 3R^4/(2y^4) for y >= R
%            { 0                                  for y < R
% t = 10 kPa
% R = 0.5 m


plot_on = 0; % turn plots on/off - debugging tool

Mesh.left_nodes = Mesh.BC_nE(find(Mesh.BC_E(:,2) == 1),1);

sxx_nodal = stress_nodal(1,Mesh.left_nodes);
sxx_L2    = stress_L2(1,Mesh.left_nodes);

t = 10e3;
R = 0.5;

sxx_exact = @(y) t*(1 + R^2./(2*y.^2) + 3*R^4./(2*y.^4));

if plot_on
    fplot(sxx_exact,[R,max(Mesh.x(:,2))])
    hold on
    plot(Mesh.x(Mesh.left_nodes,2),sxx_nodal,'o')
    plot(Mesh.x(Mesh.left_nodes,2),sxx_L2,'*')
    hold off
    xlabel('y')
    ylabel('\sigma_{xx}')
    legend('Exact','Nodal Average','L2')
end

% Step along line and determine error
y0 = R;

nsteps = length(Mesh.left_nodes) - 1;
L = max(Mesh.x(:,2))-R;
dL = L/nsteps;

err_nodal = 0;
err_L2 = 0;
err_denom = 0;

nq = 2;
[W,Q] = quadrature(nq,'GAUSS',1);

tol = 1e-6; % meshing tolerance

for i = 1:nsteps
    y1 = y0;
    y2 = y0 + dL;
    
    node1 = find(abs((Mesh.x(Mesh.left_nodes,2) - y1)) < tol);
    node2 = find(abs((Mesh.x(Mesh.left_nodes,2) - y2)) < tol);
    for q = 1:nq
        N = 0.5*[1-Q(q), 1+Q(q)];
        
        sxx_nodal_q = N*[sxx_nodal(node1); sxx_nodal(node2)];
        sxx_L2_q = N*[sxx_L2(node1); sxx_L2(node2)];
        
        yq = N*[y1; y2];
        
        err_nodal = err_nodal + (sxx_nodal_q - sxx_exact(yq))^2*W(q)*dL/2;
        err_L2 = err_L2 + (sxx_L2_q - sxx_exact(yq))^2*W(q)*dL/2;
        err_denom = err_denom + sxx_exact(yq)^2 * W(q) *dL/2;
    end
    
    y0 = y0 + dL;
end

error_nodal = sqrt(err_nodal/err_denom);
error_L2 = sqrt(err_L2/err_denom);

end

