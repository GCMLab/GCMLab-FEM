function [stress_rt_L2, stress_rt_nodal, stress_rt_an] = BoreholePressure_check(Mesh,d,p,stress_L2, stress_nodal)
%BOREHOLEPRESSURE_CHECK Calculates the error of the exact solution of a
%borehole subjected to a internal normal pressure

%   Analytical solution obtained from book:
%   Analytic Methods in Geomechanics by Kam-tim Chau (2013)
%   Link: https://www.taylorfrancis.com/books/mono/10.1201/9781315275277/analytic-methods-geomechanics-kam-tim-chau
%
%   ----------------------------------------------------------
%   Input
%   ----------------------------------------------------------
%   Mesh :                  Mesh data structure
%   d :                     displacements
%   p :                     Pressure
%   stress_L2:              Stress solution using L2 projection
%   stress_nodal:           Stress solution using nodal computation
%   ----------------------------------------------------------
%   Output
%   ----------------------------------------------------------
%   stress_rt_L2:           Stress in cylindrical coordinates using L2
%                           projection
%   stress_rt_nodal:        Stress in cylindrical coordinates using nodal
%                           computation
%   stress_rt_an:           Analytical solution of stress in cylindrical
%                           coordinates

% The displacement vector is obtained as
% u = - p R^2 / r^2 e_r 
%   where R is the radious of the hole
R = 0.5; % Taken from mesh
%   p is the pressure 
%   e_r is the unit vector in the radial direction. 

% Get cylindrical coordinates of all nodes
r_m = (Mesh.x(:,1).^2 + Mesh.x(:,2).^2).^0.5; % vector to store radial position of all points
theta_m = atan(Mesh.x(:,2)./Mesh.x(:,1)); % vector to store theta angle of all points

% Convert cartesian solution to cylindrical solution
dx = d(1:2:end); %displacements on x
dy = d(2:2:end); %displacements on y 


% Convert displacement and stress in cartesian coordinates to cylindrical
disp_rt = zeros(Mesh.nn,2); %displacement
stress_rt_L2 = zeros(Mesh.nn,3); %stress (only in e_rr and e_tt)
stress_rt_nodal = zeros(Mesh.nn,3); %stress (only in e_rr and e_tt)
for i = 1:Mesh.nn
    % Get transformation matrix
    T = [cos(theta_m(i)), sin(theta_m(i)); -sin(theta_m(i)), cos(theta_m(i))]; 
    % Displacement vector in cylindrical coordinates
    disp_rt(i,:) = (T*[dx(i); dy(i)])';
    
    % L2 projection
    % stress in cartesian coordinates
    scar = [stress_L2(1,i), stress_L2(3,i); stress_L2(3,i), stress_L2(2,i)];
    % stress in cylindrical coordinates
    stress_i = (T*scar*T')';  
    stress_rt_L2(i,:) = [stress_i(1,1), stress_i(2,2), stress_i(1,2)];
    
    % Nodal
    % stress in cartesian coordinates
    scar = [stress_nodal(1,i), stress_nodal(3,i); stress_nodal(3,i), stress_nodal(2,i)];
    % stress in cylindrical coordinates
    stress_i = (T*scar*T')';  
    stress_rt_nodal(i,:) = [stress_i(1,1), stress_i(2,2), stress_i(1,2)];
end

% Get analytical solution of stress in each node
stress_rt_an = zeros(Mesh.nn,2);
stress_rt_an(:,1) = p.*R^2./r_m.^2;
stress_rt_an(:,2) =  (-p.*R^2./r_m.^2);




        

    



