function [time_er] = DirichletTime_check(stress, Material)
%DIRICHLETTIME_CHECK Calculates the error between FEA and the manufactured
%solution
%   time_er = DIRICHLETTIME_CHECK(stress) calculates the error between FEA 
%   and the manufactured solution
%
%   ----------------------------------------------------------
%   Input
%   ----------------------------------------------------------
%   stress:                  3 x nn x nt data structure of stress at each
%                            timestep
%   ----------------------------------------------------------
%   Output
%   ----------------------------------------------------------
%   time_er:             Error related to time

% The manufactured solution is
% ux = y*sin(Omega1*t)
% uy = x*sin(Omega2*t)
% sigma_xx = 0
% sigma_yy = 0
% sigma_xy = E/2/(1+nu)* (sin(Omega1*t) + sin(Omega2*t))


global Omega1 Omega2

% Calculate analytical stress
t = linspace(0,2*pi,51);
sigma_xy_exact = Material.Prop(1).E0/2/(1+Material.Prop(1).nu)*(sin(Omega1*t) + sin(Omega2*t));

% FEA stresses
stress = permute(stress,[1 3 2]);
sigma_xy = stress(3,:,1);


% Calculate the error
time_er = sqrt(sum((sigma_xy_exact - sigma_xy).^2)/sum(sigma_xy_exact.^2));

end