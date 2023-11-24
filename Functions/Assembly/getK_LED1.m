function [K, R, Fint] = getK_LED1(~, ~, ~, ~, Fext, ~, Klin, M, d_m, dt, ~, C, alpha)
%GETK_LED1 Relevant matrices for HHT dynamic solver
%   K = GETK_LED1(...) returns the coefficient matrix for the 
%   iterative solver where the problem uses a linear elastic material
%   
%   [K, R] = GETK_LED1(...) also returns the residual vector R for the 
%   iterative solver where the problem uses a linear elastic material
%
%   [K, R, Fint] = GETK_LED1(...) also returns the internal force vector for the 
%   iterative solver where the problem uses a linear elastic material
%
%   --------------------------------------------------------------------
%   Accepted Inputs (in order)
%   --------------------------------------------------------------------
%   getK_LED1(~, ~, ~, ~, Fext, ~, Klin, M, d_m, dt, ~ ,C,alpha)
%   Mesh:       ~
%  
%   Quad:       ~
% 
%   Material:   ~
%   
%   Fintnm1:    ~
%   Fext:       External force vector at timestep n
%   Fextnm1:    ~
%   Klin:       Linear elastic stiffness matrix
%   M:          Mass matrix
%   d_m:        Structure array with the following fields
%               d:          unconverged degree of freedom vector at current timestep n and iteration
%               dnm1:       converged degree of freedom vector at timestep n-1
%               dnm2:       converged degree of freedom vector at timestep n-2
%               dnm3:       converged degree of freedom vector at timestep n-3
%   dt:         timestep size between timesteps n-1 and n
%   dtnm1:      ~
%   C:          Linear damping coefficient matrix
%   alpha:      intagration parameter


% Compute constants
gam = 1/2-alpha;
bet = (1-alpha)^2/4;

d = d_m.d;
dnm1 = d_m.dnm1;
dnm2 = d_m.dnm2;
dnm3 = d_m.dnm3;
dnm4 = d_m.dnm4;

% 2nd order accurate backwards difference approximation
vnm1 = 1/2/dt* (3*dnm1 - 4*dnm2 + dnm3);
anm1 = 1/dt^2 * (2*dnm1  - 5*dnm2 + 4*dnm3 - dnm4);

d_temp = dnm1+ dt*vnm1 + dt^2/2*(1-2*bet)*anm1;
v_temp = vnm1 + dt*(1-gam)*anm1;

% Internal forces
a = (d - d_temp)./dt^2./bet;
Fint = M*a +(1+alpha).*C*(v_temp-gam*d_temp./(dt*bet)+gam*d./(dt*bet))...
    +(1+alpha).*Klin*d - alpha.*(C*vnm1 + Klin*dnm1);

% Jacobian
K = M./dt^2./bet + (1+alpha).*gam./dt./bet.*C + (1+alpha).*Klin;

% Residual
% R = K*d - (Fext + alpha.*(C*vnm1 + K*dnm1) + M*d_temp./dt^2./bet -...
%     (1+alpha).*C*(v_temp-gam.*d_temp./dt./bet));

% Residual in terms of dnm1, vnm1, and anm1 (from maple)
R = - K*d + (Fext - ...
    (((-alpha*K*dt^2*bet - C*gam*(1 + alpha)*dt - M)/(dt^2*bet))*dnm1 + ...
    ((C*((-alpha - 1)*gam + bet)*dt - M)/(dt*bet))*vnm1 + ...
    (((2*C*(1 + alpha)*dt + 2*M)*bet - C*gam*(1 + alpha)*dt - M)/(2*bet))*anm1));

end
