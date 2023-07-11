function [Fint,d_m] = getFint_ini(BC, d_m, K, C, M, dt, alpha, TimeType)
%GETFint_ini Internal forces at initial conditions
%   Fint = GETFint_ini(BC, d_m, K, C, M, dt, alpha, TimeType) is a column vector of internal forces 
%   acting on each degree of freedom at the initial conditions
%   
%   --------------------------------------------------------------------
%   Input
%   --------------------------------------------------------------------
%   BC:     Structure array with the following fields,
%           .b                          Anonymous function of distributed
%                                       body force (size 1 x nsd)
%           .traction_force_node        Column vector of nodes with 
%                                       prescribed tractions
%           .traction_force_value       Column vector of prescribed tractions
%                                       on nodes
%           .traction_force_dof         Column vector of degrees of freedom
%                                       with prescribed tractions
%           .traction_force_dof_value   Column vector of prescribed tractions
%                                       on DOF
% 
%   d_m:    Vector with the displacements from previous timesteps
%   K:      Stiffness matrix
%   C:      Damping matrix
%   M:      Mass matrix
%   dt:     Timestep size at initial conditions
%   alpha:  Time integration parameter
%   TimeType:   0 - static, 1 - transient, 2 - dynamic

        switch TimeType
            case 0 % Static
                Fint = K*d_m.d;
            case 1 % Transient (1st order time derivative)
                Fint = (alpha*Klin+(1/dt)*C)*d_m.d + ((1-alpha)*Klin - (1/dt)*C )*d_m.dnm1;
            case 2 % Dynamic
                d_m.d = d0;  % d at timestep n-1
                d_m.d(BC.fix_disp_dof) = BC.fix_disp_value(t-dt);
                d_m.dnm1 = d0;  % d at timestep n-2
                d_m.dnm1(BC.fix_disp_dof) = BC.fix_disp_value(t-2*dt);
                d_m.dnm2 = d0;  % d at timestep n-3
                d_m.dnm2(BC.fix_disp_dof) = BC.fix_disp_value(t-3*dt);
                d_m.dnm3 = d0;  % d at timestep n-4
                d_m.dnm3(BC.fix_disp_dof) = BC.fix_disp_value(t-4*dt);
                d_m.dnm4 = d0;  % d at timestep n-5
                d_m.dnm4(BC.fix_disp_dof) = BC.fix_disp_value(t-5*dt);
                                                
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
                a = (d - d_temp)./dt.^2./bet; %Acceleration at n-1
                Fint = M*a +(1+alpha)*C*(v_temp-gam*d_temp/(dt*bet)+gam*d/(dt*bet))...
                    +(1+alpha)*Klin*d - alpha*(C*vnm1 + Klin*dnm1);
                
                % Update vectors or structures from previous timesteps
                d_m.dnm4 = d_m.dnm3;                       % d vector from timestep n-4
                d_m.dnm3 = d_m.dnm2;                       % d vector from timestep n-3
                d_m.dnm2 = d_m.dnm1;                       % d vector from timestep n-2
                d_m.dnm1 = d_m.d;                          % d vector from timestep n-1
        end
end
