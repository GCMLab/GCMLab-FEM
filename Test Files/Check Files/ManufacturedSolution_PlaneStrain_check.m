function [m_L2, m_e] = ManufacturedSolution_PlaneStrain_check(d1, d2, d3, s1, s2, s3, e1, e2, e3, Mesh1, Mesh2, Mesh3, Material, Control)
%MANUFACTUREDSOLUTION_PLANESTRAIN_CHECK calculates the convergence rates
%   m_L2 = MANUFACTUREDSOLUTION_PLANESTRAIN_CHECK(d1, d2, d3, s1, s2, s3, Mesh1,
%   Mesh2, Mesh3) calculates the rate of convergence of the L2 error norm
%   for using a manufactured solution
%
%   [m_L2, m_e] = MANUFACTUREDSOLUTION_PLANESTRAIN_CHECK(d1, d2, d3, s1, s2, s3, Mesh1,
%   Mesh2, Mesh3) also calculates the rate of convergence of the energy norm
%   for using a manufactured solution
%
%   ----------------------------------------------------------
%   Input
%   ----------------------------------------------------------
%   d1, d2, d3:             Displacement vectors from 3 runs with increasing mesh
%                           refinement
%   s1, s2, s3:             Nodal stress data from 3 runs with increasing mesh
%                           refinement
%   e1, e2, e3:             Nodal strain data from 3 runs with increasing mesh
%                           refinement
%   Mesh1, Mesh2, Mesh3:    Mesh data structure from 3 runs with increasing mesh
%                           refinement
%   1 - coarsest mesh, 2 - medium mesh, 3 - finest mesh

% Manufactured solution
%   ux = x^5 + x*y^3 - y^6
%   uy = x^5 + x*y^3 - y^6
%   exx = 5x^4 + y^3
%   eyy = 3xy^2 - 6y^5
%   exy = (3xy^2 - 6y^5 + 5x^4 + y^3) 
%   sxx = E/(1+v)/(1-2v)  * [(1-v)*exx + v*eyy ]
%   syy = E/(1+v)/(1-2v)  * [ v*exx    + (1-v)*eyy ]
%   sxy = E/(1+v)/(1-2v)  * (1-2v)/2 * exy

% Acknowledgements: Bruce Gee

plot_on = 0; % turn plots on/off - debugging tool

% Step 1 - Loop through each element and calculate the L2 and e-norms
eL2 = zeros(3,1);
eEN = zeros(3,1);
h   = zeros(3,1);

for sim = 1:3
    switch sim
        case 1
            d = d1;
            s = s1;
            e = e1;
            Mesh = Mesh1;
        case 2
            d = d2;
            s = s2;
            e = e2;
            Mesh = Mesh2;
        case 3
            d = d3;
            s = s3;
            e = e3;
            Mesh = Mesh3;
    end
    
    % Mesh size
    xI = Mesh.x(Mesh.conn(1,:),:);
    h(sim) = sqrt(polyarea(xI(:,1),xI(:,2)));
    
    % Calculate exact solutions
    x = Mesh.x(:,1);
    y = Mesh.x(:,2);
    d_exact = zeros(2*Mesh.nn,1);
    d_exact(1:2:end) = x.^5 + x.*y.^3 - y.^6;
    d_exact(2:2:end) = x.^5 + x.*y.^3 - y.^6;
    
    e_exact = zeros(3, Mesh.nn);
    s_exact = zeros(3, Mesh.nn);
    
    e_exact(1,:) = 5*x'.^4 + y'.^3;
    e_exact(2,:) = 3*x'.*y'.^2 - 6*y'.^5;
    e_exact(3,:) = ( 5*x'.^4 + y'.^3 + 3*x'.*y'.^2 - 6*y'.^5 );
    
    s_exact(1,:) = Material.Prop(1).E0/(1+Material.Prop(1).nu)/(1-2*Material.Prop(1).nu)* ((1-Material.Prop(1).nu)*e_exact(1,:) + Material.Prop(1).nu*e_exact(2,:));
    s_exact(2,:) = Material.Prop(1).E0/(1+Material.Prop(1).nu)/(1-2*Material.Prop(1).nu)* (Material.Prop(1).nu*e_exact(1,:) + (1-Material.Prop(1).nu)*e_exact(2,:));
    s_exact(3,:) = Material.Prop(1).E0/(1+Material.Prop(1).nu)/(1-2*Material.Prop(1).nu)* (1-2*Material.Prop(1).nu)/2 * e_exact(3,:);
    
    % Calculate error norms
    Quad = GlobalQuad(Mesh.nsd, Mesh.type, Control.qo);
    nq = Quad.nq;
    
    eL2_num = 0;
    eL2_den = 0;
    eEN_num = 0;
    eEN_den = 0;
    
    for i = 1:Mesh.ne
        enodes = Mesh.conn(i,:);
        xI = Mesh.x(enodes,:);
        
        for p = 1:nq
            N = Quad.Nq{p}';
            dNdxi = Quad.dNdxiq{p};
            wp = Quad.W(p);
            Je = dNdxi'*xI;
            
            % approximated displacement at quadrature point
            uxh_p = N*d(Mesh.xdofs(enodes));
            uyh_p = N*d(Mesh.ydofs(enodes));
            % exact displacement at quadrature point
            uxe_p = N*d_exact(Mesh.xdofs(enodes));
            uye_p = N*d_exact(Mesh.ydofs(enodes));
            % L2 norm
            eL2_num = eL2_num + [uxh_p - uxe_p , uyh_p - uye_p] * [uxh_p - uxe_p ; uyh_p - uye_p] * wp * det(Je);
            eL2_den = eL2_den + [uxe_p , uye_p] * [uxe_p ; uye_p] * wp * det(Je);
            
            % approximated stress and strain at quadrature point
            ehp = e(:,enodes)*N';
            shp = s(:,enodes)*N';
            % exact stress and strain
            eep = e_exact(:,enodes)*N';
            sep = s_exact(:,enodes)*N';
            % e norm
            eEN_num = eEN_num + (ehp-eep)'*(shp - sep) * wp * det(Je);
            eEN_den = eEN_den + eep'*sep * wp * det(Je);
            
            
            
        end
    end

  
    eL2(sim) = sqrt(eL2_num/eL2_den);
    eEN(sim) = sqrt(eEN_num/eEN_den);
    
end

% Determine slope of L2 norm
pL2 = polyfit(log(h), log(eL2),1);
m_L2 = pL2(1);

pEN = polyfit(log(h), log(eEN),1);
m_e = pEN(1);


% Step 2 - Calculate the slope of each curve
if plot_on
figure(1)
loglog(h,eL2,'o')
hold on
loglog(h,exp(pL2(2))*h.^pL2(1))
hold off
xlabel('Mesh size')
ylabel('L2-norm')

figure(2)
loglog(h,eEN,'o')
hold on
loglog(h,exp(pEN(2))*h.^pEN(1))
hold off
xlabel('Mesh size')
ylabel('e-norm')
end




end

