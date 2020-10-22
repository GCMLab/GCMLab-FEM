function [W,Q] = quadrature( quadorder, type, dim )
%QUADRATURE Returns a n x 1 column vector W of quadrature
%weights and a n x dim matrix of quadrature points, where n is the
%number of quadrature points.  The function is called as follows:
%   [W,Q] = quadrature( quadorder, type, dim )
%
%   ----------------------------------------------------------------------
% 	Input
%   ----------------------------------------------------------------------
%   quadorder: quadrature order. If an invalid quadrature order for
%               triangular quadrature is given, the default is 1.
%   type: type of quadrature (i.e. gaussian, triangular, etc.. );
%           The default for type is GAUSS.
%   dim: number of spatial dimensions of the problem. The default for dim 
%           is unity.
%
%   ----------------------------------------------------------------------
%   References
%   ----------------------------------------------------------------------
%   [1] Fish, J., & Belytschko, T. (2007). A first course in finite 
%       elements. Chichester; Hoboken, NJ: John Wiley & Sons.
%
%   ----------------------------------------------------------------------
%   Acknowledgements: Jack Chessa

% -------------------------------------------------------------------
% set default type to Gaussian if not specified
% -------------------------------------------------------------------
if ( nargin < 2 )
    type = 'GAUSS';
end

% -------------------------------------------------------------------
% set default dimensions to 1 if type is Gaussian, 2 otherwise
% -------------------------------------------------------------------
if nargin < 3
    if strcmp(type,'GAUSS') == 1
        dim = 1;
    else
        dim = 2;
    end
end

% -------------------------------------------------------------------
% Quadrature points for Gaussian quadrature
% -------------------------------------------------------------------
if ( strcmp(type,'GAUSS') == 1 )

    quadpoint = zeros(quadorder^dim, dim);
    quadweight = zeros(quadorder^dim, 1);

    r1pt = zeros(quadorder,1);
    r1wt = zeros(quadorder,1);

    switch ( quadorder ) 
    case 1  % Ref. [1] pg. 89
        r1pt(1) = 0.000000000000000;
        r1wt(1) = 2.000000000000000;

    case 2  % Ref. [1] pg. 89
        r1pt(1) =  0.577350269189626;   % 1/sqrt(3)                   
        r1pt(2) = -0.577350269189626;   % -1/sqrt(3)

        r1wt(1) = 1.0; 
        r1wt(2) = 1.0;  

    case 3  % Ref. [1] pg. 89
        r1pt(1) =  0.774596669241483;
        r1pt(2) = -0.774596669241483;
        r1pt(3) =  0.000000000000000;

        r1wt(1) = 0.555555555555556;
        r1wt(2) = 0.555555555555556; 
        r1wt(3) = 0.888888888888889;   

    case 4  % Ref. [1] pg. 89
        r1pt(1) =  0.861134311594053;
        r1pt(2) = -0.861134311594053;
        r1pt(3) =  0.339981043584856;
        r1pt(4) = -0.339981043584856;

        r1wt(1) = 0.347854845137454;
        r1wt(2) = 0.347854845137454; 
        r1wt(3) = 0.652145154862546;
        r1wt(4) = 0.652145154862546;  

    case 5  % Ref. [1] pg. 89
        r1pt(1) =  0.906179845938664;
        r1pt(2) = -0.906179845938664;
        r1pt(3) =  0.538469310105683;
        r1pt(4) = -0.538469310105683;
        r1pt(5) =  0.000000000000000;

        r1wt(1) = 0.236926885056189;
        r1wt(2) = 0.236926885056189;
        r1wt(3) = 0.478628670499366;
        r1wt(4) = 0.478628670499366;  
        r1wt(5) = 0.568888888888889;  

    case 6  % Ref. [1] pg. 89
        r1pt(1) =  0.932469514203152;
        r1pt(2) = -0.932469514203152;
        r1pt(3) =  0.661209386466265;
        r1pt(4) = -0.661209386466265;
        r1pt(5) =  0.238619186003152;
        r1pt(6) = -0.238619186003152;

        r1wt(1) = 0.171324492379170;
        r1wt(2) = 0.171324492379170;
        r1wt(3) = 0.360761573048139;
        r1wt(4) = 0.360761573048139;   
        r1wt(5) = 0.467913934572691; 
        r1wt(6) = 0.467913934572691;

    case 7
        r1pt(1) =  0.949107912342759;
        r1pt(2) = -0.949107912342759;
        r1pt(3) =  0.741531185599394;
        r1pt(4) = -0.741531185599394;
        r1pt(5) =  0.405845151377397;
        r1pt(6) = -0.405845151377397;
        r1pt(7) =  0.000000000000000;

        r1wt(1) = 0.129484966168870;
        r1wt(2) = 0.129484966168870;
        r1wt(3) = 0.279705391489277;
        r1wt(4) = 0.279705391489277;
        r1wt(5) = 0.381830050505119;
        r1wt(6) = 0.381830050505119;
        r1wt(7) = 0.417959183673469;

    case 8
        r1pt(1) =  0.960289856497536;
        r1pt(2) = -0.960289856497536;
        r1pt(3) =  0.796666477413627;
        r1pt(4) = -0.796666477413627;
        r1pt(5) =  0.525532409916329;
        r1pt(6) = -0.525532409916329;
        r1pt(7) =  0.183434642495650;
        r1pt(8) = -0.183434642495650;

        r1wt(1) = 0.101228536290376;
        r1wt(2) = 0.101228536290376;
        r1wt(3) = 0.222381034453374;
        r1wt(4) = 0.222381034453374;
        r1wt(5) = 0.313706645877887;
        r1wt(6) = 0.313706645877887;
        r1wt(7) = 0.362683783378362;
        r1wt(8) = 0.362683783378362;

    otherwise
        r1pt = roots(legpol(quadorder));
        r1pt = sort(r1pt);
        P  = legendre(quadorder,r1pt);
        for i = 1:quadorder
            r1wt(i) = 2/(P(2,i)^2);
        end

    end  % end of quadorder switch

    n = 1;  % counter

    if ( dim == 1 ) 
        quadpoint = r1pt;           
        quadweight = r1wt; 

    elseif ( dim == 2 )     % Ref. [1] pg. 180
        for i = 1:quadorder
            for j = 1:quadorder
                quadpoint(n,:) = [r1pt(i), r1pt(j)];           
                quadweight(n) = r1wt(i)*r1wt(j); 
                n = n+1;
            end
        end

    else % dim == 3
        for i = 1:quadorder
            for j = 1:quadorder
                for k = 1:quadorder
                    quadpoint(n,:) = [r1pt(i), r1pt(j), r1pt(k)];           
                    quadweight(n) = r1wt(i)*r1wt(j)*r1wt(k); 
                    n = n+1;
                end
            end
        end

    end

    Q = quadpoint;
    W = quadweight;
    % END OF GAUSSIAN QUADRATURE DEFINITION

% -------------------------------------------------------------------
% Triangular quadrature
% -------------------------------------------------------------------
elseif ( strcmp(type,'TRIANGULAR') == 1 ) 

    if ( dim == 3 )  %%% TETRAHEDRA

        if (quadorder ~= 1 &&  quadorder ~= 2 &&  quadorder ~= 3 && ...
                quadorder ~= 5) 
            % check for valid quadrature order
            disp(['Incorect quadrature order for triangular quadrature.'...
                    'Quadrature order of 1 used.']);
            quadorder = 1;
        end

        if  ( quadorder == 1 )  % Ref. [1] pg. 185
            quadpoint = [ 0.25 0.25 0.25 ];
            quadweight = 1;

        elseif ( quadorder == 2 ) 
            quadpoint = [   0.58541020  0.13819660  0.13819660;
                            0.13819660  0.58541020  0.13819660;
                            0.13819660  0.13819660  0.58541020;
                            0.13819660  0.13819660  0.13819660  ];
            quadweight = [1; 1; 1; 1]/4;

        elseif ( quadorder == 3 ) 
            quadpoint = [   0.25  0.25  0.25;
                            1/2   1/6   1/6;
                            1/6   1/2   1/6;
                            1/6   1/6   1/2;
                            1/6   1/6   1/6 ];
            quadweight = [-4/5 9/20 9/20 9/20 9/20]';

        elseif (quadorder == 5 )
            a = 0.25;
            b1 = (7-sqrt(15))/34;
            b2 = (7+sqrt(15))/34;
            c1 = (13+3*sqrt(15))/34;
            c2 = (13-3*sqrt(15))/34;
            wa = 112/5670;
            w1 = (2665+14*sqrt(15))/226800;
            w2 = (2665-14*sqrt(15))/226800;
            d = (5 - sqrt(15))/20;
            e = (5 + sqrt(15))/20;
            wde = 5/567;

            quadpoint = [ a  a  a ;
                          b1 b1 b1;
                          b1 b1 c1;
                          b1 c1 b1;
                          c1 b1 b1;
                          b2 b2 b2;
                          b2 b2 c2;
                          b2 c2 b2;
                          c2 b2 b2;
                          d  d  e ;
                          d  e  d ;
                          e  d  d ;
                          d  e  e ;
                          e  d  e ;
                          e  e  d ];
            quadweight = 6*[wa w1 w1 w1 w1 w2 w2 w2 w2 wde ...
                            wde wde wde wde wde]';          

        end

        Q = quadpoint;
        W = quadweight/6;

% -------------------------------------------------------------------
% TRIANGLES
% -------------------------------------------------------------------
    else  

        if ( quadorder == 1 )   % set quad points and quadweights
            quadpoint = [ 0.3333333333333, 0.3333333333333 ];
            quadweight = 1;

        elseif ( quadorder == 2 ) 
            quadpoint = zeros( 3, 2 );
            quadweight = zeros( 3, 1 );

            quadpoint(1,:) = [ 0.1666666666667, 0.1666666666667 ];
            quadpoint(2,:) = [ 0.6666666666667, 0.1666666666667 ];
            quadpoint(3,:) = [ 0.1666666666667, 0.6666666666667 ]; 

            quadweight(1) = 0.3333333333333; 
            quadweight(2) = 0.3333333333333; 
            quadweight(3) = 0.3333333333333;   

        elseif (quadorder <= 5 )
            % TODO: Find references for this
            a = (6+sqrt(15))/21;
            b = 4/7-a;
            A = 2*( (155+sqrt(15))/2400 );
            Aby2 = A/2;
            B = 2* 31/240-A  ;
            Bby2 = B/2;
            quadpoint(1,:) = [1/3,1/3];
            quadpoint(2,:) = [a,a];
            quadpoint(3,:) = [1-2*a,a];
            quadpoint(4,:) = [a,1-2*a];
            quadpoint(5,:) = [b,b];
            quadpoint(6,:) = [1-2*b,b];
            quadpoint(7,:) = [b,1-2*b];

            quadweight(1) = 18/80;
            quadweight(2) = A;
            quadweight(3) = A;
            quadweight(4) = A;
            quadweight(5) = B;
            quadweight(6) = B;
            quadweight(7) = B;

        elseif ( quadorder == 6 )
            a = 0.063089014491502;
            b = 0.249286745170910;
            c = 0.310352451033785;
            d = 0.053145049844816;

            quadpoint(1,:) = [a,a];
            quadpoint(2,:) = [1-2*a,a];
            quadpoint(3,:) = [a,1-2*a];
            quadpoint(4,:) = [b,b];
            quadpoint(5,:) = [1-2*b,b];
            quadpoint(6,:) = [b,1-2*b];
            quadpoint(7,:) = [c,d];
            quadpoint(8,:) = [d,c];
            quadpoint(9,:) = [1-(c+d),c];
            quadpoint(10,:) = [1-(c+d),d];
            quadpoint(11,:) = [c,1-(c+d)];
            quadpoint(12,:) = [d,1-(c+d)];

            quadweight(1) = 0.025422453185103;
            quadweight(2) = 0.025422453185103;
            quadweight(3) = 0.025422453185103;
            quadweight(4) = 0.058393137863189;
            quadweight(5) = 0.058393137863189;
            quadweight(6) = 0.058393137863189;
            quadweight(7) = 0.041425537809187;
            quadweight(8) = 0.041425537809187;
            quadweight(9) = 0.041425537809187;
            quadweight(10) = 0.041425537809187;
            quadweight(11) = 0.041425537809187;
            quadweight(12) = 0.041425537809187;

            quadweight = quadweight*2;

        elseif ( quadorder == 8 ) 
            quadpoint = zeros( 13, 2 );
            quadweight = zeros( 13, 1 );

            quadpoint(1 ,:) = [ 0.0651301029022, 0.0651301029022 ];
            quadpoint(2 ,:) = [ 0.8697397941956, 0.0651301029022 ];
            quadpoint(3 ,:) = [ 0.0651301029022, 0.8697397941956 ];
            quadpoint(4 ,:) = [ 0.3128654960049, 0.0486903154253 ];
            quadpoint(5 ,:) = [ 0.6384441885698, 0.3128654960049 ];
            quadpoint(6 ,:) = [ 0.0486903154253, 0.6384441885698 ];
            quadpoint(7 ,:) = [ 0.6384441885698, 0.0486903154253 ];
            quadpoint(8 ,:) = [ 0.3128654960049, 0.6384441885698 ];
            quadpoint(9 ,:) = [ 0.0486903154253, 0.3128654960049 ];
            quadpoint(10,:) = [ 0.2603459660790, 0.2603459660790 ];
            quadpoint(11,:) = [ 0.4793080678419, 0.2603459660790 ];
            quadpoint(12,:) = [ 0.2603459660790, 0.4793080678419 ];
            quadpoint(13,:) = [ 0.3333333333333, 0.3333333333333 ];

            quadweight(1 ) =  0.0533472356088;
            quadweight(2 ) =  0.0533472356088; 
            quadweight(3 ) =  0.0533472356088;
            quadweight(4 ) =  0.0771137608903;
            quadweight(5 ) =  0.0771137608903;
            quadweight(6 ) =  0.0771137608903;
            quadweight(7 ) =  0.0771137608903;
            quadweight(8 ) =  0.0771137608903;
            quadweight(9 ) =  0.0771137608903;
            quadweight(10) =  0.1756152576332; 
            quadweight(11) =  0.1756152576332; 
            quadweight(12) =  0.1756152576332;
            quadweight(13) = -0.1495700444677; 

        elseif ( quadorder > 8 )
            [quadweight,quadpoint] = GaussQuadTri(quadorder);
            quadweight = quadweight*2;
            
        end
        
        Q = quadpoint;
        W = quadweight/2;   % TODO: Figure out why divided quadweight by 2
        
    end

end  % end of TRIANGULAR initialization

% END OF FUNCTION

%% 
function poly = legpol(n)
% Returns the coefficients of the Legendre polynomial of order n

    switch n
    case 0,
        poly = 1;
    case 1,
        poly = [1 0];
    otherwise,
        poly2 = 1;
        poly1 = [1 0];
        for m = 1:n-1
            poly = ((2*m+1)/(m+1))*[poly1 0] - m/(m+1)*[0 0 poly2];
            poly2 = poly1;
            poly1 = poly;
        end
    end
end

end