clear; close all; clc;

% Loading conditions
P = 4.5;            % point load [kN]
L = 1.41;           % length between supports [m]
Ltot = 1.5;         % total beam length [m]

% Aluminium constants
rho = 2710;         % density [kg/m3]
E = 71;             % Young's modulus [kN/mm2]
nu = 0.33;          % Poisson's ratio
fy = 90;            % yield stress [N/mm2]
tauy = fy/sqrt(3);  % yield shear stress [N/mm2]
fu = 130;           % ultimate stress [N/mm2]
eps_u = 0.12;       % ultimate strain [%]

% Section dimensions [mm]
b = 110; h = 160;     % width & height
tf = 1; tw = 1;       % flange & web thickness
ta = 1.6;               % angle connection thickness
c1 = 19; c2 = 12.7;     % angle size

h_est = 0; iter = 0;
while abs(h - h_est) > 0.01
    %% Centroid height, h_cen [mm]
    % Total cross-sectional area, A
    area = 2*(b*tf + tw*(h-2*tf)... 
        + (c1^2 - (c1-ta)^2) + (c2^2 - (c2-ta)^2));
    
    % Total moment about bottom edge, M
    M = tf*b*(h-0.5*tf) + tf*b*(0.5*tf)...
        + 2*tw*(h-2*tf)*(tf + 0.5*(h-2*tf))...
        + 2*ta*(c1-ta)*(h-tf-0.5*ta) + 2*ta*c1*(h-tf-0.5*c1)...
        + 2*ta*(c2-ta)*(tf+0.5*ta) + 2*ta*c2*(tf+0.5*c2);
    
    % Centroid height from bottom edge
    hcen = M/area + tf;
    
    %% Second moment of area of section, Iyy [mm4]
    
    % Flange
    Iyy_flange = 2*b*tf^3/12 + b*tf*((hcen-0.5*tf)^2 + (h-hcen-0.5*tf)^2);
    
    % Web
    Iyy_web = 2*((h-2*tf)*tw^3/12 + tw*(h-2*tf)*(0.5*h-hcen)^2);
    
    % Angles (using centreline approximations)
    Iyy_angles = 2*(ta*c1^3/12 + ta*c1*(h-hcen-tf-0.5*c1)^2)...
               + 2*((c1-ta)*ta^3/12 + ta*(c1-ta)*(h-hcen-tf-0.5*ta)^2)...
               + 2*(ta*c2^3/12 + ta*c2*(hcen-tf-0.5*c2)^2)...
               + 2*((c2-ta)*ta^3/12 + ta*(c2-ta)*(hcen-tf-0.5*ta)^2);
    
    % Iyy total
    Iyy = Iyy_flange + Iyy_web + Iyy_angles;
    
    %% Check yield stress
    My = P*L/4;
    h_est = fy*Iyy/(My*1e6);
    
    step_size = 0.5;
    h = (1+step_size)*h - step_size*h_est;
    iter = iter + 1;
end
fprintf('In %i iterations, ideal h = %.2f mm \n', iter, h)


%% Torsion constant, J [mm4]
% using centerline approximations & neglecting angles
J = 4 * ((b-tw)*(h-tf))^2 / ((2*(h-tw)/tw) + (2*(b-tf)/tf)); 


%% Checks
% ------ MATERIAL YIELDING (BENDING)------
sigmaxx = My*hcen/Iyy;
if sigmaxx < fy
    disp('Material yielding (bending) check OK.');
else
    disp('Material yielding (bending) check FAILED.');
end

% ------ MATERIAL YIELDING (SHEAR)------
Vz = P/2 *10^3;   % Maximum shear force (midpoint) [N]

% Shear flow of doubly-symmetric box section
q_a = tf*(b/2-c1) * (h+tf-hcen+tf/2) * Vz/Iyy;
q_b = q_a + (tf*(c1-ta)*(h+tf-hcen+tf/2) + (c1-ta)*ta*(h+tf-hcen-ta/2)) * Vz/Iyy;
q_c = q_b + ((c1-ta)*tw + ta*(c1-ta)) * (h+tf-hcen-ta-(c1-ta)/2) * Vz/Iyy;

% Finding maximum shear flow (at the web)
x = [0:0.1:(h-c1-c2)];
q_web = q_c + x*tw .* (h+tf-c1-x/2-hcen) .* Vz/Iyy;
[q_d, I] = max(q_web);
z_maxq = h+tf-c1-x(I); % Location where max shear flow occurs (from bottom)

q_e = q_web(end);
q_f = q_e + ((c2-ta)*tw + ta*(c2-ta)) * -(hcen-tf-ta-(c2-ta)/2) * Vz/Iyy;
q_g = q_f + (tf*(c2-ta)*-(hcen-tf/2) + (c2-ta)*ta*-(hcen-tf-ta/2)) * Vz/Iyy;


% Check for maximum shear stress, taux_max [N/mm^2]
taux_topf = max([q_a, q_b])/tf;               % Max shear stress at top flange
taux_web = max([q_b, q_c, q_d, q_e, q_f])/tw; % Max shear stress at web
taux_botf = max([q_f, q_g])/tf;               % Max shear stress at bot flange

taux_max = max([taux_topf, taux_web, taux_botf]);

% Check if beam will yield in shear
if taux_max < tauy
    disp('Material yielding (shear) check OK.');
else
    disp('Material yielding (shear) check FAILED.');
end

% ------ PLASTIC HINGE FORMATION ------
% Determine plastic neutral axis location from bottom, hplas
hplas = -(ta*(-2*c1+2*c2)-h*tw-2*tf*tw)/(2*tw);
% Determine plastic section modulus, Zp
Zp = b*tf*(h+tf+tf/2-hplas) +  2*((h+tf-hplas)*tw*((h+tf-hplas)/2)) + ...
    2*(c1*ta*(h+tf-c1/2-hplas) + (c1-ta)*ta*(h+tf-ta/2-hplas)) + ...
    b*tf*(hplas-tf/2) + 2*((hplas-tf)*tw*((hplas-tf)/2)) + ...
    2*(c2*ta*(hplas-tf-c2/2) + (c2-ta)*ta*(hplas-tf-ta/2));
% Determine plastic moment, Mp [kNm]
Mp = fy*Zp/10^6;

% Check for formation of plastic hinge
if My < Mp
    disp('Plastic hinge formation check OK.')
else
    disp('Plastic hinge formation check FAILED.')
end




%% Cost
mass = (area*1e-6)*Ltot*rho;   % mass [kg]

