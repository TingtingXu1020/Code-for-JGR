clear
clc
close all
dbstop error
%% Input variables
% phase numbers, the last one is the matrix
Numb = 1; % Numb is the set of biotite inclusions (along the same direction)
phasen = Numb + 38; % In total 37 different crack orientations
psi = zeros(phasen, 1);
theta = zeros(phasen, 1);
phi = zeros(phasen, 1);

% Parameters Analysis
 psi_range = [0; pi/6; pi/4; pi/3; pi/2];
 theta_range = [0; pi/6; pi/4; pi/3; pi/2];
%psi_range = [pi/12; 5*pi/12];
% theta_range = [pi/12; 5*pi/12];

Rc_range = [0.5; 1; 2; 3];
rhoc_range = [0.01; 0.05; 0.5; 1; 5];
depth_range = [1;10;30;100];

% % Depth
for i = 3% : length(depth_range)
    psi(1) = psi_range(1);
    theta(1) = theta_range(1);
    R_c = Rc_range(4);
    rho_c = rhoc_range(3);
    depth = depth_range(i);
    
    [omega_save, Sigma, Eps, HC_save] = MainCode(psi, theta, phi, ...
        R_c, rho_c,depth);
end
% 
% % R_c
% for i = 1 : length(Rc_range)
%     psi(1) = psi_range(1);
%     theta(1) = theta_range(1);
%     R_c = Rc_range(i);
%     rho_c = rhoc_range(3);
%     depth = depth_range(3);
%     
%     [omega_save, Sigma, Eps, HC_save] = DamageGraniteMain(psi, theta, phi, ...
%         R_c, rho_c,depth);
% end
% 
% % rho_c
% for i = 1 : length(rhoc_range)
%     psi(1) = psi_range(1);
%     theta(1) = theta_range(1);
%     R_c = Rc_range(4);
%     rho_c = rhoc_range(i);
%     depth = depth_range(3);
%     
%     [omega_save, Sigma, Eps, HC_save] = DamageGraniteMain(psi, theta, phi, ...
%         R_c, rho_c,depth);
% end

% psi and theta
% for i = 1 : length(psi_range)
%     for j = 1 : length(theta_range)
%         psi(1) = psi_range(i);
%         theta(1) = theta_range(j);
%         R_c = Rc_range(4);
%         rho_c = rhoc_range(3);
%         depth = depth_range(3);
%         
%         [omega_save, Sigma, Eps, HC_save] = DamageGraniteMain(psi, theta, phi, ...
%             R_c, rho_c,depth);
%     end
% end
