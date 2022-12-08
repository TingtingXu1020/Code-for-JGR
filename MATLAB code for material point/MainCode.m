function [omega_save, Sigma, Eps, HC_save] = MainCode(psi, theta, phi, ...
    R_c, rho_c,depth)
%dbstop warning
%dbstop if naninf
%% Input variables
% phase numbers, the last one is the matrix
Numb = 1; % Numb is the set of biotite inclusions (along the same direction)
phasen = Numb + 38; % In total 37 different crack orientations
% Interval in each time step  (s)
inter = 3600*24*30*12;
% Weathering time (year)
wtime = (3600*24*30*12/inter)*1e3;
%T = 0; % Time (s)
BC = 1; % Proportional Loading
pro_bio = 0.15; % Abundance of biotite

%% Material properties
% Bulk modulus and shear modulus for Boitite and Vermiculite
bbm = 76700;
bsm = 41600;
Everm = 14100;
Enu = 0.33;
vbm = Everm / 3 / (1-2*Enu);
vsm = Everm / 2 / (1 + Enu);

% % Poisson's ratio
% POISONB = (3*bbm-2*bsm) / 2 / (3*bbm+bsm);
% POISONV = (3*vbm-2*vsm) / 2 / (3*vbm+vsm);

% Matrix Properties
sm = zeros(phasen, 1);
bm = zeros(phasen, 1);
msm = 3.13e4;
mbm = 6.07e4;
sm(phasen) = msm;
bm(phasen) = mbm;
%POISONM = (3*mbm-2*msm)/2/(3*mbm+msm);
HC = zeros(3, 3, 3, 3);

%% Crack parameters
% crack orientation
crackorien = rotation37(phasen - (Numb + 1));
crackdirection = orientation37(phasen - (Numb + 1));
orientensor = zeros(3, 3, 37);
for i = 1 : 37
    orientensor( :, :, i) = crackdirection(:, i) * crackdirection(:,i)';
end

% crack weights
crackweight = weight37(37);

% Euler's angle
% For crack inclusion. there is no spinal rotation (phi=0)
psi((Numb + 1) : phasen - 1) = crackorien(1,:);
theta((Numb + 1) : phasen - 1) = crackorien(2,:);
phi((Numb + 1) : phasen - 1) = crackorien(3,:);
Q = zeros(3, 3, phasen);
for i = 1 : phasen
    Q(:, :, i) = transmatrixo(psi(i), theta(i), phi(i));
end

% Identity tensor
[~, fourthI, ~, fourthJ] = Identity(); % fourthK = fourthI - fourthJ;

%% Boundary conditions
% External  Stress (Proportional Loading)
ks = 0.5; % lateral stress coefficient
sig_v = 2e3 * 9.81e-3 * depth * 1e-3; % MPa
tstress_aim = [-sig_v * ks 0 0; 0 -sig_v * ks 0; 0 0 -sig_v ];
tstress = zeros(3, 3);
tstrain = zeros(3, 3);

% Volume fraction and geometry
Vrev = 1;
pr = zeros(phasen, 1);
ao = ones(phasen, 1);
co = ones(phasen, 1);
X_b = 1/3; % c/a
for i = 1 : Numb
    ao(i) = 1e-3; % m
    N_bio = Vrev * pro_bio/ao(i) ^3 *(3/4/pi/X_b);
    co(i) = X_b * ao(i);
    pr(i) = pro_bio/Numb;
end

N_crack = 1e5; % Crack weight for each direction
X_c = 1e-8; % the aspect ratio of crack should be very small
ao((Numb + 1) : phasen - 1) = 1e-4;
co((Numb + 1) : phasen - 1) = ao((Numb + 1) : phasen - 1) * X_c;
ao(39) = 1 + 1e-3;
co(39) = 1 + 1e-3;

% a,c are initialized
a = ao;
c = co;
rho = N_crack .* a(1 + Numb : end - 1).^3 /Vrev;
%rho_save(1:37, 1 : wtime) = 0;
pro_cracks = 4/3*pi*X_c*rho/Vrev; % volume fraction of cracks are small
pr(Numb+1 : phasen-1) = pro_cracks;
pr(phasen) = 1 - sum(pr(1:phasen-1));
drho = zeros(37, 1);

%% Save matrices
% Stiffness Tensor
pC = zeros(3, 3, 3, 3, phasen);

% Strain
pstrain = zeros(3, 3, phasen);
Lstrain = zeros(3, 3, phasen); % Local strain
wr = zeros(Numb,1);
Letastrain = zeros(3, 3, Numb); % Local eigenstrain
Getastrain = zeros(3, 3, Numb); %  Global eigenstrain
%leps = zeros(3, 3, 37); % Local strain of crack inclusion
neps = zeros(37, 1);  % Normal strain on the crack surface

% Influence tensor
Indi = zeros(37, 1);

% Damage matrix
omega = zeros(3, 3);
for num = 1 : 37
    omega = omega + crackweight(num)* rho(num)*orientensor(:, :, num);
end

% Loading matrix
tinc = 1000;
omega_save_L = zeros(3, 3, tinc);
Sigma_L = zeros(3, 3, tinc);
Eps_L = zeros(3, 3, tinc);
HC_save_L(1:3, 1:3, 1:3, 1:3, 1:tinc) = 0;

% Time evolution matrix
HC_save(1:3, 1:3, 1:3, 1:3, 1:wtime) = 0;
omega_save = zeros(3, 3, wtime);
Sigma = zeros(3, 3, wtime);
Eps = zeros(3, 3, wtime);
First = zeros(37, wtime);
Second = zeros(37, wtime);
Third = zeros(37, wtime);

%% Mechanical loading

for inc = 1 : tinc
    dstress = tstress_aim/tinc;
    omega_save_L(:, :, inc) = omega;
    Sigma_L(:, :, inc) = tstress;
    Eps_L(:, :, inc) = tstrain;
    HC_save_L(:, :, :, :, inc) = HC;
    %% Weathering
    % Weathering rate
    for i = 1 : Numb
        wr(i) = 0;
        Letastrain(3, 3, i) = 0.4*wr(i)/1.4;
        Getastrain(:, :, i) = transpose(Q(:, :, i))*Letastrain(:, :, i)*Q(:, :, i) ;
        sm(i) = (1 - wr(i)) * bsm + vsm*wr(i);
        bm(i) = (1 - wr(i)) * bbm + vbm*wr(i);
        pr(i) = 4/3*pi * N_bio/Numb * a(i)^2 *c(i);
    end
    
    pro_cracks = 4/3*pi*X_c*rho/Vrev; % volume fraction of cracks are small
    pr(Numb+1 : phasen-1) = pro_cracks;
    pr(phasen) = 1 - sum(pr(1:phasen-1));
    
    % Stiffness tensors for Biotite and Matrix
    for i = 1 : Numb
        pC(:, :, :, :, i) = matrixc(bm(i), sm(i));
    end
    pC(:, :, :, :, phasen)  = matrixc(bm(phasen), sm(phasen));
    
    % Stiffness tensors for cracks
    for i = (Numb + 1) : phasen - 1
        pC(:, :, :, :, i) =  Indi(i - Numb)*3*mbm*fourthJ;
    end
    
    %% First step of homogenization -- matrix and cracks
    % Hill's tensors
    Closeset = []; %find(Indi > 0) + Numb;
    Openset = [setdiff((2:38)', Closeset); phasen];
    
    pP = zeros(3, 3, 3, 3, phasen);
    for j = 1 : size(Openset, 1)
        i = Openset(j);
        TEMP = Hill(a(i), c(i), bm(phasen), sm(phasen)); % Since the matrix is isotropic, there is an explicit formulation
        pP(:, :, :, :, i)  = Transform(TEMP, Q(:, :, i)'); % Local to Global
    end
    
    % Concentration tensors
    pAo = zeros(3, 3, 3, 3, phasen);
    for j = 1 : size(Openset, 1)
        i = Openset(j);
        pAo(:, :, :, :, i) = inversegeneral(fourthI + doubledotff(pP(:, :, :, :, i), (pC(:, :, :, :, i) ...
            -pC(:, :, :, :, phasen)))); % The matrix stiffness here is the reference stiffness tensor
    end
    pAo(:, :, :, :, phasen) = fourthI; % For matrix itself
    
    % Average of infinite conentration tensor
    sAo= zeros(3, 3, 3, 3);
    for j = 1 : size(Openset, 1)
        i = Openset(j);
        sAo = sAo + pr(i)*pAo(:, :, :, :, i);
    end
    
    pA = zeros(3, 3, 3, 3, phasen);
    for j = 1 : size(Openset, 1)
        i = Openset(j);
        pA(:, :, :, :, i) = doubledotff(pAo(:, :, :, :, i), inversegeneral(sAo));
    end
    
    HC_Do = zeros(3, 3, 3, 3); % Damaged homogenized tensor with open crack
    for j = 1 : size(Openset, 1)
        i = Openset(j);
        HC_Do = HC_Do + pr(i)*doubledotff(pC(:, :, :, :, i), pA(:, :, :, :, i));
    end
    
    %% closed cracks
    for j = 1 : size(Closeset, 1)
        i = Closeset(j);
        TEMP = Hill_P([a(i), a(i), c(i)], 16, 16, HC_Do);
        pP(:, :, :, :, i)  = Transform(TEMP, Q(:, :, i)'); % Local to Global
    end
    %TEMP = Hill_P([a(phasen), a(phasen), c(phasen)], 16, 16, HC_Do);
    %pPHCDo  = Transform(TEMP, Q(:, :, i)'); % Local to Global
    
    % Concentration tensors
    for j = 1 : size(Closeset, 1)
        i = Closeset(j);
        pAo(:, :, :, :, i) = inversegeneral(fourthI + doubledotff(pP(:, :, :, :, i), (pC(:, :, :, :, i) ...
            -pC(:, :, :, :, phasen)))); % The matrix stiffness here is the reference stiffness tensor
    end
    pAoHCDo = fourthI; % For matrix itself
    
    prHCDo = sum(pr(Openset));
    % Average of infinite conentration tensor
    sAo= zeros(3, 3, 3, 3);
    for j = 1 : size(Closeset, 1)
        i = Closeset(j);
        sAo = sAo + pr(i)*pAo(:, :, :, :, i);
    end
    sAo = sAo + prHCDo*pAoHCDo;
    
    for j = 1 : size(Closeset, 1)
        i = Closeset(j);
        pA(:, :, :, :, i) = doubledotff(pAo(:, :, :, :, i), inversegeneral(sAo));
    end
    PAHCDo = doubledotff(pAoHCDo, inversegeneral(sAo));
    
    
    HC_D = zeros(3, 3, 3, 3); % Damaged homogenized tensor with open crack
    for j = 1 : size(Closeset, 1)
        i = Closeset(j);
        HC_D = HC_D + pr(i)*doubledotff(pC(:, :, :, :, i), pA(:, :, :, :, i));
    end
    HC_D = HC_D + prHCDo*doubledotff(HC_Do, PAHCDo);
    
    %% Second step of homogenization -- damaged matrix and biotites
    pP2 = zeros(3, 3, 3, 3, Numb + 1);
    PrepP2 = pP2;
    for i = 1 : Numb
        TEMP = Hill_P([a(1), a(1), c(1)], 16, 16, HC_D);
        pP2(:, :, :, :, i)  = Transform(TEMP, Q(:, :, i)'); % Local to GlobaL
    end
    % damaged matrix Hill's tensor (numerical intergration)
    TEMP = Hill_P([a(phasen), a(phasen), c(phasen)], 16, 16, HC_D);
    pP2(:, :, :, :, Numb + 1) = Transform(TEMP, Q(:, :, phasen)'); % Local to Global
    dpP2 = pP2 - PrepP2;
    
    
    pAo2 = zeros(3, 3, 3, 3, Numb + 1);
    for i = 1 : Numb
        pAo2(:, :, :, :, i) = inversegeneral(fourthI + doubledotff(pP2(:, :, :, :, i), (pC(:, :, :, :, i) ...
            -HC_D)));
    end
    pAo2(:, :, :, :, Numb + 1) = fourthI; % For damaged matrix
    
    pr2 = zeros(Numb + 1, 1);
    for i = 1 : Numb
        pr2(i) = pr(i);
    end
    pr2(end) = 1- sum(pr2(1 : Numb));
    
    % Average of infinite conentration tensor
    sAo2= zeros(3, 3, 3, 3);
    for i = 1: Numb + 1
        sAo2 = sAo2 + pr2(i)*pAo2(:, :, :, :, i);
    end
    
    pA2 = zeros(3, 3, 3, 3, Numb + 1);
    for i = 1 : Numb + 1
        pA2(:, :, :, :, i) = doubledotff(pAo2(:, :, :, :, i), inversegeneral(sAo2));
    end
    
    pC2 = zeros(3, 3, 3, 3, Numb + 1);
    pC2(:, :, :, :, 1:Numb) = pC(:, :, :, :, 1:Numb);
    pC2(:, :, :, :, Numb + 1) = HC_D;
    
    HC = zeros(3, 3, 3, 3); % Damaged homogenized tensor
    for i = 1 : Numb + 1
        HC = HC + pr2(i)*doubledotff(pC2(:, :, :, :, i), pA2(:, :, :, :, i));
    end
    
    % Influence tensor
    aveAoP = zeros(3, 3, 3, 3);
    for i = 1 : Numb + 1
        aveAoP = aveAoP + pr2(i)*doubledotff(pAo2(:, :, :, :, i), pP2(:, :, :, :, i));
    end
    
    avedCAoP = zeros(3, 3, 3, 3);
    for i = 1 : Numb + 1
        avedCAoP = avedCAoP + ...
            pr2(i)*(doubledotff(doubledotff((HC-pC2(:, :, :, :, i)), ...
            pAo2(:, :, :, :, i)), pP2(:, :, :, :, i)));
    end
    
    D = zeros(3, 3, 3, 3, Numb+1, Numb);
    for i = 1 : Numb + 1
        for j = 1: Numb
            TEMP1 = KronD(i, j)*fourthI - pr2(j)*pA2(:, :, :, :, i);
            TEMP2 = doubledotff(doubledotff(TEMP1, pAo2(:, :, :, :, j)), pP2(:, :, :, :, j));
            TEMP3 = doubledotff(pA2(:, :, :, :, i), aveAoP) - doubledotff(pAo2(:, :, :, :, i), pP2(:, :, :, :,i));
            TEMP4 = inversegeneral(avedCAoP);
            TEMP5 = transposefour(fourthI- pA2(:, :, :, :, j)) + ...
                doubledotff(doubledotff((HC - pC2(:, :, :, :, j)), pAo2(:, :, :, :, j)), pP2(:, :, :, :, j));
            TEMP6 = TEMP2 + doubledotff(doubledotff(TEMP3, TEMP4), pr2(j)*TEMP5);
            D(:, :, :, :, i, j) = doubledotff(TEMP6, pC2(:, :, :, :,  j));
        end
    end
    
    % REV Strain
    eigenstress = zeros(3, 3, Numb);
    avePIA = zeros(3, 3);
    for i = 1 : Numb
        eigenstress(:, :, i) = (-1)*doubledotft(pC2(:, :, :, :, i), Getastrain(:, :, i));
        avePIA = avePIA + pr2(i)*doubledottf(eigenstress(:, :, i), pA2(:, :, :, :, i));
    end
    
    % For incremental strain and stress calculation
    prestress = tstress;
    tstress = prestress + dstress;
    tstrain = doubledotft(inversegeneral(HC), (tstress - avePIA));
    
    % Phase Strain in step 1
    for j = Numb + 1 : phasen
        pstrain(:, :, j) = doubledotft(pA(:, :, :, :, j), tstrain);
        Lstrain(:, :, j) = Transform(pstrain(:, :, j), Q(:, :, j)); % Global to local
    end
    
    % Phase Strain in step 2
    pstrain2 = zeros(3, 3, Numb + 1);
    Lstrain2 = zeros(3, 3, Numb + 1);
    Deta = zeros(3, 3, phasen);
    for j = 1 : Numb + 1
        for i = 1: Numb % Only Numb inclusions have eigenstrain
            Deta(:, :, j) = doubledotft(D(:, :, :, :, j, i), Getastrain(:, :, i));
        end
        pstrain2(:, :, j) = doubledotft(pA2(:, :, :, :, j), tstrain) + Deta(:, :, j);
        Lstrain2(:, :, j) = Transform(pstrain2(:, :, j), Q(:, :, j)); % Global to local
    end
    
    %% Crack Propagation
    T_crack = zeros(3, 3, 3, 3, 37);
    for num = 1 : 37
        %                 for i = 1: 3
        %                     for j = 1: 3
        %                         for k = 1: 3
        %                             for l = 1: 3
        %                                 if Indi(num) == 1 % closed crack
        %                                     T_crack(i, j, k, l, num) = T_crack(i, j, k, l, num) + ...
        %                                         4*(1 - POISONM)/pi/(2-POISONM)*(1/2*...
        %                                         ((orientensor(j, l, num) * ...
        %                                         KronD(i, k) +  orientensor(j, k, num) ...
        %                                         * KronD(i, l)) ...
        %                                         + (orientensor(i, k, num) *...
        %                                         KronD(j, l) +  orientensor(i, l, num) ...
        %                                         * KronD(j, k))) - 2*crackdirection(i, num) ...
        %                                         * crackdirection(j, num) * crackdirection(k, num)...
        %                                         * crackdirection(l, num));
        %                                 else % open crack
        %                                     T_crack(i, j, k, l, num) = T_crack(i, j, k, l, num)...
        %                                         + 4*(1 - POISONM)/pi...
        %                                         *(POISONM/(1-2*POISONM)*...
        %                                         (orientensor(i, j, num)...
        %                                         * KronD(k, l)) + 1/(2-POISONM) *1/2*...
        %                                         ((orientensor(j, l, num) * ...
        %                                         KronD(i, k) +...
        %                                         orientensor(j, k, num) * ...
        %                                         KronD(i, l)) +...
        %                                         (orientensor(i, k, num) * ...
        %                                         KronD(j, l) +  ...
        %                                         orientensor(i, l, num) * ...
        %                                         KronD(j, k)))...
        %                                         - POISONM/(2-POISONM)*crackdirection(i, num) * ...
        %                                         crackdirection(j, num) * crackdirection(k, num) *...
        %                                         crackdirection(l, num));
        %                                 end
        %                             end
        %                         end
        %                     end
        %                 end
        T_crack(:, :, :, :, num)=X_c* inversegeneral(fourthI + doubledotff(pP(:, :, :, :, num+1), (pC(:, :, :, :, num +1) - pC(:, :, :, :, 39))));
    end
    
    
    H_crack = zeros(3, 3, 3, 3, 37);
    A_partial_rho = zeros(3, 3, 3, 3, Numb, 37);
    D_partial_rho = zeros(3, 3, 3, 3, Numb, Numb);
    g_crack = zeros(37, 1);
    f_crack = zeros(37, 1);
    for num = 1 : 37
        H_crack(:, :, :, :, num) = 4 * pi/3 * doubledotff(T_crack(:, :, :, :, num), inversegeneral(pC(:, :, :, :, phasen)));
        CD_partial_rho = -doubledotff(doubledotff(HC_D, H_crack(:, :, :, :, num)), HC_D);
        if drho(num) > 0
            TEMP1 = doubledotff(dpP2(:, :, :, :, 1)/drho(num), HC_D) ...
                + doubledotff(pP2(:, :, :, :, 1), CD_partial_rho);
        else
            TEMP1 = doubledotff(pP2(:, :, :, :, 1), CD_partial_rho);
        end
        Abo_partial_rho  = doubledotff(doubledotff(pAo2(:, :, :, :, 1), TEMP1), pAo2(:, :, :, :, 1));
        invsAo_partial_rho = - doubledotff(doubledotff(inversegeneral(sAo2), pr2(1)*Abo_partial_rho),...
            inversegeneral(sAo2));
        AD_partial_rho = invsAo_partial_rho;
        Ab_partial_rho = doubledotff(Abo_partial_rho, inversegeneral(sAo2)) +...
            doubledotff(pAo(:, :, :, :, 1), invsAo_partial_rho);
        HCD_partial_rho = -doubledotff(doubledotff(HC_D, H_crack(:, :, :, :, num)), HC_D);
        HC_partial_rho = pr2(2)*(doubledotff(HCD_partial_rho, pA2(:, :, :, :, 2)) + ...
            doubledotff(HC_D, AD_partial_rho)) + pr2(1)*Ab_partial_rho;
        
        for i = 1 : Numb
            A_partial_rho(:, :, :, :, i, num) = -4*pi/3*doubledotff(doubledotff(pA(:, :, :, :, i), ...
                T_crack(:, :, :, :, num)), ...
                inversegeneral(sAo));
            for j = 1 : Numb
                %TEMP1 = KronD(i, j)*fourthI - pr(j)*pA(:, :, :, :, i);
                %TEMPA = doubledotff(doubledotff(TEMP1, pAo(:, :, :, :, j)), pP(:, :, :, :, j));
                TEMPB = doubledotff(pA(:, :, :, :, i), aveAoP) - doubledotff(pAo(:, :, :, :, i), pP(:, :, :, :,i));
                TEMPC = inversegeneral(avedCAoP);
                TEMPD = transposefour(fourthI- pA(:, :, :, :, j)) + ...
                    doubledotff(doubledotff((HC - pC(:, :, :, :, j)), pAo(:, :, :, :, j)), pP(:, :, :, :, j));
                %TEMP6 = TEMP2 + doubledotff(doubledotff(TEMP3, TEMP4), pr(j)*TEMP5);
                TEMPA_partial_rho = - pr(j)*doubledotff(doubledotff(A_partial_rho(:, :, :, :, i, num), pAo(:, :, :, :, j)), ...
                    pP(:, :, :, :, j));
                TEMPB_partial_rho = doubledotff(pA(:, :, :, :, i), doubledotff(4*pi/3*T_crack(:, :, :, :, num), ...
                    pP(:, :, :, :, num+Numb))) + doubledotff(A_partial_rho(:, :, :, :, i, num), aveAoP);
                TEMPCinv_partial_rho = doubledotff(doubledotff((HC - pC(:, :, :, :, num + Numb)), 4*pi/3*T_crack(:, :, :, :, num)), ...
                    pP(:, :, :, :, num + Numb));
                TEMPC_partial_rho = -1*doubledotff(doubledotff(TEMPC, TEMPCinv_partial_rho), TEMPC);
                TEMPD_partial_rho = - transposefour(A_partial_rho(:, :, :, :, i, num)) + doubledotff(doubledotff(...
                    HC_partial_rho, pAo(:, :, :, :, j)), pP(:, :, :, :, num+Numb));
                D_partial_rho(:, :, :, :, i, j) = TEMPA_partial_rho + doubledotff(doubledotff(TEMPB_partial_rho, TEMPC), ...
                    pr(j)*TEMPD) + doubledotff(doubledotff(TEMPB, TEMPC_partial_rho), ...
                    pr(j)*TEMPD) + doubledotff(doubledotff(TEMPB, TEMPC), ...
                    pr(j)*TEMPD_partial_rho);
            end
        end
        
        % Energy release rate
        for i = 1: Numb
            for j = 1 : Numb
                g_crack(num) =  -1/2*doubledottt(doubledottf((tstrain), HC_partial_rho), (tstrain))...
                    + pr2(i)*doubledottt(doubledottf(eigenstress(:, :, i), Ab_partial_rho), tstrain)...
                    -1/2*pr(i)*doubledottt(doubledottf(doubledottf(eigenstress(:, :, i), D_partial_rho(:, :, :, :, i, j)), inversegeneral(pC(:, :, :, :, j))),...
                    eigenstress(:, :, j));
            end
        end
        if g_crack(num) > R_c
            g_crack(num) = R_c;
        end
        neps(num) = doubledottt(doubledotft(pA(:, :, :, :, num+Numb), tstrain), orientensor(:, :, num));
        
        if neps(num)   > 0
            Indi(num) = 0;
        else
            Indi(num) = 1;
        end
        errorf = 1;
        % Damage criterion
        Prerho = rho(num);
        f_crack(num) = g_crack(num) - (R_c*rho(num)*exp(1 - rho(num)/rho_c))/rho_c;
        if f_crack(num) > 0 && neps(num) > 0
            if abs(g_crack(num) - R_c) < 1e-8
                TEMPrho = rho_c;
            else
                f = @(x) g_crack(num) - (R_c*x*exp(1 - x/rho_c))/rho_c;
                partialf_x = @(x) (R_c*exp(-(x - rho_c)/rho_c)*(x - rho_c))/rho_c^2;
                % Newton's method (Set initial guess)
                TEMPrho = rho(num);
                iter = 1;
                while abs(errorf) > 1e-6
                    iter = iter+1;
                    if iter > 5000
                        error('Too many iterations')
                    end
                    deltaa = - real(f(TEMPrho)/partialf_x(TEMPrho));
                    TEMPrho = TEMPrho + deltaa;
                    errorf = real(f(TEMPrho));
                end
            end
            rho(num) = TEMPrho;
        end
    end
    drho(num) = rho(num) - Prerho;
    %  Recalculate the stiffness tensors for cracks
    for i = (Numb + 1) : phasen - 1
        pC(:, :, :, :, i) =  Indi(i - Numb)*3*mbm*fourthJ;
    end
    %% Save Matrices
    
    omega = zeros(3,3);
    for num = 1 : 37
        omega = omega + crackweight(num)*rho(num)* orientensor(:, :, num);
    end
end

omega_save(:, :, 1) = omega;

%% Weathering time evolution
for t = 1 : wtime
    omega_save(:, :, t) = omega;
    Sigma(:, :, t) = tstress;
    Eps(:, :, t) = tstrain;
    HC_save(:, :, :, :, t) = HC;
    %%
    % Weathering rate
    T = inter*t;
    %wrct = 10^(-12.32)*(T/(24*30*3600*12))^(-0.603);
    Surface = ((a(1)^(3.2) + 2*(a(1)*c(1))^(1.6))/3)^(1/1.6)*4*pi;
    Volume = 4/3*pi*a(1)^2*c(1);
    vm = 2.10e-4;
    V_bw = vm*Surface*10^(-12.32)*2.5189*T^(0.397)/((24*30*3600*12)^(-0.603));
    
    %V_bw = vm*Surface*wrct*T;
    for i = 1 : Numb
        wr(i) = V_bw/Volume;
        Letastrain(3, 3, i) = 0.4*wr(i)/1.4;
        Getastrain(:, :, i) = transpose(Q(:, :, i))*Letastrain(:, :, i)*Q(:, :, i) ;
        sm(i) = (1 - wr(i)) * bsm + vsm*wr(i);
        bm(i) = (1 - wr(i)) * bbm + vbm*wr(i);
        pr(i) = 4/3*pi * N_bio/Numb * a(i)^2 *c(i);
    end
    
    pro_cracks = 4/3*pi*X_c*rho/Vrev; % volume fraction of cracks are small
    pr(Numb+1 : phasen-1) = pro_cracks;
    pr(phasen) = 1 - sum(pr(1:phasen-1));
    
    % Stiffness tensors for Biotite and Matrix
    for i = 1 : Numb
        pC(:, :, :, :, i) = matrixc(bm(i), sm(i));
    end
    pC(:, :, :, :, phasen)  = matrixc(bm(phasen), sm(phasen));
    
    % Stiffness tensors for cracks
    for i = (Numb + 1) : phasen - 1
        pC(:, :, :, :, i) =  Indi(i - Numb)*3*mbm*fourthJ;
    end
    
    %% First step of homogenization -- matrix and cracks
    % Hill's tensors
    Closeset = []; %find(Indi > 0) + Numb;
    Openset = [setdiff((2:38)', Closeset); phasen];
    
    pP = zeros(3, 3, 3, 3, phasen);
    for j = 1 : size(Openset, 1)
        i = Openset(j);
        TEMP = Hill(a(i), c(i), bm(phasen), sm(phasen)); % Since the matrix is isotropic, there is an explicit formulation
        pP(:, :, :, :, i)  = Transform(TEMP, Q(:, :, i)'); % Local to Global
    end
    
    % Concentration tensors
    pAo = zeros(3, 3, 3, 3, phasen);
    for j = 1 : size(Openset, 1)
        i = Openset(j);
        pAo(:, :, :, :, i) = inversegeneral(fourthI + doubledotff(pP(:, :, :, :, i), (pC(:, :, :, :, i) ...
            -pC(:, :, :, :, phasen)))); % The matrix stiffness here is the reference stiffness tensor
    end
    pAo(:, :, :, :, phasen) = fourthI; % For matrix itself
    
    % Average of infinite conentration tensor
    sAo= zeros(3, 3, 3, 3);
    for j = 1 : size(Openset, 1)
        i = Openset(j);
        sAo = sAo + pr(i)*pAo(:, :, :, :, i);
    end
    
    pA = zeros(3, 3, 3, 3, phasen);
    for j = 1 : size(Openset, 1)
        i = Openset(j);
        pA(:, :, :, :, i) = doubledotff(pAo(:, :, :, :, i), inversegeneral(sAo));
    end
    
    HC_Do = zeros(3, 3, 3, 3); % Damaged homogenized tensor with open crack
    for j = 1 : size(Openset, 1)
        i = Openset(j);
        HC_Do = HC_Do + pr(i)*doubledotff(pC(:, :, :, :, i), pA(:, :, :, :, i));
    end
    
    %% closed cracks
    for j = 1 : size(Closeset, 1)
        i = Closeset(j);
        TEMP = Hill_P([a(i), a(i), c(i)], 16, 16, HC_Do);
        pP(:, :, :, :, i)  = Transform(TEMP, Q(:, :, i)'); % Local to Global
    end
    %TEMP = Hill_P([a(phasen), a(phasen), c(phasen)], 16, 16, HC_Do);
    %pPHCDo  = Transform(TEMP, Q(:, :, i)'); % Local to Global
    
    % Concentration tensors
    for j = 1 : size(Closeset, 1)
        i = Closeset(j);
        pAo(:, :, :, :, i) = inversegeneral(fourthI + doubledotff(pP(:, :, :, :, i), (pC(:, :, :, :, i) ...
            -pC(:, :, :, :, phasen)))); % The matrix stiffness here is the reference stiffness tensor
    end
    pAoHCDo = fourthI; % For matrix itself
    
    prHCDo = sum(pr(Openset));
    % Average of infinite conentration tensor
    sAo= zeros(3, 3, 3, 3);
    for j = 1 : size(Closeset, 1)
        i = Closeset(j);
        sAo = sAo + pr(i)*pAo(:, :, :, :, i);
    end
    sAo = sAo + prHCDo*pAoHCDo;
    
    for j = 1 : size(Closeset, 1)
        i = Closeset(j);
        pA(:, :, :, :, i) = doubledotff(pAo(:, :, :, :, i), inversegeneral(sAo));
    end
    PAHCDo = doubledotff(pAoHCDo, inversegeneral(sAo));
    
    
    HC_D = zeros(3, 3, 3, 3); % Damaged homogenized tensor with open crack
    for j = 1 : size(Closeset, 1)
        i = Closeset(j);
        HC_D = HC_D + pr(i)*doubledotff(pC(:, :, :, :, i), pA(:, :, :, :, i));
    end
    HC_D = HC_D + prHCDo*doubledotff(HC_Do, PAHCDo);
    
    %% Second step of homogenization -- damaged matrix and biotites
    pP2 = zeros(3, 3, 3, 3, Numb + 1);
    PrepP2 = pP2;
    for i = 1 : Numb
        TEMP = Hill_P([a(1), a(1), c(1)], 16, 16, HC_D);
        pP2(:, :, :, :, i)  = Transform(TEMP, Q(:, :, i)'); % Local to GlobaL
    end
    % damaged matrix Hill's tensor (numerical intergration)
    TEMP = Hill_P([a(phasen), a(phasen), c(phasen)], 16, 16, HC_D);
    pP2(:, :, :, :, Numb + 1) = Transform(TEMP, Q(:, :, phasen)'); % Local to Global
    dpP2 = pP2 - PrepP2;
    
    
    pAo2 = zeros(3, 3, 3, 3, Numb + 1);
    for i = 1 : Numb
        pAo2(:, :, :, :, i) = inversegeneral(fourthI + doubledotff(pP2(:, :, :, :, i), (pC(:, :, :, :, i) ...
            -HC_D)));
    end
    pAo2(:, :, :, :, Numb + 1) = fourthI; % For damaged matrix
    
    pr2 = zeros(Numb + 1, 1);
    for i = 1 : Numb
        pr2(i) = pr(i);
    end
    pr2(end) = 1- sum(pr2(1 : Numb));
    
    % Average of infinite conentration tensor
    sAo2= zeros(3, 3, 3, 3);
    for i = 1: Numb + 1
        sAo2 = sAo2 + pr2(i)*pAo2(:, :, :, :, i);
    end
    
    pA2 = zeros(3, 3, 3, 3, Numb + 1);
    for i = 1 : Numb + 1
        pA2(:, :, :, :, i) = doubledotff(pAo2(:, :, :, :, i), inversegeneral(sAo2));
    end
    
    pC2 = zeros(3, 3, 3, 3, Numb + 1);
    pC2(:, :, :, :, 1:Numb) = pC(:, :, :, :, 1:Numb);
    pC2(:, :, :, :, Numb + 1) = HC_D;
    
    HC = zeros(3, 3, 3, 3); % Damaged homogenized tensor
    for i = 1 : Numb + 1
        HC = HC + pr2(i)*doubledotff(pC2(:, :, :, :, i), pA2(:, :, :, :, i));
    end
    
    % Influence tensor
    aveAoP = zeros(3, 3, 3, 3);
    for i = 1 : Numb + 1
        aveAoP = aveAoP + pr2(i)*doubledotff(pAo2(:, :, :, :, i), pP2(:, :, :, :, i));
    end
    
    avedCAoP = zeros(3, 3, 3, 3);
    for i = 1 : Numb + 1
        avedCAoP = avedCAoP + ...
            pr2(i)*(doubledotff(doubledotff((HC-pC2(:, :, :, :, i)), ...
            pAo2(:, :, :, :, i)), pP2(:, :, :, :, i)));
    end
    
    D = zeros(3, 3, 3, 3, Numb+1, Numb);
    for i = 1 : Numb + 1
        for j = 1: Numb
            TEMP1 = KronD(i, j)*fourthI - pr2(j)*pA2(:, :, :, :, i);
            TEMP2 = doubledotff(doubledotff(TEMP1, pAo2(:, :, :, :, j)), pP2(:, :, :, :, j));
            TEMP3 = doubledotff(pA2(:, :, :, :, i), aveAoP) - doubledotff(pAo2(:, :, :, :, i), pP2(:, :, :, :,i));
            TEMP4 = inversegeneral(avedCAoP);
            TEMP5 = transposefour(fourthI- pA2(:, :, :, :, j)) + ...
                doubledotff(doubledotff((HC - pC2(:, :, :, :, j)), pAo2(:, :, :, :, j)), pP2(:, :, :, :, j));
            TEMP6 = TEMP2 + doubledotff(doubledotff(TEMP3, TEMP4), pr2(j)*TEMP5);
            D(:, :, :, :, i, j) = doubledotff(TEMP6, pC2(:, :, :, :,  j));
        end
    end
    
    % REV Strain
    eigenstress = zeros(3, 3, Numb);
    avePIA = zeros(3, 3);
    for i = 1 : Numb
        eigenstress(:, :, i) = (-1)*doubledotft(pC2(:, :, :, :, i), Getastrain(:, :, i));
        avePIA = avePIA + pr2(i)*doubledottf(eigenstress(:, :, i), pA2(:, :, :, :, i));
    end
    
    % For incremental strain and stress calculation
    prestress = tstress;
    switch BC
        case 1
            tstrain = doubledotft(inversegeneral(HC), (tstress - avePIA));
        case 2
            tempsig = tstress - avePIA;
            ST = tensor2matrix(inversegeneral(HC));
            tempsig(1,1) = (ST(1,2)*ST(2,3)-ST(2,2)*ST(1,3))/(ST(1,1)*ST(2,2)-ST(1,2)*ST(2,1))...
                *tempsig(3,3);
            tempsig(2,2) = (ST(2,3)*ST(1,1)-ST(2,1)*ST(1,3))/(ST(2,1)*ST(1,2)-ST(1,1)*ST(2,2))...
                *tempsig(3,3);
            tstress = tempsig + avePIA;
            tstrain = doubledotft(inversegeneral(HC), (tstress - avePIA));
            tstrain(1,1) = 0;
            tstrain(2,2) = 0;
        case 3
            tstrainrate = -1e-8;
            tstrain(3, 3) = tstrain(3, 3) + tstrainrate;
            ST = tensor2matrix((HC));
            tstrain(1, 1) = ((ST(2,2)*ST(1,3)-ST(1,2)*ST(2,3))* tstrain(3, 3) + ...
                ST(2,2)*(avePIA(1,1) - prestress(1,1)) + ST(1,2)*(prestress(2,2) -avePIA(2,2)))...
                /(ST(2,1)*ST(1,2)-ST(1,1)*ST(2,2));
            tstrain(2, 2) = ((ST(2,1)*ST(1,3)-ST(1,1)*ST(2,3))* tstrain(3, 3) +...
                ST(2,1)*(avePIA(1,1) - prestress(1,1)) + ST(1,1)*(prestress(2,2) -avePIA(2,2)))...
                /(ST(1,1)*ST(2,2)-ST(2,1)*ST(1,2)) ;
            tstress = doubledotft(HC, tstrain) + avePIA;
    end

    % Phase Strain in step 1
    for j = Numb + 1 : phasen
        pstrain(:, :, j) = doubledotft(pA(:, :, :, :, j), tstrain);
        Lstrain(:, :, j) = Transform(pstrain(:, :, j), Q(:, :, j)); % Global to local
    end
    
    % Phase Strain in step 2
    pstrain2 = zeros(3, 3, Numb + 1);
    Lstrain2 = zeros(3, 3, Numb + 1);
    Deta = zeros(3, 3, phasen);
    for j = 1 : Numb + 1
        for i = 1: Numb % Only Numb inclusions have eigenstrain
            Deta(:, :, j) = doubledotft(D(:, :, :, :, j, i), Getastrain(:, :, i));
        end
        pstrain2(:, :, j) = doubledotft(pA2(:, :, :, :, j), tstrain) + Deta(:, :, j);
        Lstrain2(:, :, j) = Transform(pstrain2(:, :, j), Q(:, :, j)); % Global to local
    end
    
    %% Crack Propagation
    T_crack = zeros(3, 3, 3, 3, 37);
    for num = 1 : 37
        %                 for i = 1: 3
        %                     for j = 1: 3
        %                         for k = 1: 3
        %                             for l = 1: 3
        %                                 if Indi(num) == 1 % closed crack
        %                                     T_crack(i, j, k, l, num) = T_crack(i, j, k, l, num) + ...
        %                                         4*(1 - POISONM)/pi/(2-POISONM)*(1/2*...
        %                                         ((orientensor(j, l, num) * ...
        %                                         KronD(i, k) +  orientensor(j, k, num) ...
        %                                         * KronD(i, l)) ...
        %                                         + (orientensor(i, k, num) *...
        %                                         KronD(j, l) +  orientensor(i, l, num) ...
        %                                         * KronD(j, k))) - 2*crackdirection(i, num) ...
        %                                         * crackdirection(j, num) * crackdirection(k, num)...
        %                                         * crackdirection(l, num));
        %                                 else % open crack
        %                                     T_crack(i, j, k, l, num) = T_crack(i, j, k, l, num)...
        %                                         + 4*(1 - POISONM)/pi...
        %                                         *(POISONM/(1-2*POISONM)*...
        %                                         (orientensor(i, j, num)...
        %                                         * KronD(k, l)) + 1/(2-POISONM) *1/2*...
        %                                         ((orientensor(j, l, num) * ...
        %                                         KronD(i, k) +...
        %                                         orientensor(j, k, num) * ...
        %                                         KronD(i, l)) +...
        %                                         (orientensor(i, k, num) * ...
        %                                         KronD(j, l) +  ...
        %                                         orientensor(i, l, num) * ...
        %                                         KronD(j, k)))...
        %                                         - POISONM/(2-POISONM)*crackdirection(i, num) * ...
        %                                         crackdirection(j, num) * crackdirection(k, num) *...
        %                                         crackdirection(l, num));
        %                                 end
        %                             end
        %                         end
        %                     end
        %                 end
        T_crack(:, :, :, :, num)=X_c* inversegeneral(fourthI + doubledotff(pP(:, :, :, :, num+1), (pC(:, :, :, :, num +1) - pC(:, :, :, :, 39))));
    end
    
    
    H_crack = zeros(3, 3, 3, 3, 37);
    A_partial_rho = zeros(3, 3, 3, 3, Numb, 37);
    D_partial_rho = zeros(3, 3, 3, 3, Numb, Numb);
    g_crack = zeros(37, 1);
    f_crack = zeros(37, 1);
    for num = 1 : 37
        H_crack(:, :, :, :, num) = 4 * pi/3 * doubledotff(T_crack(:, :, :, :, num), inversegeneral(pC(:, :, :, :, phasen)));
        CD_partial_rho = -doubledotff(doubledotff(HC_D, H_crack(:, :, :, :, num)), HC_D);
        if drho(num) > 0
            TEMP1 = doubledotff(dpP2(:, :, :, :, 1)/drho(num), HC_D) ...
                + doubledotff(pP2(:, :, :, :, 1), CD_partial_rho);
        else
            TEMP1 = doubledotff(pP2(:, :, :, :, 1), CD_partial_rho);
        end
        Abo_partial_rho  = doubledotff(doubledotff(pAo2(:, :, :, :, 1), TEMP1), pAo2(:, :, :, :, 1));
        invsAo_partial_rho = - doubledotff(doubledotff(inversegeneral(sAo2), pr2(1)*Abo_partial_rho),...
            inversegeneral(sAo2));
        AD_partial_rho = invsAo_partial_rho;
        Ab_partial_rho = doubledotff(Abo_partial_rho, inversegeneral(sAo2)) +...
            doubledotff(pAo(:, :, :, :, 1), invsAo_partial_rho);
        HCD_partial_rho = -doubledotff(doubledotff(HC_D, H_crack(:, :, :, :, num)), HC_D);
        HC_partial_rho = pr2(2)*(doubledotff(HCD_partial_rho, pA2(:, :, :, :, 2)) + ...
            doubledotff(HC_D, AD_partial_rho)) + pr2(1)*Ab_partial_rho;
        
        for i = 1 : Numb
            A_partial_rho(:, :, :, :, i, num) = -4*pi/3*doubledotff(doubledotff(pA(:, :, :, :, i), ...
                T_crack(:, :, :, :, num)), ...
                inversegeneral(sAo));
            for j = 1 : Numb
                %TEMP1 = KronD(i, j)*fourthI - pr(j)*pA(:, :, :, :, i);
                %TEMPA = doubledotff(doubledotff(TEMP1, pAo(:, :, :, :, j)), pP(:, :, :, :, j));
                TEMPB = doubledotff(pA(:, :, :, :, i), aveAoP) - doubledotff(pAo(:, :, :, :, i), pP(:, :, :, :,i));
                TEMPC = inversegeneral(avedCAoP);
                TEMPD = transposefour(fourthI- pA(:, :, :, :, j)) + ...
                    doubledotff(doubledotff((HC - pC(:, :, :, :, j)), pAo(:, :, :, :, j)), pP(:, :, :, :, j));
                %TEMP6 = TEMP2 + doubledotff(doubledotff(TEMP3, TEMP4), pr(j)*TEMP5);
                TEMPA_partial_rho = - pr(j)*doubledotff(doubledotff(A_partial_rho(:, :, :, :, i, num), pAo(:, :, :, :, j)), ...
                    pP(:, :, :, :, j));
                TEMPB_partial_rho = doubledotff(pA(:, :, :, :, i), doubledotff(4*pi/3*T_crack(:, :, :, :, num), ...
                    pP(:, :, :, :, num+Numb))) + doubledotff(A_partial_rho(:, :, :, :, i, num), aveAoP);
                TEMPCinv_partial_rho = doubledotff(doubledotff((HC - pC(:, :, :, :, num + Numb)), 4*pi/3*T_crack(:, :, :, :, num)), ...
                    pP(:, :, :, :, num + Numb));
                TEMPC_partial_rho = -1*doubledotff(doubledotff(TEMPC, TEMPCinv_partial_rho), TEMPC);
                TEMPD_partial_rho = - transposefour(A_partial_rho(:, :, :, :, i, num)) + doubledotff(doubledotff(...
                    HC_partial_rho, pAo(:, :, :, :, j)), pP(:, :, :, :, num+Numb));
                D_partial_rho(:, :, :, :, i, j) = TEMPA_partial_rho + doubledotff(doubledotff(TEMPB_partial_rho, TEMPC), ...
                    pr(j)*TEMPD) + doubledotff(doubledotff(TEMPB, TEMPC_partial_rho), ...
                    pr(j)*TEMPD) + doubledotff(doubledotff(TEMPB, TEMPC), ...
                    pr(j)*TEMPD_partial_rho);
            end
        end
        
        % Energy release rate
        for i = 1: Numb
            for j = 1 : Numb
                g_crack(num) =  -1/2*doubledottt(doubledottf((tstrain), HC_partial_rho), (tstrain))...
                    + pr2(i)*doubledottt(doubledottf(eigenstress(:, :, i), Ab_partial_rho), tstrain)...
                    -1/2*pr(i)*doubledottt(doubledottf(doubledottf(eigenstress(:, :, i), D_partial_rho(:, :, :, :, i, j)), inversegeneral(pC(:, :, :, :, j))),...
                    eigenstress(:, :, j));
            end
        end
        First(num, t) =  -1/2*doubledottt(doubledottf((tstrain), HC_partial_rho), (tstrain));
        Second(num, t) =  pr2(i)*doubledottt(doubledottf(eigenstress(:, :, i), Ab_partial_rho), tstrain);
        Third(num, t) = -1/2*pr(i)*doubledottt(doubledottf(doubledottf(eigenstress(:, :, i), D_partial_rho(:, :, :, :, i, j)), inversegeneral(pC(:, :, :, :, j))),...
            eigenstress(:, :, j));
                
        if g_crack(num) > R_c
            g_crack(num) = R_c;
        end
        neps(num) = doubledottt(doubledotft(pA(:, :, :, :, num+Numb), tstrain), orientensor(:, :, num));
        
        if neps(num)   > 0
            Indi(num) = 0;
        else
            Indi(num) = 1;
        end
        errorf = 1;
        % Damage criterion
        Prerho = rho(num);
        f_crack(num) = g_crack(num) - (R_c*rho(num)*exp(1 - rho(num)/rho_c))/rho_c;
        if f_crack(num) > 0 && neps(num) > 0
            if abs(g_crack(num) - R_c) < 1e-8
                TEMPrho = rho_c;
            else
                f = @(x) g_crack(num) - (R_c*x*exp(1 - x/rho_c))/rho_c;
                partialf_x = @(x) (R_c*exp(-(x - rho_c)/rho_c)*(x - rho_c))/rho_c^2;
                % Newton's method (Set initial guess)
                TEMPrho = rho(num);
                iter = 1;
                while abs(errorf) > 1e-6
                    iter = iter+1;
                    if iter > 5000
                        error('Too many iterations')
                    end
                    deltaa = - real(f(TEMPrho)/partialf_x(TEMPrho));
                    TEMPrho = TEMPrho + deltaa;
                    errorf = real(f(TEMPrho));
                end
            end
            rho(num) = TEMPrho;
        end
    end
    drho(num) = rho(num) - Prerho;
    %  Recalculate the stiffness tensors for cracks
    for i = (Numb + 1) : phasen - 1
        pC(:, :, :, :, i) =  Indi(i - Numb)*3*mbm*fourthJ;
    end
    %% Save Matrices
    omega = zeros(3,3);
    for num = 1 : 37
        omega = omega + crackweight(num)*rho(num)* orientensor(:, :, num);
    end
end

save(['BC',num2str(BC),'theta',num2str(theta(1)*180/pi), 'psi',num2str(psi(1)*180/pi),'D',num2str(depth),...
    'R_c',num2str(R_c),'rho_c',num2str(rho_c), 'omega_save.mat'],'omega_save');
save(['BC',num2str(BC),'theta',num2str(theta(1)*180/pi), 'psi',num2str(psi(1)*180/pi),'D',num2str(depth),...
    'R_c',num2str(R_c),'rho_c',num2str(rho_c), 'Sigma.mat'],'Sigma');
save(['BC',num2str(BC),'theta',num2str(theta(1)*180/pi), 'psi',num2str(psi(1)*180/pi),'D',num2str(depth),...
    'R_c',num2str(R_c),'rho_c',num2str(rho_c), 'Eps.mat'],'Eps');
save(['BC',num2str(BC),'theta',num2str(theta(1)*180/pi), 'psi',num2str(psi(1)*180/pi),'D',num2str(depth),...
    'R_c',num2str(R_c),'rho_c',num2str(rho_c), 'HC_save.mat'],'HC_save');

%%
close all
time = size(omega_save, 3);
omega_1 = zeros(t, 1);
omega_2 = zeros(t, 1);
omega_3 = zeros(t, 1);
for i = 1:time
    omega_1(i) = omega_save(1, 1, i );
    omega_2(i) = omega_save(2, 2, i );
    omega_3(i) = omega_save(3, 3, i );
end
figure(1)
semilogx((1:time)/(24*30*3600*12/inter), (omega_1),'b', 'LineWidth',3);
hold on;
semilogx((1:time)/(24*30*3600*12/inter), (omega_2),'r','LineWidth',3);
hold on;
semilogx((1:time)/(24*30*3600*12/inter), (omega_3),'k','LineWidth',3);
hold on;
legend('\Omega_{11}','\Omega_{22}','\Omega_{33}','Location','Best','Orientation','vertical')
xlabel('Time','FontSize',12)
ylabel('Damage','FontSize',12)
set(gca,'Fontname', 'Times New Roman','FontSize', 12);
grid on
set(gca,'FontSize',12)


Epsv =  zeros(t, 1);
for i = 1:time
    Epsv(i) = trace(Eps(:, :, i ));
end
figure(2)
semilogx((1:time)/(24*30*3600/inter), Epsv,'b', 'LineWidth',3);

figure(3)
set(gcf, 'Units', 'centimeter', 'Position', [5 12 12.6 12] )
semilogx(1 : time, First(3, :),'Color', [0.2990    0.4717    0.7293], 'LineWidth',3);
hold on
semilogx(1 : time, Second(3, :),'Color', [0.0996    0.2587    0.4163], 'LineWidth',3);
hold on
semilogx(1 : time, Third(3, :), 'Color', [0.5724    0.6524    0.7896], 'LineWidth',3);
legend('Macroscopic Strain Term','Eigenstress Term','Coupled Term','Location','Best','Orientation','vertical')
xlabel('Time (years)','FontSize',12)
ylabel('Driving Force Components','FontSize',12)
set(gca,'Fontname', 'Times New Roman','FontSize', 12);
set(gca,'FontSize',12)
set(gcf,'renderer','Painters')
print('-depsc2', 'C:\Users\txu90\Dropbox (GaTech)\Vector-based figure\Aligned')

figure(4)
set(gcf, 'Units', 'centimeter', 'Position', [5 12 12.6 12] )
semilogx(1 : time, First(2, :),'Color', [0.2990    0.4717    0.7293], 'LineWidth',3);
hold on
semilogx(1 : time, Second(2, :),'Color', [0.0996    0.2587    0.4163], 'LineWidth',3);
hold on
semilogx(1 : time, Third(2, :), 'Color', [0.5724    0.6524    0.7896], 'LineWidth',3);
legend('Macroscopic Strain Term','Eigenstress Term','Coupled Term','Location','Best','Orientation','vertical')
xlabel('Time (years)','FontSize',12)
ylabel('Driving Force Components','FontSize',12)
set(gca,'Fontname', 'Times New Roman','FontSize', 12);
set(gca,'FontSize',12)
set(gcf,'renderer','Painters')
print('-depsc2', 'C:\Users\txu90\Dropbox (GaTech)\Vector-based figure\NonAligned')
end