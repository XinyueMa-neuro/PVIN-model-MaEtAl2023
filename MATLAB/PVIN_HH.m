function dydt=PVIN_HH(t, y, r, mytype)
% ====================================================
% ODEs for the parvalbumin-expressing interneuron (PVIN) model
% used in the manuscript in preparation:
% "Ma, X., Miraucourt, L., Qiu, H., Sharif-Naeini, R., Khadra, A. (2023). 
% Calcium buffering tunes intrinsic excitability of spinal dorsal horn 
% parvalbumin-expressing interneurons: A computational model."
%
% ====================================================
% >>> INPUTs
% t: time point unit in ms
% y: [V, h, n1, n3, Cai] 
%    mV, 1,  1,  1, uM
% r: [Bt, gSK, ksk, Iapp]
%     uM,   nS, uM,   pA
% mytype: current stimulation protocol
%         'step': stable step current of value r(end)
%         'ramp': ramp current increasing at a rate of r(end) nA/s
%         'syn':  synaptic current. 
%                 See the above publication "Ma et al 2023" for details
%     
% ====================================================
% >>> OUTPUTs
% dydt: [dVdt, dhdt, dn1dt, dn3dt, dCaidt];
% ====================================================

    dydt = zeros(5,1);

    % Model Parameters
    gNa = 300; VNa = 58; gKv1=15; gKv3 = 180;  VK=-80; gCa = 8; VCa=68; 
    Vleak = -50; gleak = 8; gSK = 10;  Cm = 30; pgamma = 0.01;
    % - INa
    Vm = -17.5; Sm = -11.4;
    Aah = 0.0025; Sah = 10.0; Vah = 23;
    Abh = 0.094; Sbh = -5.5; Vbh = -31;
    % - IKv1
    Aan1=0.0020000000000; Van1=-30.000000000000000; San1=-9.000000000000000;
    Abn1=0.0170000000000; Vbn1=-35.000000000000000; Sbn1=5.900000000000000;
    % - IKv3
    Aan3=1.9800000000000; Van3=96.00000000000000; San3=-12.600000000000;
    Abn3=0.34000000000; Vbn3=-36.000000000000000; Sbn3=10.5;
    % - ICa
    Va=3; Sa=-10.4;
    % - ISK
    nk = 5; ksk = 0.8; KD=0.1;
    % Ca dynamics
    F=0.096485332100000; mArea=3000.000000000000000;
    d=0.100000000000000; Car=0.070000000000000; Brest=10;  
    xppfile = 'trial26';

    Iapp = r(end);
    if length(r) == 2
        Bt = r(1); 
    elseif length(r) == 3
        Bt = r(1); gSK = r(2);
    elseif length(r) == 4
        Bt = r(1); gSK = r(2); ksk = r(3);
    else 
        return;
    end

    switch mytype
        case 'step' %% STEP
            Iapp = r(end);
        case 'ramp' %% RAMP
            Iramp = r(end); 
            Tramp = 1000; fre_data = 10; 
            Ihold = -100; Ihyper = -50;
            tspan = 0:1/fre_data:Tramp;
            current = linspace(Ihyper, Iramp, length(tspan)) + Ihold;
            Iapp = interp1(tspan, current, t, 'nearest');
        case 'syn' %% synaptic input
            global gAMPA gNMDA tStim
            % -- NMDA Mg block
            MgE = 1; Vexc = 0;
            mgbl = 1./(1+MgE/3.57*exp(-0.062*y(1)));
            gSyn = mgbl .* interp1(tStim*1000, gNMDA, t, 'nearest') + interp1(tStim*1000, gAMPA, t, 'nearest');
            Iapp = - gSyn .* ( y(1) - Vexc );
    end

    % INa
    mmax=1./(1+exp(( y(1)-Vm )/Sm));
    ah=Aah./exp(( y(1)-Vah )/Sah );
    bh=Abh*( y(1)-Vbh )./(1-exp(( y(1)-Vbh )/Sbh));
    INa = gNa.*mmax.^3.*y(2).*( y(1)-VNa );

    % IKv1
    an1=Aan1*( y(1)-Van1 )./(1-exp(( y(1)-Van1 )/San1));
    bn1=Abn1./exp(( y(1)-Vbn1 )/Sbn1);
    IKv1 = gKv1.*y(3).^4.*(y(1)-VK);
    
    % IKv3
    an3=Aan3*( y(1)-Van3 )./(1-exp(( y(1)-Van3 )/San3));
    bn3=Abn3./exp(( y(1)-Vbn3 )/Sbn3);
    IKv3 = gKv3.*y(4).^2.*(y(1)-VK);

    % ICa
    amax=1./(1+exp((y(1)-Va)./Sa));
    ICa = gCa*amax.^2.*(y(1)-VCa);

    % ISK
    k = y(5).^nk ./(ksk.^nk+y(5).^nk);
    ISK = gSK.*k.^1.*(y(1)-VK);

    % Ileak
    Ileak = gleak.*(y(1)-Vleak);

    dydt(1) = (-Ileak-INa-IKv1-IKv3-ICa-ISK+Iapp)/Cm; % v
    dydt(2) = ah*(1-y(2))-bh*y(2); % h
    dydt(3) = an1*(1-y(3))-bn1*y(3); % n1
    dydt(4) = an3*(1-y(4))-bn3*y(4); % n3
    dydt(5) = ( - ICa/2/F/mArea/d  - pgamma*(y(5)-Car) )/(1+Bt/KD); % Cai
    

end