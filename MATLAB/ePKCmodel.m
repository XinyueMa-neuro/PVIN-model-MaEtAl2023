function dydt=ePKCmodel(t,y,r,mytype)
    % adopted from: https://senselab.med.yale.edu/ModelDB/showmodel.cshtml?model=267056&file=/SDHmodel/cells.py#tabs-2
    % unit:
    %   g - S/cm^2
    %   Cm - uF/cm^2
    % Here we only coded the soma
    dydt = zeros(10,1);
    
    % -- applied current
    switch mytype
        case 'step'
            Iapp = r; 
        case 'syn'
            global gAMPA gNMDA gGABA gGlycin tStim
            % -- NMDA Mg block
            MgE = 1; Vexc = 0; Vinh = -70;
            mgbl = 1./(1+MgE/3.57*exp(-0.062*y(1)));
            gSyn_exc = mgbl .* interp1(tStim*1000, gNMDA, t, 'nearest') ...
                   + interp1(tStim*1000, gAMPA, t, 'nearest');
            gSyn_inh = interp1(tStim*1000, gGABA+gGlycin, t, 'nearest');
            Iapp = - gSyn_exc .* ( y(1) - Vexc ) - gSyn_inh .* ( y(1) - Vinh );
    end
    
    % -- soma: cylindrical compartments (cm)
    L = 20; diam = 20; Area = L*diam*pi; Cm = 1;
    VNa = 50; VK = -70; Vleak = -65; celsius = 23;  
    
    %% -- B_Na: Hodgkin - Huxley squid sodium channel
    gbna = 0.0001652;
    gbna =  0.1625;
    IBNa = gbna*y(2).^3*y(3)*(y(1)-VNa);
    [minf_BNa, mtau_BNa, hinf_BNa, htau_BNa] = gateIBNa(y(1), celsius);
    
    %% -- HH2: Fast Na+ and K+ currents responsible for action potentials
    gna = 0.08548; gk = 0.0043;
    gna = 85.48; gk = 4.3; 

    INa = gna*y(4).^3*y(5)*(y(1)-VNa);
    IK = gk*y(6).^4*(y(1)-VK);
    v = y(1);
    [minf, mtau, hinf, htau] = gateINa(y(1), celsius);
    [ninf, ntau] = gateIK(y(1), celsius);
    
    %% -- bordka: Borg-Graham type generic K-A channel for a Sympathetic Preganglionic Neuron
    gkabar = 0.01090 * 1.0; gkabar = 10.9*1e-3; gkabar =  10.9;
    gkabar = 10.9;
  
    IA = gkabar*y(7)*y(8)*(y(1)-VK);
    
    [nainf, natau, lainf, latau] = gateIA(y(1), celsius);
    
    %% -- KDRI:  Hodgkin - Huxley k channel
    gkdrbar = 0.0001110; gkdrbar = 0.111*1e-3; gkdrbar = 3.111;
    IKDR = gkdrbar*y(9).^4*y(10)*(y(1)-VK);
    [nkdrinf, nkdrtau, hkdrinf, hkdrtau] = gateIKDR(y(1), celsius);
    
    %% -- iKCa
    IKCa = 0;
    
    %% -- Ileak
    gleak = 0.96e-06;  gleak = 0.96 * 1e-3; gleak = 0.96;
    Ileak = gleak*(y(1)-Vleak);
    
%     Iapp/Area*1e-9
    dydt(1) = (-IBNa-INa-IK-IA-IKDR-IKCa-Ileak+Iapp/Area*1e2)/Cm * 1e3;
    dydt(2) = ( minf_BNa-y(2) )/mtau_BNa;
    dydt(3) = ( hinf_BNa-y(3) )/htau_BNa;
    dydt(4) = ( minf-y(4) )/mtau;
    dydt(5) = ( hinf-y(5) )/htau;
    dydt(6) = ( ninf-y(6) )/ntau; 
    dydt(7) = ( nainf-y(7) )/natau;
    dydt(8) = ( lainf-y(8) )/latau;    
    dydt(9) = ( nkdrinf-y(9) )/nkdrtau;
    dydt(10) = ( hkdrinf-y(10) )/hkdrtau;
end

function [minf, mtau, hinf, htau] = gateIBNa(v, celsius)
%     tadj = 3^((celsius - 23) / 10); 
    tadj = 1;

    ma = 0.182 * trap(-v + 7 - 35, 9);
    mb = 0.124 * trap(v - 7 + 35, 9);
    minf = ma ./ (ma+mb);
    mtau = 1 ./  (ma+mb) / tadj;
    
    hinf = 1 ./ (1 + exp(( v + 75 - 11) / 9));
    ha = 0.061 * trap(-v + 13 - 48, 3) + 0.0166;
    hb = 0.0018 * trap(v - 13 + 84, 18);
    htau = 1 ./ (ha+hb) / tadj;
end

function [minf, mtau, hinf, htau] = gateINa(v, celsius)
    tadj = 3.0 ^ ( (celsius-36)/ 10 );
%     tadj = 1;
    vtraub = -50.2; vh = 5;
    vtraub = -63.0; vh = -2;
    v2 = v - vtraub;
	ma = 0.32 * (13-v2) ./ ( exp((13-v2)/4) - 1);
	mb = 0.28 * (v2-40) ./ ( exp((v2-40)/5) - 1);
    minf = ma ./ (ma+mb);
    mtau = 1 ./  (ma+mb) / tadj;
    
    ha = 0.128 * exp((17-v2-vh)/18);
    hb = 4 ./ ( 1 + exp((40-v2-vh)/5) );
    hinf = ha ./ (ha+hb);
    htau = 1./(ha+hb) / tadj;
end

function [ninf, ntau] = gateIK(v, celsius)
    tadj = 3.0 ^ ( (celsius-36)/ 10 );
    vtraub = -50.2; 
    v2 = v - vtraub;
    na = 0.032 * (15-v2) ./ ( exp((15-v2)/5) - 1);
    nb = 0.5 * exp((10-v2)/40);
    ninf = na ./ (na+nb);
    ntau = 1 ./ (na+nb) / tadj;
end

function [nainf, natau, lainf, latau] = gateIA(v, celsius)
% 	celsius = 40;% this value is corrected so that IA can delay the neural excitation
	
    vhalfn=-45; 
    vhalfl=-67;
	vhalfm=-67; vhalfm = -76;
	vhalfk=-45; vhalfk = -15;
	a0l=0.023;
	a0n=0.04;
	zetan=-4; zetan = -3;
	zetal=2; 
	gmn=0.45;
	gml=1;
	zetam=4;
	zetak=-5; zetak = -0.4;
    q10=3^((celsius-30)/10);
    alpn = exp(1.e-3*zetan*(v-vhalfn)*9.648e4/(8.315*(273.16+celsius)));
    alpk = exp(1.e-3*zetak*(v-vhalfk)*9.648e4/(8.315*(273.16+celsius)));
    betn = exp(1.e-3*zetan*gmn*(v-vhalfn)*9.648e4/(8.315*(273.16+celsius)));
    alpl = exp(1.e-3*zetal*(v-vhalfl)*9.648e4/(8.315*(273.16+celsius)));
    alpm = exp(1.e-3*zetam*(v-vhalfm)*9.648e4/(8.315*(273.16+celsius)));
    betl = exp(1.e-3*zetal*gml*(v-vhalfl)*9.648e4/(8.315*(273.16+celsius)));
    
    nainf = 1/(1 + alpk) * 1.2;
    natau = betn/(q10*a0n*(1 + alpn)) *0.1;
    lainf = 1/(1 + alpm);
    latau = betl/(q10*a0l*(1 + alpl)) * 1.2;
end

function [nkdrinf, nkdrtau, hkdrinf, hkdrtau] = gateIKDR(v, celsius)
    tau_factor = 1;
    nkdra = tau_factor * 1 * .035*trap(-v - 15, 9);
    nkdrb = tau_factor * 1 * .014*exp((-v + 12)/46);
    hkdra = tau_factor * 0.0083*(1 ./(exp((v + 20)/10)+1) + 1);
    hkdrb = tau_factor * 0.0083 ./(exp((-v - 20)/10)+1);
    
    nkdrinf = nkdra ./ (nkdra+nkdrb);
    nkdrtau = 1 ./ (nkdra+nkdrb);
    hkdrinf = hkdra ./ (hkdra+hkdrb);
    hkdrtau = 1 ./ (hkdra+hkdrb);
end


function trap = trap(x,y)
	if abs(x/y) < 1e-6
		trap = y*(1 - x./y/2);
    else
		trap = x./(exp(x./y) - 1);
	end
end