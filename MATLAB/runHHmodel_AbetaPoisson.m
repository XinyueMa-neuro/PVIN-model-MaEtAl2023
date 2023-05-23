function [t_model, v, gSyn] = runHHmodel_AbetaPoisson(r, mytype)
% This is a 3-s long synaptic current protocol
%
% ====================================================
% >>> INPUTs
% r: [Bt, gSK, ksk, Iapp]
%     uM,   nS, uM,   pA
% mytype: current stimulation protocol
%         'syn':  synaptic current. 
%                 See the above publication "Ma et al 2023" for details
%
% ====================================================
% >>> OUTPUTs
% 't_model' (ms)  timepoints
% 'v'       (mV)  voltage trace
% "gSyn"    (nS)  synaptic conductance [AMPA; NMDA]
%
    opts = odeset('RelTol',1e-4,'AbsTol',1e-5);

    rtest = r(end); % -- test current

    % Initial conditions
    y0 = [ -69.9853123603117	0.999439610028142	0.000453598845522021	0.000112234556488220	0.0295977873620370];
    v = []; gSyn = [];
    TSyn = 3; % simulation duration (unit in second)
    dt = 0.005; % integration time step (unit in ms)
    Ihold = 0; % hold voltage at -70 mV

    % Synaptic Input
    global gAMPA gNMDA tStim tdelay_NMDA tdelay_AMPA
    abeta_fre = r(end); 
    dt_Syn = 0.0001; % a smaller timestep to improve the synaptic array resolution
    tStim = 0:dt_Syn:TSyn;
    [gAMPA, gNMDA, tdelay_AMPA, tdelay_NMDA] = mySynInput(abeta_fre, tStim);
    
    % - I=0
    trest1 = 500;
    r(end) = 0 + Ihold; 
    sol = ode45(@(t,y)PVIN_HH(t,y,r,'step'), 0:dt:trest1, y0,opts);
    y0 = sol.y(:,end);
    v = [v, interp1(sol.x, sol.y(1,:), 0:dt:trest1, 'nearest')];
    gSyn(1,:) = zeros(size(0:dt:trest1)); % AMPA
    gSyn(2,:) = zeros(size(0:dt:trest1)); % NMDA
    
    % - test stimuli
    r(end) = rtest; ttest = TSyn*1000;
    sol = ode45(@(t,y)PVIN_HH(t,y,r,mytype), 0:dt:ttest, y0,opts);
    y0 = sol.y(:,end);
    v = [v(1:end-1), interp1(sol.x, sol.y(1,:), 0:dt:ttest, 'nearest')];
    gSyn = [gSyn(:,1:end-1), [interp1(tStim*1000, gAMPA, 0:dt:ttest, 'nearest'); interp1(tStim*1000, gNMDA, 0:dt:ttest, 'nearest')] ]; 
    
    % - I=0
    r(end) = 0 + Ihold; trest2 = 500;
    sol = ode45(@(t,y)PVIN_HH(t,y,r,'step'), 0:dt:trest2, y0,opts);
    v = [v(1:end-1), interp1(sol.x, sol.y(1,:), 0:dt:trest2, 'nearest')];
    gSyn = [gSyn, zeros(2, size(0:dt:trest2,2)-1)];

    t_model = 0:dt:(trest1+trest2+ttest);

end 