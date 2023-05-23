function [t_model, v_model,current] = runHHmodel_STEP(r, mytype, dt)
% The simulation used the STEP current protocol same as that of experimental data
%
% ====================================================
% >>> INPUTs
% r: [Bt, gSK, ksk, Iapp]
%     uM,   nS, uM,   pA
% mytype: current stimulation protocol
%         'step': stable step current of value r(end)
% dt: time step (ms)
%
% ====================================================
% -- OUTPUT
% 't_model' (ms)  timepoints
% 'v'       (mV)  voltage trace
% ====================================================

    rtest = r(end); % -- test current
    interval = [200 1000 300];
    t_model = 0:dt:sum(interval);
    opts = odeset('RelTol',1e-4,'AbsTol',1e-5);

    % -- Initial conditions
    y0 = [ -69.9853123603117	0.999439610028142	0.000453598845522021	0.000112234556488220 0.0295977873620370];   
    v_model = []; 
    holdI = -100; % hold voltage at -70 mV
    
    % pre-holding: get the steady state
    r(end) = 0 + holdI;
    sol = ode45(@(t,y)PVIN_HH(t,y,r,'step'), 0:dt:500, y0, opts);
    y0 = sol.y(:,end); 

    % Simulation
    % I=0
    r(end) = 0 + holdI;
    sol = ode45(@(t,y)PVIN_HH(t,y,r,'step'), 0:dt:interval(1), y0, opts);
    y0 = sol.y(:,end);
    v_model = [v_model, interp1(sol.x, sol.y(1,:), 0:dt:interval(1))];
    
    % test stimuli
    if strcmp(mytype, 'step'), r(end) = rtest + holdI; else, r(end) = rtest; end
    sol = ode45(@(t,y)PVIN_HH(t,y,r,mytype), 0:dt:interval(2), y0, opts);
    y0 = sol.y(:,end);
    v_model = [v_model(1:end-1), interp1(sol.x, sol.y(1,:), 0:dt:interval(2))];

    % I=0
    r(end) = 0 + holdI;
    sol = ode45(@(t,y)PVIN_HH(t,y,r,'step'), 0:dt:interval(3), y0, opts);
    v_model = [v_model(1:end-1), interp1(sol.x, sol.y(1,:), 0:dt:interval(3))];

    % -- current
    current = zeros(size(t_model)); 
    current(t_model>=interval(1) & t_model<=interval(2)) = rtest;
end 