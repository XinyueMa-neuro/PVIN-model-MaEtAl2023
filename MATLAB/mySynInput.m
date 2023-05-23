function [gsyn_AMPA, gsyn_NMDA, tdelay_AMPA, tdelay_NMDA] = mySynInput(abeta_fre, tStim)
% generate time series of synaptic conductances at customized firing rate
%
% >>> INPUTs
% abeta_fre: Poisson process rate
% tStim: (s) time array
%
% >>> OUTPUTs
% gsyn_AMPA     (nS) AMPA conductance: Abeta -> iPVIN
% gsyn_NMDA     (nS) NMDA conductance: Abeta -> iPVIN
% tdelay_AMPA   (s) start timing of the synaptic conductance of AMPA
% tdelay_NMDA   (s) start timing of the synaptic conductance of NMDA
%

    Vexc = 0;
    nTrials = 1;
    weight_ampa = 2; rtau_ampa = 0.1; dtau_ampa = 5;
    weight_nmda = 3; rtau_nmda = 2; dtau_nmda = 100; 
    
    % -- independent Poisson distributed synaptic input
    ifspike = poissonSpikeGen ( abeta_fre , tStim , nTrials );
    tdelay_AMPA = tStim(ifspike); 
    tdelay_NMDA = tStim(ifspike); 
    gsyn_AMPA = genSyn(tStim, tdelay_AMPA, weight_ampa, rtau_ampa, dtau_ampa);
    gsyn_NMDA = genSyn(tStim, tdelay_NMDA, weight_nmda, rtau_nmda, dtau_nmda);
    
end




