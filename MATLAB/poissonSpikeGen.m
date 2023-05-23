function [ ifspike ] = poissonSpikeGen ( fr , tt , nTrials )
% Generate Poisson-distributed timing of spike trains
%
% >>> INPUTs
% fr: firing rate 
% tt: time array 
% nTrials: number of trials
%
% >>> OUTPUTs
% ifspike: timing of spike
%
    dt = tt(2)-tt(1); % s
    nBins = length(tt) ;
    ifspike = rand ( nTrials , nBins ) < fr * dt;
end