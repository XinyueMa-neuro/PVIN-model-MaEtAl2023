function gsyn = genSyn(tStim, tdelay, weight, rtau, dtau)
% -- generate a temporal synaptic current
    % time constant msec --> sec
    rtau = rtau/1000; dtau = dtau/1000;
    fsyn = factor_syn(dtau, rtau);
    exp2syn = zeros(size(tStim));
    for ii = 1:length(tdelay)
        tt = tStim - tdelay(ii);
        grise = -exp(-tt/rtau); grise(tt<0) = 0;
        gdecay = exp(-tt/dtau); gdecay(tt<0) = 0;
        g = grise + gdecay;
        exp2syn = exp2syn + g;
    end
    gsyn = weight * fsyn * exp2syn;
end