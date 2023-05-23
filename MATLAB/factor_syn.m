function f_syn = factor_syn(tau_d, tau_r)
% generate a factor "f_syn" so that the exponential equations peaks at 1
    t_peak = log(tau_r/tau_d) / (1/tau_d - 1/tau_r);
    f_syn = exp(t_peak/tau_d) * tau_d / (tau_d-tau_r);
end