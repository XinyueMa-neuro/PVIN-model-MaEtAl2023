%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example code to reproduce Poisson-distributed synaptic current-stimulated voltage response of PVINs
% shown in Fig. 7 in the publication:
% "Ma, X., Miraucourt, L., Qiu, H., Sharif-Naeini, R., Khadra, A. (2023). 
% Calcium buffering tunes intrinsic excitability of spinal dorsal horn 
% parvalbumin-expressing interneurons: A computational model."
%
%---------------------------------------------
% Tested Under MATLAB Version: 9.12.0 (R2022a)
% Time-stamp: <2023-Jan-17> 
%---------------------------------------------
%
% Xinyue Ma
% Email: xinyue.ma@mail.mcgill.ca
% Integrated Program in Neuroscience
% McGill University
% Montreal, QC, H3A 1A1 
% Canada
%
%-------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% parameter setting
% paramter setting
itype = input('>>> Input the condition (1-naive PVIN | 2-CCI PVIN)\n');

switch itype
    case 1 %% naive condition
        Bt = 90; sgn = 'Naive';
    case 2 %% CCI 
        Bt = 10; sgn = 'CCI';
    otherwise
        disp('error: unrecognized condition input')
        return
end
abeta_fre = input('>>> Input the presynaptic input firing rate (Hz):\n');

% itype = 1; % 1-naive PVIN | 2-CCI PVIN
% abeta_fre = 5;
% % a beta fiber firing frequency
% % 5 Hz - light touch | 15 Hz - maximal firing rate

rng(2);

r = [Bt, abeta_fre]; dt = 0.0005; % ms

% -- simulation
[t_model, v_model, gSyn] = runHHmodel_AbetaPoisson(r, 'syn');

% -- visualization
figure('Position',[0,0,700,500])
subplot 211
plot(t_model, v_model,'k'); 
xlabel('T (ms)'); ylabel('V (mV)'); 
title(['PVIN model ([B_{tot}]_i=',num2str(Bt),'\muM): synaptic current stimulation (',num2str(abeta_fre),'Hz)']);
subplot 212
plot(t_model, gSyn); 
legend({'g_{AMPA}','g_{NMDA}'},'Location','best')
xlabel('T (ms)'); ylabel('g_{syn} (nS)'); 