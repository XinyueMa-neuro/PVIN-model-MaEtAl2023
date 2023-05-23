%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example code to reproduce the step current-stimulated voltage response of PVINs 
% shown in Figs. 2-4 in the manuscript in preparation:
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

% paramter setting
itype = input('>>> Input the condition (1-naive PVIN | 2-CCI PVIN)\n');

% itype = 1; % 1-naive PVIN | 2-CCI PVIN
% iapp = 200; % applied step current value

switch itype
    case 1 %% naive condition
        Bt = 90; sgn = 'Naive';
    case 2 %% CCI 
        Bt = 10; sgn = 'CCI';
    otherwise
        disp('error: unrecognized condition input')
        return
end
iapp = input('>>> Input applied step current value (pA):\n');

%% -- model simulation
dt = 0.001;
clear t_model v_model

% - STEP protocol
r = [Bt, iapp];
[t_model, v_model,current] = runHHmodel_STEP(r,'step',dt);

% -- visualization
figure('Position',[0,0,700,500])
subplot 211
plot(t_model, v_model,'k'); 
xlabel('T (ms)'); ylabel('V (mV)'); 
title(['PVIN model ([B_{tot}]_i=',num2str(Bt),'\muM): step current stimulation']);
subplot 212
plot(t_model, current,'k'); 
xlabel('T (ms)'); ylabel('I_{app} (pA)'); 