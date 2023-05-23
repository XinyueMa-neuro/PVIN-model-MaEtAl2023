%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example code to reproduce the two-parameter bifurcation diagram 
% shown in Fig. 6 in the manuscript in preparation:
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
% visualization setting
PICTURE_WIDTH = 17.6; % cm
PICTURE_HEIGHT = PICTURE_WIDTH*1; % always play with the figure height
fig6 = figure('Units','centimeters','Position',[3 3 PICTURE_WIDTH PICTURE_HEIGHT]);
%
r = [0  60]; 
global cds
p=[0.6, r(end)]'; ap=[1]; 
init = [-55.2239090010468	0.978799646708047	0.0587410956247739		1.21580141205003e-05]';
[x0,v0]=init_EP_EP(@PVIN_Cai,init,p,ap);
opt=contset;
opt=contset(opt,'MaxStepSize',0.025);
opt=contset(opt,'MaxNumPoints',10000);
opt=contset(opt,'Singularities',1);
opt=contset(opt,'Userfunctions',0);
[x,v,s,h,f]=cont(@equilibrium,x0,[],opt);
[x,v,s,h,f]=cont(x,v,s,h,f,cds);

% extend HB2
x1=x(1:4, s(2).index);
p=[x(end,s(2).index); r(end)];
[x0,v0]=init_H_H(@PVIN_Cai,x1,p,[1 2]);
opt=contset;
opt=contset(opt,'Singularities',1);
opt=contset(opt,'MaxStepsize',0.1);
opt=contset(opt,'MaxNumPoints',10000);
opt=contset(opt,'backward',0);
[x4,v4,s4,h4,f4]=cont(@hopf,x0,v0,opt);
pHB2_1 = plot(x4(5,:),x4(6,:),'Color','k','LineStyle',"-",'LineWidth',1); hold on; plotlabel(s4, x4)
GH1 = x4([5,6],s4(2).index);
GH2 = x4([5,6],s4(3).index);
opt=contset(opt,'backward',1); opt=contset(opt,'MaxNumPoints',10000);
[x4,v4,s4,h4,f4]=cont(@hopf,x0,v0,opt); 
pHB2_2 = plot(x4(5,x4(5,:)>-0.5),x4(6,x4(5,:)>-0.5),'Color','k','LineStyle',"-",'LineWidth',1,'HandleVisibility',"off"); plotlabel(s4, x4);
% The HB2 curve on the left side is V<-60 and is not covered in codim1
% bifurcation diagram, thus is also not shown in the codim2 diagram

% ================ LP =================
% extend LP
x1=x(1:4, s(2).index);
p=[x(end,s(2).index); r(end)];
[x0,v0]=init_LP_LP(@PVIN_Cai,x1,p,[1 2]);
opt=contset;
opt=contset(opt,'Singularities',1);
opt=contset(opt,'MaxStepsize',0.1);
opt=contset(opt,'MaxNumPoints',10000);
opt=contset(opt,'backward',0);
[x4,v4,s4,h4,f4]=cont(@limitpoint,x0,v0,opt);
CP = x4([5,6],s4(2).index);
pLP_1 = plot(x4(5,:),x4(6,:),'Color','k','LineStyle',"-.",'LineWidth',1); plotlabel(s4, x4)
opt=contset(opt,'backward',1);
opt=contset(opt,'MaxNumPoints',3000);
[x4,v4,s4,h4,f4]=cont(@limitpoint,x0,v0,opt);
BT = x4([5,6],s4(2).index);
pLP_2 = plot(x4(5,:),x4(6,:),'Color','k','LineStyle',"-.",'LineWidth',1,'HandleVisibility',"off"); plotlabel(s4, x4)
% -- visualization
xlabel('[Ca^{2+}]_i (\muM)','FontSize',10); ylabel('I_{app} (pA)','FontSize',10)

xlim([0, 1.1]);
ylim([32,580]); yticks([35 55 100 280 500])
hAxis = gca;
hAxis.YScale = 'log';
hAxis.FontSize = 14;

% ===========================================
% Marker on 2D two-par bifurcation
% ============= Fill the region =============
% color
SE_area = [1,1,0.5];
UE_area = [1,1,1];
SP_area = [0.6 1 0.1];

% -- LP1 line
[~, ind] = sort(pLP_1.YData);
% the lower line
LP1_x = pLP_1.XData(ind); 
LP1_y = pLP_1.YData(ind); 

LP_x = [flip(pLP_2.XData) pLP_1.XData];
LP_y = [flip(pLP_2.YData) pLP_1.YData];
LP_x = LP_x(LP_y > 0.1);
LP_y = LP_y(LP_y > 0.1);

% y at x=0 on the lower LP line 
[~, LP_0_ind] = min(abs(LP1_x));
LP_0_x = LP1_x(LP_0_ind);
LP_0_y = LP1_y(LP_0_ind);

% -- HB2 line
HB_x = [flip(pHB2_2.XData) pHB2_1.XData];
HB_y = [flip(pHB2_2.YData) pHB2_1.YData];
HB_x = HB_x(HB_y > 0.1);
HB_y = HB_y(HB_y > 0.1);

% y at the top of HB line
[~, HB_max_ind] = max(HB_y);
HB_max_x = HB_x(HB_max_ind);
HB_max_y = HB_y(HB_max_ind);
% y at I=0 on the HB 
HB_x_low = pHB2_2.XData;
HB_y_low = pHB2_2.YData;

[~, HB_0_ind] = min(abs(HB_x_low));
HB_0_x = HB_x_low(HB_0_ind);
HB_0_y = HB_y_low(HB_0_ind);

%
ind = HB_y<GH1(2) & HB_x>-0.1;
subHB_x = HB_x(ind);
subHB_y = HB_y(ind);

ind = HB_y>GH1(2) & HB_x>-0.1;
supHB_x = HB_x(ind);
supHB_y = HB_y(ind);

p_sub_HB = plot(subHB_x, subHB_y, 'Color',[1,0,1],'LineWidth',2,'LineStyle','-');
p_sup_HB = plot(supHB_x, supHB_y, 'Color',[1,0,1],'LineWidth',2,'LineStyle','--');
p_LP_line = plot(LP_x, LP_y, 'Color',[0,0,0],'LineWidth',1,'LineStyle','-');
uistack(p_LP_line,'bottom')
uistack(p_sub_HB,'bottom')
uistack(p_sup_HB,'bottom')

% -- stable equilibrium
set(gca, 'Color', SE_area);

% -- stable periodic        
p_SP = patch([-0.65 -0.5 HB_x -1],  [1 1.2 HB_y 1], ...
             SP_area,'EdgeColor','none','FaceAlpha', 0.5);

% -- unstable equilibrium     
p_UE = patch([ 0.01 LP_x 0.01], [1 LP_y 1], ...
            UE_area,'EdgeColor','none','FaceAlpha', 1);

uistack(p_UE,'bottom')
uistack(p_SP,'bottom')

% % ======== Not Physiological Relevant ==========
% patch([-1 0 0 -1], ...
%       [400, 400, 0, 0], ...
%       [0, 0, 0],'EdgeColor','none','FaceAlpha', 0.5);

% ============= state marker =============
text(0.05, 570, 'I','FontSize',14,'Color','k','HorizontalAlignment','Center',"FontWeight","bold"); % SE
text(0.05, 35, 'II','FontSize',14,'Color','k','HorizontalAlignment','Center',"FontWeight","bold"); % US
text(0.05, 400, 'III','FontSize',14,'Color','k','HorizontalAlignment','Center',"FontWeight","bold"); % SP

px = [-1, 2]; tx = 0.635;

% ============= firing activity divider =============
% plot(px, [BT(2), BT(2)], ':k', 'LineWidth', 1.2);
plot(px, [CP(2), CP(2)], ':k', 'LineWidth', 1.5);
plot(px, [GH1(2), GH1(2)], ':k', 'LineWidth', 1.5);
plot(px, [HB_0_y, HB_0_y], ':k', 'LineWidth', 1.5);
%         plot(px, [GH2(2), GH2(2)], ':k', 'LineWidth', 1.2);
%         plot(px, [GH3(2), GH3(2)], ':k', 'LineWidth', 1.2);
%         plot(px, [LP_0_y, LP_0_y], ':k', 'LineWidth', 1.2);
plot(px, [HB_max_y, HB_max_y], ':k', 'LineWidth', 1.5);

text(tx, 33.5,"quiescence",'FontSize',12,'Color','k','HorizontalAlignment','left',"FontWeight","bold",'FontName','Helvetica')
text(tx, BT(2)-8,['SNIC type firing'],'FontSize',12,'Color','k','HorizontalAlignment','left',"FontWeight","bold",'FontName','Helvetica')
text(tx, CP(2)+100,['elliptic type firing'],'FontSize',12,'Color','k','HorizontalAlignment','left',"FontWeight","bold",'FontName','Helvetica')
text(tx, GH1(2)+150,['parabolic type firing'],'FontSize',12,'Color','k','HorizontalAlignment','left',"FontWeight","bold",'FontName','Helvetica')
% text(tx, GH2(2)+15,"",'FontSize',9,'Color','k','HorizontalAlignment','left',"FontWeight","bold",'FontName','Helvetica')
% text(tx, HB_0_y,"SubHopf/Circle type",'FontSize',12,'Color','k','HorizontalAlignment','left',"FontWeight","bold",'FontName','Helvetica')

% ============= state transition =====================
annotation('doublearrow',[0.83, 0.83],[0.708 0.9081],'Head1Style','plain','Head2Style','plain', ...
    'Head1Length',6,'Head1Width',6, 'Head2Length',6, 'Head2Width',6);
text(1.045, 130, 'On-Off Switching','FontSize',13,'Color','k','HorizontalAlignment','Center',"FontWeight","bold",'FontName','Helvetica','Rotation',-90)
annotation('doublearrow',[0.83, 0.83],[0.2821,0.701],'Head1Style','plain','Head2Style','plain', ...
    'Head1Length',6,'Head1Width',6, 'Head2Length',6, 'Head2Width',6);
text(1.045, 382, 'Damped Oscillation','FontSize',13,'Color','k','HorizontalAlignment','Center',"FontWeight","bold",'FontName','Helvetica','Rotation',-90)


function [mytext, mydot] = plotlabel3(s4, x4)
%     skew = 0.03*[d(2)-d(1) d(4)-d(3) d(6)-d(5)];
    xindex=cat(1,s4.index); xindex=xindex(2:end-1);
    s4(1).label=''; s4(end).label=''; 
    
    % for label of more than 1 occurrance, add number
    if length(unique({s4(2:end-1).label})) ~= length({s4(2:end-1).label})
        [label,~,group] = unique({s4(2:end-1).label});
        for ii = 1:length(label)
            if sum(group==ii)>1
                ind = find( strcmp({s4.label}, label(ii)) == 1);
                ilabel = 1;
                for jj = ind
                    s4(jj).label = [s4(jj).label, num2str(ilabel)];
                    ilabel = ilabel + 1;
                end
            end
        end
    end
    
    mytext = text( x4(5,xindex)+0.1,x4(1,xindex), x4(6,xindex), {s4(2:end-1).label},'FontSize',8,'FontName','Helvetica');
    mydot = plot3( x4(5,xindex),x4(1,xindex), x4(6,xindex),'ro','MarkerFaceColor','r','MarkerSize',4);
end
function [mytext, mydot] = plotlabel(s4, x4)
%     skew = 0.03*[d(2)-d(1) d(4)-d(3) d(6)-d(5)];
    xindex=cat(1,s4.index); 
    xindex=xindex(2:end-1); s4 = s4(2:end-1);
%     s4(1).label=''; s4(end).label=''; 
    
    % exclude bifurcation point out of bound
    ind_inbound = find(x4(6,xindex)<0);
    if any(ind_inbound)
%         for ii = ind_inbound
%             s4(ii).label = '';
%         end
        s4(ind_inbound) = [];
        xindex(ind_inbound) = [];
    end
    
     
    % for label of more than 1 occurrance, add number
    if length(unique({s4.label})) ~= length({s4.label})
        [label,~,group] = unique({s4.label});
        for ii = 1:length(label)
            if sum(group==ii)>1
                ind = find( strcmp({s4.label}, label(ii)) == 1);
                ilabel = 1;
                for jj = ind
                    s4(jj).label = [s4(jj).label, num2str(ilabel)];
                    ilabel = ilabel + 1;
                end
            end
        end
    end
    
    mytext = text( x4(5,xindex)-0.06,x4(6,xindex), {s4.label},'FontSize',8,'FontName','Helvetica');
    mydot = plot( x4(5,xindex),x4(6,xindex),'ro','MarkerFaceColor','r','MarkerSize',4);
end
