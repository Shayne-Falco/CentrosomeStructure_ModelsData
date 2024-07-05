clear all
close all
% LabelsCells = {'ABa', 'ABp', 'EMS', 'P2', 'Unclear(0)','Unclear(1)'}; 
% LabelsPatterns = {'ABa P2', 'ABa EMS', 'ABp P2', 'ABP P2', 'Unclear(0)','Unclear(1)'}; 
% 
% ABaP2 = 20;
% ABaEMS = 16;
% ABpP2 = 12;
% ABpEMS = 14;
% Unclear_0 = 3;
% Unclear_1 = 4;
% 
% Proba = [ABaP2,ABaEMS,ABpP2,ABpEMS];
% Proba = Proba/sum(Proba);
% ABa = ABaEMS+ABaP2;
% ABp = ABpEMS+ABpP2;
% EMS = ABaEMS+ABpEMS;
% P2 = ABaP2+ABpP2;

unknownl1s = 0:0.01:1;
unknownl2s = 0:0.01:1;
PossibleL1s = zeros([length(unknownl1s) length(unknownl2s)]);
PossibleL2s = zeros([length(unknownl1s) length(unknownl2s)]);
PossibleL3s = zeros([length(unknownl1s) length(unknownl2s)]);
PossibleL4s = zeros([length(unknownl1s) length(unknownl2s)]);
PossibleL5s = zeros([length(unknownl1s) length(unknownl2s)]);
j=1;
XTitle = 'l1';
YTitle = 'l4';
for l1 = unknownl1s
    i=1;
    for l4 = unknownl2s
        %% ABa P2 - 15 ABa EMS - 35 ABp P2 - 35 ABp EMS - 15
%         TitleText = 'ABa/P2-15 ABa/EMS-35 ABp/P2-35 ABp/EMS-15';
%         PossibleL1s(i,j) = l1;
%         PossibleL2s(i,j) = -(.5000000000*(-1.+(-2.*l1+2.)*l4))/l1;
%         PossibleL3s(i,j) = -(1.*(-3.-(1.*(-20.*l1+20.))*l4*((-10.*l1+10.)*l4+3.*l1-5.)/((20.*l1-20.)*l4-10.*l1+10.)))/(10.-(10.*(-2.*l1+2.))*l4);
%         PossibleL4s(i,j) = l4;
%         PossibleL5s(i,j) = -(1.*((-10.*l1+10.)*l4+3.*l1-5.))/((20.*l1-20.)*l4-10.*l1+10.);
        %% ABa P2 - 35 ABa EMS - 35 ABp P2 - 15 ABp EMS - 15
%         TitleText = 'ABa/P2-35 ABa/EMS-35 ABp/P2-15 ABp/EMS-15';
%         PossibleL1s(i,j) = l1;
%         PossibleL2s(i,j) = -(.1000000000*(-7.+(-10.*l1+10.)*l4))/l1;
%         PossibleL3s(i,j) = -(1.*(-7.+(-10.*l1+10.)*l4))/(14.-(2.*(-10.*l1+10.))*l4);
%         PossibleL4s(i,j) = l4;
%         PossibleL5s(i,j) = .5;
        %% ABa P2 - 15 ABa EMS - 35 ABp P2 - 15 ABp EMS - 35
%         TitleText = 'ABa/P2-15 ABa/EMS-35 ABp/P2-15 ABp/EMS-35';
%         PossibleL1s(i,j) = l1;
%         PossibleL2s(i,j) = -(.5000000000*(-1.+(-2.*l1+2.)*l4))/l1;
%         PossibleL3s(i,j) = -(1.*(-3.+(-6.*l1+6.)*l4))/(10.-(10.*(-2.*l1+2.))*l4);
%         PossibleL4s(i,j) = l4;
%         PossibleL5s(i,j) = 0.3;
%         %% ABa P2 - 35 ABa EMS - 15 ABp P2 - 35 ABp EMS - 15
%         TitleText = 'ABa/P2-35 ABa/EMS-15 ABp/P2-35 ABp/EMS-15';
%         PossibleL1s(i,j) = l1;
%         PossibleL2s(i,j) = -(.5000000000*(-1.+(-2.*l1+2.)*l4))/l1;
%         PossibleL3s(i,j) = -(1.*(-7.+(-14.*l1+14.)*l4))/(10.-(10.*(-2.*l1+2.))*l4);
%         PossibleL4s(i,j) = l4;
%         PossibleL5s(i,j) = .7000000000;
        %% ABp P2 - 35 ABp EMS - 15
%         TitleText = 'ABa/P2-25 ABa/EMS-25 ABp/P2-35 ABp/EMS-15';
%         YTitle = 'l5';
%         PossibleL1s(i,j) = l1;
%         PossibleL2s(i,j) = -(1.*(-1.-(1.*(-4.*l1+4.))*l4*((-10.*l1+10.)*l4+5.*l1-6.)/((20.*l1-20.)*l4-12.*l1+12.)))/(2.400000000-(.8000000000*(-5.*l1+5.))*l4);
%         PossibleL3s(i,j) = -(.2000000000*(-3.+(-5.*l1+5.)*l4))/l1;
%         PossibleL4s(i,j) = -(1.*((-10.*l1+10.)*l4+5.*l1-6.))/((20.*l1-20.)*l4-12.*l1+12.);
%         PossibleL5s(i,j) = l4;
        %% ABa P2 - 35 ABa EMS - 15
%         TitleText = 'ABa/P2-35 ABa/EMS-15 ABp/P2-25 ABp/EMS-25';
%         PossibleL1s(i,j) = l1;
%         PossibleL2s(i,j) = -(.5000000000*(-1.+(-2.*l1+2.)*l4))/l1;
%         PossibleL3s(i,j) = -(1.*(-7.-(1.*(-20.*l1+20.))*l4*((-12.*l1+12.)*l4+7.*l1-6.)/((20.*l1-20.)*l4-10.*l1+10.)))/(10.-(10.*(-2.*l1+2.))*l4);
%         PossibleL4s(i,j) = l4;
%         PossibleL5s(i,j) = -(1.*((-12.*l1+12.)*l4+7.*l1-6.))/((20.*l1-20.)*l4-10.*l1+10.);
        %% Uniofrm Distribution
        TitleText = 'ABa/P2-25 ABa/EMS-25 ABp/P2-25 ABp/EMS-25';
        PossibleL1s(i,j) = l1;
        PossibleL2s(i,j) = -(.5*(-1.+(-2.*l1+2.)*l4))/l1;
        PossibleL3s(i,j) = -(1.*(-1.+(-2.*l1+2.)*l4))/(2.-(2.*(-2.*l1+2.))*l4);
        PossibleL4s(i,j) = l4;
        PossibleL5s(i,j) = .5;
        %% my Distribution
%         TitleText = "My Data";
%         PossibleL1s(i,j) = l1;
%         PossibleL2s(i,j) = -(0.2000000000e-1*(-29.+(-50.*l1+50.)*l4))/l1;
%         PossibleL3s(i,j) = -(1.*(-8.-(1.*(-25.*l1+25.))*l4*((-650.*l1+650.)*l4+400.*l1-377.)/((1250.*l1-1250.)*l4-725.*l1+725.)))/(14.50000000-(.5000000000*(-50.*l1+50.))*l4);
%         PossibleL4s(i,j) = l4;
%         PossibleL5s(i,j) = -(1.*((-650.*l1+650.)*l4+400.*l1-377.))/((1250.*l1-1250.)*l4-725.*l1+725.);
        i=i+1;
    end
    j=j+1;
end
DidntPass = PossibleL1s>1 | PossibleL1s<0 | PossibleL2s>1 | PossibleL2s<0 | PossibleL3s>1 | PossibleL3s<0 | PossibleL4s>1 | PossibleL4s<0 | PossibleL5s>1 | PossibleL5s<0;
PossibleL1s(DidntPass) = NaN;
PossibleL2s(DidntPass) = NaN;
PossibleL3s(DidntPass) = NaN;
PossibleL4s(DidntPass) = NaN;
PossibleL5s(DidntPass) = NaN;

% figure
% subplot(3,2,1)
% surf(unknownl1s,unknownl2s,PossibleL2s,'linestyle','none')
% xlabel(XTitle)
% ylabel(YTitle)
% zlabel('l2')
% view(2)
% subplot(3,2,2)
% surf(unknownl1s,unknownl2s,PossibleL3s,'linestyle','none')
% xlabel(XTitle)
% ylabel(YTitle)
% zlabel('l3')
% view(2)
% subplot(3,2,3)
% surf(unknownl1s,unknownl2s,PossibleL5s,'linestyle','none')
% xlabel(XTitle)
% ylabel(YTitle)
% zlabel('l5')
% view(2)

figure
ProbaABam = PossibleL1s.*PossibleL2s;
ProbaABad = (1.-PossibleL1s).*PossibleL4s;
ProbABpm = PossibleL1s.*(1.-PossibleL2s);
ProbABpd = (1.-PossibleL1s).*(1.-PossibleL4s);
ProbaEMSm = (1.-PossibleL1s).*(1.-PossibleL5s);
ProbaEMSd = PossibleL1s.*(1.-PossibleL3s);
ProbaP2m = (1.-PossibleL1s).*PossibleL5s;
ProbaP2d = PossibleL1s.*PossibleL3s;

subplot(2,2,1)
surf(unknownl1s,unknownl2s,ProbaABam-ProbaABad,'linestyle','none')
shading interp
xlabel(XTitle)
ylabel(YTitle)
title('ABa')
xlim([0 1])
ylim([0 1])
colorbar()
view(2)

subplot(2,2,2)
surf(unknownl1s,unknownl2s,ProbABpm-ProbABpd,'linestyle','none')
xlabel(XTitle)
ylabel(YTitle)
title('ABp')
xlim([0 1])
ylim([0 1])
colorbar()
view(2)

subplot(2,2,3)
surf(unknownl1s,unknownl2s,ProbaEMSm-ProbaEMSd,'linestyle','none')
xlabel(XTitle)
ylabel(YTitle)
title('EMS')
xlim([0 1])
ylim([0 1])
colorbar()
view(2)

subplot(2,2,4)
surf(unknownl1s,unknownl2s,ProbaP2m-ProbaP2d,'linestyle','none')
xlabel(XTitle)
ylabel(YTitle)
title('P_2')
xlim([0 1])
ylim([0 1])
colorbar()
view(2)

sgtitle([TitleText, " Probability w/ M_0 - Probability w/ D_0 "])

figure
Proba1m = PossibleL1s.*PossibleL2s.*PossibleL3s;
Proba1d = (1.-PossibleL1s).*PossibleL4s.*PossibleL5s;
Proba2m = PossibleL1s.*PossibleL2s.*(1.-PossibleL3s);
Proba2d = (1.-PossibleL1s).*PossibleL4s.*(1.-PossibleL5s);
Proba3m = PossibleL1s.*(1.-PossibleL2s).*PossibleL3s;
Proba3d = (1.-PossibleL1s).*(1.-PossibleL4s).*PossibleL5s;
Proba4m = PossibleL1s.*(1.-PossibleL2s).*(1.-PossibleL3s);
Proba4d = (1.-PossibleL1s).*(1.-PossibleL4s).*(1.-PossibleL5s);

subplot(2,2,1)
surf(unknownl1s,unknownl2s,Proba1m+Proba1d,'linestyle','none')
xlabel(XTitle)
ylabel(YTitle)
zlabel('P1')
xlim([0 1])
ylim([0 1])

subplot(2,2,2)
surf(unknownl1s,unknownl2s,Proba2m+Proba2d,'linestyle','none')
xlabel(XTitle)
ylabel(YTitle)
zlabel('P2')
xlim([0 1])
ylim([0 1])

subplot(2,2,3)
surf(unknownl1s,unknownl2s,Proba3m+Proba3d,'linestyle','none')
xlabel(XTitle)
ylabel(YTitle)
zlabel('P3')
xlim([0 1])
ylim([0 1])

subplot(2,2,4)
surf(unknownl1s,unknownl2s,Proba4m+Proba4d,'linestyle','none')
xlabel(XTitle)
ylabel(YTitle)
zlabel('P4')
xlim([0 1])
ylim([0 1])


figure
subplot(2,3,1)
surf(unknownl1s,unknownl2s,PossibleL1s,'linestyle','none')
title(XTitle)
xlabel(XTitle)
ylabel(YTitle)
xlim([0 1])
ylim([0 1])
zlim([0 1])
view(2)
subplot(2,3,2)
surf(unknownl1s,unknownl2s,PossibleL2s,'linestyle','none')
title('l2')
xlabel(XTitle)
ylabel(YTitle)
xlim([0 1])
ylim([0 1])
zlim([0 1])
view(2)
subplot(2,3,3)
surf(unknownl1s,unknownl2s,PossibleL3s,'linestyle','none')
title('l3')
xlabel(XTitle)
ylabel(YTitle)
xlim([0 1])
ylim([0 1])
zlim([0 1])
view(2)
subplot(2,3,4)
surf(unknownl1s,unknownl2s,PossibleL4s,'linestyle','none')
title(YTitle)
xlabel(XTitle)
ylabel(YTitle)
xlim([0 1])
ylim([0 1])
zlim([0 1])
view(2)
subplot(2,3,5)
surf(unknownl1s,unknownl2s,PossibleL5s,'linestyle','none')
title('l5')
xlabel(XTitle)
ylabel(YTitle)
xlim([0 1])
ylim([0 1])
zlim([0 1])
view(2)
subplot(2,3,6)
colorbar()
sgtitle(TitleText)

figure
subplot(2,3,1)
histogram(PossibleL1s(:),21,'BinLimits',[0,1],'Normalization','probability')
title('l1')

subplot(2,3,2)
histogram(PossibleL2s(:),21,'BinLimits',[0,1],'Normalization','probability')
title('l2')

subplot(2,3,3)
histogram(PossibleL3s(:),21,'BinLimits',[0,1],'Normalization','probability')
title('l3')

subplot(2,3,4)
histogram(PossibleL4s(:),21,'BinLimits',[0,1],'Normalization','probability')
title('l4')

subplot(2,3,5)
histogram(PossibleL5s(:),21,'BinLimits',[0,1],'Normalization','probability')
title('l5')
sgtitle(TitleText)