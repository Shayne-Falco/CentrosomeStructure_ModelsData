close all
clear all

%% Load all MMH sets ran

running=1;
ToPlot=0;
mmhStyle=1;
deltaT=2; %must be a divisor of 10: 10,5,2,1
recoveryTime = 50/(deltaT)+1;

LeadBleach=1;
BleachTime=-170;

S=load('AverageLeadRecovery.mat');
AverageLeadRecovery = S.AverageLeadRecovery(1:deltaT:end);

S=load('AverageLaggRecovery.mat');
AverageLaggRecovery = S.AverageLaggRecovery(1:deltaT:end);

S=load('PhotobleachCorrection.mat');
DataDisplay = S.CorrectforPhotobleach;

t = DataDisplay(:,1);
TimeFrame = t>-331 & t<1;
DataIntT = DataDisplay(TimeFrame,1);
Cytos = DataDisplay(TimeFrame,10:5:end);
Leads = DataDisplay(TimeFrame,8:5:end);
Laggs = DataDisplay(TimeFrame,9:5:end);
Totls = DataDisplay(TimeFrame,11:5:end);
for i = 1:size(Totls,2)
    AvgTot = mean(Totls(:,i));
    Totls(:,i) = Totls(:,i)./AvgTot;
    Cytos(:,i) = Cytos(:,i)./AvgTot;
    Leads(:,i) = Leads(:,i)./AvgTot;
    Laggs(:,i) = Laggs(:,i)./AvgTot;
end
AvgTotl = mean(Totls,2);
DataIntCyto = mean(Cytos,2);
DataIntLead = mean(Leads,2);
DataIntLagg = mean(Laggs,2);
lengthDataInt=length(DataIntLagg);
titlesCon={'a1','a2','a3','a4','a5','a6','a7','n1','n2','t','Lead Error','Lagger Error','Cytoplasm Error','Lead Frap Error','Lagg Frap Error'}; %Constant and Conserved

DataInt = [DataIntLead DataIntLagg DataIntCyto];
t_max=lengthDataInt;

S = load(append('FRAPFitToBoth_mmh',num2str(mmhStyle),'_Both.mat'));
AllParameters=S.AnalysisParameters;
count = size(AllParameters,1);

avg_t_Half=[];
std_t_half=[];
avg_Vmax=[];
avg_Koff=[];
std_Vmax=[];
std_Koff=[];
avg_RecoveryCurve=[];

t_half=[];
PassedAnalysis = [];
Vmax=[];
Koff=[];
DistanceMetric = [];
RecoveryCurves=zeros(count,-BleachTime/deltaT+1);
Recruitment = zeros(count,t_max,6);
FinalRecovery = [BleachTime:deltaT:BleachTime+recoveryTime*deltaT]';
FinalLeadRecruitment = [-330:deltaT:0]';
FinalLaggRecruitment = [-330:deltaT:0]';
FinalCytoRecruitment = [-330:deltaT:0]';
FinalTotlRecruitment = [-330:deltaT:0]';
for ParameterCounter=1:count
    FRAPcurveF=[];

    %% Parameters

    B=0;
    F=1-B;
    ConstantCy=0; %0 for conservation, else that is Cy value
    %% setup the model

    C1f = DataIntLead(1); %leading MTOC initial value
    C2f = DataIntLagg(2); %Lagging MTOC initial value
    Cyf = DataIntCyto(1);
    C_0 = C1f+C2f+Cyf; %

    C1b=0; %leading MTOC initial value
    C2b=0; %Lagging MTOC initial value
    Cyb=0;
    hits=AllParameters(ParameterCounter,:);
    array=hits(1:9);
    Parameters=[array,C_0,F,B];
    %Parameters=[a_1,a_2,a_3,a_4,a_5,a_6,a_7,n1,n2,C_0,F,B];

    %% ODE solving

    y0=[C1b,C1f,C2b,C2f,Cyb,Cyf];
    y=zeros(t_max,6);
    %if ConstantCy

    y(1,:)=y0;
    %end
    running=1;

    timer_FRAP=0;
    t=1;
    time=0;
    recoverytime=NaN;


    while running
        [~,y1] = ode15s(@(t,y0)ODE_DL_FRAP_Dynamic(t,Parameters,y0),[t-.1*deltaT,t],y0);
        y1=y1(end,:);
        t=t+1;
        time=time+.1*deltaT;
        if time>=t_max-1.1
            running=0;
        end
        y0=y1;

        if DataIntT(1)+time*10>=BleachTime && (B==0)
            if LeadBleach
                B=y0(2);
            else
                B=y0(4);
            end

            F=1-B;
            Parameters(11)=F;
            Parameters(12)=B;
            if LeadBleach
                y0=[1,0,y0(3:end)/F]; %leading FRAP
            else
                y0=[y0(1:2)/F,1,0,y0(5:6)/F];
            end
        end

        y(t,:)=y0.*[B,F,B,F,B,F];
        if B>0
            if LeadBleach
                FRAPcurveF=[FRAPcurveF,y0(2)/B];
            else
                FRAPcurveF=[FRAPcurveF,y0(4)/B];
            end
        end

    end
    Recruitment(ParameterCounter,1:t,1:6)=y;
    
    subplot(2,2,1)
    plot(-330:deltaT:0,y(:,1)+y(:,2),'black','LineWidth',1)
    hold on

    subplot(2,2,2)
    plot(-330:deltaT:0,y(:,3)+y(:,4),'black','LineWidth',1)
    hold on


    subplot(2,2,3)
    plot(-330:deltaT:0,y(:,5)+y(:,6),'black','LineWidth',1)
    hold on

    subplot(2,2,4)
    plot(-330:deltaT:0,sum(y(:,1:6)'),'black','LineWidth',1)
    hold on

    FinalRecovery = [FinalRecovery FRAPcurveF(1:recoveryTime+1)'];
    FinalLeadRecruitment = [FinalLeadRecruitment y(:,1)+y(:,2)];
    FinalLaggRecruitment = [FinalLaggRecruitment y(:,3)+y(:,4)];
    FinalCytoRecruitment = [FinalCytoRecruitment y(:,5)+y(:,6)];
    FinalTotlRecruitment = [FinalTotlRecruitment sum(y(:,1:6)')'];

end

filename = append('Paper_Data/Curves/FRAP_mmh',num2str(mmhStyle),'/LeadRecruit.dat');
writematrix(FinalLeadRecruitment,string(filename),'filetype', 'text','Delimiter','tab')
filename = append('Paper_Data/Curves/FRAP_mmh',num2str(mmhStyle),'/LaggRecruit.dat');
writematrix(FinalLeadRecruitment,string(filename),'filetype', 'text','Delimiter','tab')
filename = append('Paper_Data/Curves/FRAP_mmh',num2str(mmhStyle),'/CytoRecruit.dat');
writematrix(FinalLeadRecruitment,string(filename),'filetype', 'text','Delimiter','tab')
filename = append('Paper_Data/Curves/FRAP_mmh',num2str(mmhStyle),'/TotlRecruit.dat');
writematrix(FinalLeadRecruitment,string(filename),'filetype', 'text','Delimiter','tab')
if LeadBleach
    filename = append('Paper_Data/Curves/FRAP_mmh',num2str(mmhStyle),'/LeadRecvovery.dat');
else
    filename = append('Paper_Data/Curves/FRAP_mmh',num2str(mmhStyle),'/LaggRecvovery.dat');
end
writematrix(FinalLeadRecruitment,string(filename),'filetype', 'text','Delimiter','tab')
subplot(2,2,1)
plot(DataIntT,DataIntLead,'r-','LineWidth',4)
%plot(DataIntT,y(:,1)+y(:,2),'black','LineWidth',1)
hold on
ylim([0 .2])
xlim([-330 0])
%title(sprintf('Leading; %0.4f%',DistanceSimData(1)))
title('Leading Compartment/Centrosome ')
%legend('Normalized AIR-1::GFP signal','Simulation')

subplot(2,2,2)

plot(DataIntT,DataIntLagg,'r','LineWidth',4)
%plot(DataIntT,y(:,3)+y(:,4),'black','LineWidth',1)
hold on
ylim([0 .2])
xlim([-330 0])
title('Lagging Compartment/Centrosome ')
%title(sprintf('Lagging; %0.4f%',DistanceSimData(2)))

subplot(2,2,3)

plot(DataIntT,DataIntCyto,'r','LineWidth',4)
%plot(DataIntT,y(:,5)+y(:,6),'black','LineWidth',1)
hold on
ylim([.75 .95])
xlim([-330 0])
title('Cytoplasm')
xlabel('time (sec) from NEBD')
%title(sprintf('Cyto; %0.4f%',DistanceSimData(3)))

subplot(2,2,4)

plot(DataIntT,AvgTotl,'r','LineWidth',4)
%plot(DataIntT,y(:,5)+y(:,6),'black','LineWidth',1)
hold on
ylim([.95 1.15])
xlim([-330 0])
title('Total')
xlabel('time (sec) from NEBD')
%title(sprintf('Cyto; %0.4f%',DistanceSimData(3)))
