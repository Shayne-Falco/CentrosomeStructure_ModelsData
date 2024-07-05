%%
close all
%clear all
%
%

%% Load all MMH sets ran
count=1;
running=1;
ToPlot=0;
BleachTime=-250;
recoveryTime = 5;
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

DataInt = [DataIntLead DataIntLagg DataIntCyto];
t_max=lengthDataInt;


LeadBleach=0; % 0 = Lagg, 1 == Lead is bleached


t_half=NaN;
Vmax=NaN;
Koff=NaN;
FRAPcurveF=[];

%% Parameters
B=0;
F=1-B;
%% setup the model

C1f = DataIntLead(1); %leading MTOC initial value
C2f = DataIntLagg(2); %Lagging MTOC initial value
Cyf = DataIntCyto(1);
C_0 = C1f+C2f+Cyf; %

C1b=0; %leading MTOC initial value
C2b=0; %Lagging MTOC initial value
Cyb=0;

Parameters=[array,C_0,F,B];
%Parameters=[a_1,a_2,a_3,a_4,a_5,a_6,a_7,n1,n2,C_0,F,B];

%% ODE solving
y0=[C1b,C1f,C2b,C2f,Cyb,Cyf];
y=zeros(t_max,6);

y(1,:)=y0;
running=1;
timer_FRAP=0;
t=1;
recoverytime=NaN;


while running
    [~,y1] = ode15s(@(t,y0)ODE_DL_FRAP_Dynamic(t,Parameters,y0),[t-1,t],y0);
    y1=y1(end,:);
    t=t+1;
    if t==t_max
        running=0;
    end
    y0=y1;
    if DataIntT(t)==BleachTime
        
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
FRAPcurveF
DistanceSimData = sum(SimulationDataDistanceNormalized(DataInt,[y(1:t,1)+y(1:t,2),y(1:t,3)+y(1:t,4),y(1:t,5)+y(1:t,6)]));

ft = fittype( 'a-a*exp(-x/b)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Robust = 'LAR';

if recoveryTime
    recoverytime=(find(FRAPcurveF(1:recoveryTime+1)-FRAPcurveF(recoveryTime+1)/2>0,1,'first')-1)*10;
    [xData, yData] = prepareCurveData(  0:10:recoveryTime*10,FRAPcurveF(1:recoveryTime+1) );
    opts.StartPoint = [FRAPcurveF(recoveryTime+1) -log(0.5)*recoverytime];
else
    recoverytime=(find(FRAPcurveF-FRAPcurveF(end)/2>0,1,'first')-1)*10;
    [xData, yData] = prepareCurveData(  0:10:10*length(FRAPcurveF)-1,100*FRAPcurveF );
    opts.StartPoint = [FRAPcurveF(recoveryTime)*100 -log(0.5)*recoverytime];
end


%% Fit: 'Vmax-Vmax exp(-t/t_off)'.
coeffvals=[NaN,NaN];
% Fit model to data.

% try
    [fitresult, gof] = fit( xData, yData, ft, opts );
    
    coeffvals = coeffvalues(fitresult)
% end

%             recoverytime=(find(FRAPcurveF-FRAPcurveF(end)/2>0,1,'first')-1)*10;
if isfinite(coeffvals(1)) && coeffvals(1)<1000
    t_half=[t_half -log(0.5)*coeffvals(2)];
    Vmax=[Vmax,coeffvals(1)/100];
    Koff=[Koff,coeffvals(2)];

else
    t_half=[t_half NaN];
    Vmax=[Vmax,NaN];
    Koff=[Koff,NaN];

end
titles = {'C1B','C1F','C1','C2B','C2F','C2','CyB','CyF','Cy'};
j=1;
figure(1)
for i = 1:9
    subplot(3,3,i)
    if i == 3
        plot(DataIntT,y(:,1)+y(:,2),DataIntT,DataIntLead)
    elseif i == 6
        plot(DataIntT,y(:,3)+y(:,4),DataIntT,DataIntLagg)
    elseif i == 9
        plot(DataIntT,y(:,5)+y(:,6),DataIntT,DataIntCyto)
    else
        plot(DataIntT,y(:,j))
        j=j+1;
    end
    title(titles(i))
    xlim([DataIntT(1) DataIntT(end)])
end
figure(2)
plot(BleachTime:10:10*(length(FRAPcurveF)-1)+BleachTime,FRAPcurveF)

