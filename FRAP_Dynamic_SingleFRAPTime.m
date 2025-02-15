%%

%close all
%clear all
clearvars -except PassedAnalysis
%
%

%% Settings to change
photobleachModel = true;
LeadBleach=1;
StylesToPlot=1; %which mmh target density, 1=exp, 2=uniform, 3=none
UsedPassedAnalysis = 1; 
deltaT=1; %must be a divisor of 10: 10,5,2,1
%There is a vector "PassedAnalysis" That saves the successful parameter sets each run. 
%   Run the code first with this false then with true to only plot the
%   successful partameter sets.
%Keep this set to true when switching bleached compartments to get sets
%   that pass both bleaching compartmets.



%% Load all MMH sets ran
ToPlotDy = false;
count=1;
running=1;
%ToPlot=1;
recoveryTime = 50/(deltaT)+1;
Thalf_cutoff = 35;
Vmax_cutoff = 135;
BleachTime=-170;
LeadErrorCutoff=.3;
LaggErrorCutoff=.3;
CytoErrorCutoff=.35;


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

count=0;
FileExists=[];
passed=[];
FilesTried = 0;
for run=1:10000
    for Style = StylesToPlot
        if photobleachModel
            myfilename=['MMH_FRAP_Dynamic_Full_V2_style_',num2str(Style),'/Run_',num2str(run), '.mat'];
        else
            myfilename=['MMH_Dynamic_style_',num2str(Style),'/Run_',num2str(run), '.mat'];
        end
        if isfile(myfilename)
            FilesTried = FilesTried + 1;
            load(myfilename)
            if photobleachModel
                Leaderrors=hits(end,end-4);
                Laggerrors=hits(end,end-3);
                CytoErrors=hits(end,end-2);
            else
                Leaderrors=hits(end,end-2);
                Laggerrors=hits(end,end-1);
                CytoErrors=hits(end,end);
            end
            if UsedPassedAnalysis
                if Leaderrors<LeadErrorCutoff && Laggerrors<LaggErrorCutoff && CytoErrors < CytoErrorCutoff && ismember(FilesTried,PassedAnalysis)
                    if ToPlotDy
                        f = figure(1);
                        f.Visible = 'off';
                        for i = 1:15
                            subplot(5,3,i)
                            plot(1:length(hits(30:end,i)),hits(30:end,i))
                            title(titlesCon(i))
                        end
                        saveas(f,append('DyFRAPFit_Lagg/',string(run)),'png');
                        close(f)
                    end
                    count=count+1;
                    results{count}=hits;
                    FileExists=[FileExists run];
                    passed=[passed FilesTried];
                end
            
            elseif Leaderrors<LeadErrorCutoff && Laggerrors<LaggErrorCutoff && CytoErrors < CytoErrorCutoff
                if ToPlotDy
                    f = figure(1);
                    f.Visible = 'off';
                    for i = 1:15
                        subplot(5,3,i)
                        plot(1:length(hits(30:end,i)),hits(30:end,i))
                        title(titlesCon(i))
                    end
                    saveas(f,append('DyFRAPFit_Lagg/',string(run)),'png');
                    close(f)
                end
                count=count+1;
                results{count}=hits;
                FileExists=[FileExists run];
                passed=[passed FilesTried];
            end
        end
    end
end
FileExists


%AllParameters=[];
%FinalParameters=zeros(count,10);
% passed=[];

% for resu=1:count
%     hits=results{resu};
%     passed=[passed resu];
%
% end

avg_t_Half=[];
std_t_half=[];
avg_Vmax=[];
avg_Koff=[];
std_Vmax=[];
std_Koff=[];
avg_RecoveryCurve=[];



t_half=[];
PassedAnalysis = [];
AnalysisParameters = [];
Dups = [];
Vmax=[];
Koff=[];
DistanceMetric = [];
RecoveryCurves=zeros(count,-BleachTime/deltaT+1);
for ParameterCounter=1:count
    %resu=passed(ParameterCounter);
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
    hits=results{ParameterCounter};
    array=hits(end,1:9);
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
        [~,y1] = ode15s(@(t,y0)ODE_DL_FRAP_Dynamic(t,Parameters,y0),[time-.1*deltaT,time],y0);
        y1=y1(end,:);
        t=t+1;
        time=time+.1*deltaT;
        if time>=t_max-1
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
%     DistanceSimData = SimulationDataDistanceNormalizedFRAP(DataInt,[y(1:t,1)+y(1:t,2),y(1:t,3)+y(1:t,4),y(1:t,5)+y(1:t,6)],FRAPcurveF,recoveryTime,LeadBleach,AverageLeadRecovery,AverageLaggRecovery)
    DistanceSimData = sum(SimulationDataDistanceNormalized(DataInt,[y(1:10/deltaT:t,1)+y(1:10/deltaT:t,2),y(1:10/deltaT:t,3)+y(1:10/deltaT:t,4),y(1:10/deltaT:t,5)+y(1:10/deltaT:t,6)]));
    
    ft = fittype( 'a-a*exp(-x/b)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Robust = 'LAR';
    
    if recoveryTime
        recoverytime=(find(FRAPcurveF(1:recoveryTime+1)-FRAPcurveF(recoveryTime+1)/2>0,1,'first')-1);
        [xData, yData] = prepareCurveData(  0:deltaT:recoveryTime*deltaT,100*FRAPcurveF(1:recoveryTime+1) );
        opts.StartPoint = [FRAPcurveF(recoveryTime+1)*100 -log(0.5)*recoverytime];
    else
        recoverytime=(find(FRAPcurveF-FRAPcurveF(end)/2>0,1,'first')-1)*10;
        [xData, yData] = prepareCurveData(  0:10:10*length(FRAPcurveF)-1,100*FRAPcurveF );
        opts.StartPoint = [FRAPcurveF(end)*100 -log(0.5)*recoverytime];
    end
    
    
    %% Fit: 'Vmax-Vmax exp(-t/t_off)'.
    coeffvals=[NaN,NaN];
    % Fit model to data.
    
    try
        [fitresult, gof] = fit( xData, yData, ft, opts );
        
        coeffvals = coeffvalues(fitresult);
    end
    
    %             recoverytime=(find(FRAPcurveF-FRAPcurveF(end)/2>0,1,'first')-1)*10;
    if coeffvals(1)<Vmax_cutoff && coeffvals(2) < Thalf_cutoff && coeffvals(2)>10
        if ~ismember(array(1),AnalysisParameters)
            t_half=[t_half -log(0.5)*coeffvals(2)];
            Vmax=[Vmax,coeffvals(1)/100];
            Koff=[Koff,coeffvals(2)];
            RecoveryCurves(ParameterCounter,1:length(FRAPcurveF))=FRAPcurveF;
            PassedAnalysis = [PassedAnalysis passed(ParameterCounter)];
            AnalysisParameters(ParameterCounter,:)=[array,-log(0.5)*coeffvals(2),coeffvals(1)/100,1,FRAPcurveF];
        else
            Dups = [Dups passed(ParameterCounter)];
        end
    else
        t_half=[t_half NaN];
        Vmax=[Vmax,NaN];
        Koff=[Koff,NaN];
        RecoveryCurves(ParameterCounter,1:length(FRAPcurveF))=FRAPcurveF;
        AnalysisParameters(ParameterCounter,:)=[array,-log(0.5)*coeffvals(2),coeffvals(1)/100,0,FRAPcurveF];
     end
    
    
    
    
end
PassedTests = size(t_half,2);
%avg_RecoveryCurve=[avg_RecoveryCurve;mean(RecoveryCurves(isfinite(RecoveryCurves(:,1))&RecoveryCurves(:,2)>0,:))];
avg_RecoveryCurve=[avg_RecoveryCurve;RecoveryCurves(isfinite(RecoveryCurves(:,1))&real(RecoveryCurves(:,2))>0,:)];
avg_t_Half=[avg_t_Half;BleachTime, mean(t_half,'omitnan')];
std_t_half=[std_t_half,std(t_half,'omitnan')];
avg_Vmax=[avg_Vmax;BleachTime,mean(Vmax,'omitnan')];
avg_Koff=[avg_Koff;BleachTime,mean(Koff,'omitnan')];
std_Vmax=[std_Vmax,std(Vmax,'omitnan')];
std_Koff=[std_Koff,std(Koff,'omitnan')];

%
figure
subplot(2,2,1)
%errorbar(avg_t_Half(:,1),avg_t_Half(:,2),std_t_half)
%plot(ones(PassedTests,1)+(rand(PassedTests,1)-.5)/10,t_half,'*')
histogram(t_half,'NumBins',10)
% xlabel('Time of Photobleaching')
% ylabel('Half Recovery Time')
xlabel('Half Recovery Time')
% xlim([BleachTime-10 BleachTime+10])
if LeadBleach
    title('half recovery time of leading compartment')
else
    title('half recovery time of lagging compartment')
end
subplot(2,2,2)
% errorbar(avg_Vmax(:,1),avg_Vmax(:,2),std_Vmax)
histogram(Vmax,'NumBins',10)
% xlabel('Time of Photobleaching')
% ylabel('Vmax from Recovery Curve Fit')
xlabel('Vmax from Recovery Curve Fit')
% xlim([BleachTime-10 BleachTime+10])
if LeadBleach
    title('Vmax of leading compartment')
else
    title('Vmax of lagging compartment')
end
subplot(2,2,3)
% errorbar(avg_Koff(:,1),avg_Koff(:,2),std_Koff)
histogram(Koff,'NumBins',10)
% xlabel('Time of Photobleaching')
% ylabel('K_{off} from Recovery Curve Fit')
xlabel('K_{off} from Recovery Curve Fit')
% xlim([BleachTime-10 BleachTime+10])
if LeadBleach
    title('K_{off} of leading compartment')
else
    title('K_{off} of lagging compartment')
end
subplot(2,2,4)
plot(BleachTime:deltaT:-1,RecoveryCurves(:,1:length(FRAPcurveF)))
xlim([BleachTime-10 BleachTime+recoveryTime*deltaT])
% legend(string(BleachTime))
xlabel('Time post Photobleaching')
title('Average Recovery Curve')
%ylabel('K_{off} from Recovery Curve Fit')
%xlim([BleachTimes(1)-10 BleachTimes(end)+10])
%     if LeadBleach
%         title('K_{off} of leading compartment')
%     else
%         title('K_{off} of lagging compartment')
%     end
%     if LeadBleach
%         savefig('Lead Bleach_to0')
%     else
%         savefig('Lag Bleach_to0')
%     end
PassedAnalysis


maximumsPlot = [1,.25,10,10,10,.12,.12,3,3,60,2];
% if ToPlot
%     titles={'a1','a2','a3','a4','a5','a6','a7','n1','n2','Thalf','Vmax'};
%     for i=1:10
%         for j=i+1:11
%             f = figure('visible','off');
%             colormap flag
%             x=AnalysisParameters(:,i);
%             y=AnalysisParameters(:,j);
%             scatter(x,y,[],AnalysisParameters(:,12))
%             xlabel(titles(i))
%             ylabel(titles(j))
%             if LeadBleach
%                 title('Black = Leading recovery similar to data')
%             else
%                 title('Black = Lagging recovery similar to data')
%             end
%             if i==8 && j==9
%                 hold on
%                 plot(0:.1:3,0:.1:3)
%                 legend('Data',sprintf('n_1>n_2; %0.2f%%', 100*sum(x>y)/length(y)),'Location','best');
%             end
%             if i == 8 || i == 9
%                 xlim([0 3])
%             else
%                 xlim([0 min(maximumsPlot(i),max(x))+.1])
%             end
%             if j == 8 || j == 9
%                 ylim([0 3])
%             else
%                 ylim([0 min(maximumsPlot(j),max(y))+.1])
%             end
%             if LeadBleach
%                 saveas(f,append('FRAPcorFRAPFit/',string(titles(i)),string(titles(j)),'Lead'),'png');
%             else
%                 saveas(f,append('FRAPcorFRAPFit/',string(titles(i)),string(titles(j)),'Lagg'),'png');
%             end
%             
%         end
%     end
% end

S=load('AverageLaggRecovery.mat');
AverageLaggRecovery = S.AverageLaggRecovery;

S=load('AverageLeadRecovery.mat');
AverageLeadRecovery = S.AverageLeadRecovery;

figure
for ParameterCounter=1:count
    if AnalysisParameters(ParameterCounter,12)==1
        plot(BleachTime:deltaT:-1,AnalysisParameters(ParameterCounter,13:end),'black')
    else
        plot(BleachTime:deltaT:-1,AnalysisParameters(ParameterCounter,13:end),'r')
    end
    hold on
    xlim([BleachTime-10 BleachTime+recoveryTime*deltaT])
    
end
if LeadBleach
    plot(BleachTime:deltaT:BleachTime+50,AverageLeadRecovery(1:deltaT:end),'r','LineWidth', 2)
    title('Leading Frap, Black recovery similar to data')
else
    plot(BleachTime:deltaT:BleachTime+50,AverageLaggRecovery(1:deltaT:end),'r','LineWidth', 2)
    title('Lagging Frap, Black recovery similar to data')
end