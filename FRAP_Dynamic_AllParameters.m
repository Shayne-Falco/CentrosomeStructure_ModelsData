%%
close all
clear all
%


%% Load all MMH sets ran
count=1;
running=1;
ToPlot=1;
recoveryTime=6;
RunParameters=1000;
BleachTimes=-200:10:-140;
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
LeadErrorCutoff=.35;
LaggErrorCutoff=.35;
CytoErrorCutoff=.45;
count=0;
FileExists=[];
passed=[];
FilesTried = 0;
for run=1:10000
    for Style = 1
        myfilename=['MMH_FRAP_Dynamic_Full_V2_style_',num2str(Style),'/Run_',num2str(run), '.mat'];
        if isfile(myfilename)
            FilesTried = FilesTried + 1;
            load(myfilename)
            Leaderrors=hits(end,end-2);
            Laggerrors=hits(end,end-1);
            CytoErrors=hits(end,end);
            
            if Leaderrors<LeadErrorCutoff && Laggerrors<LaggErrorCutoff && CytoErrors < CytoErrorCutoff %&& ismember(FilesTried,PassedAnalysis)
                count=count+1;
                results{count}=hits;
                FileExists=[FileExists run];
                passed=[passed FilesTried];
            end
        end
    end
end
FileExists


AllParameters=[];
FinalParameters=zeros(count,10);
% passed=[];

% for resu=1:count
%     hits=results{resu};
%     passed=[passed resu];
%
% end

for LeadBleach=[1,0]
    avg_t_Half=[];
    std_t_half=[];
    avg_Vmax=[];
    avg_Koff=[];
    std_Vmax=[];
    std_Koff=[];
    avg_RecoveryCurve=[];
    all_t_half=zeros(count,length(BleachTimes));
    all_Vmax=zeros(count,length(BleachTimes));
    all_Koff=zeros(count,length(BleachTimes));
    all_DistanceMetric=zeros(count,length(BleachTimes));
    for BleachCounter=1:length(BleachTimes)
        BleachTime=BleachTimes(BleachCounter);
        
        t_half=[];
        Vmax=[];
        Koff=[];
        DistanceMetric = [];
        RecoveryCurves=zeros(count,t_max);
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
            recoverytime=NaN;
            
            
            while running
                [~,y1] = ode15s(@(t,y0)ODE_DL_FRAP_Dynamic(t,Parameters,y0),[t-1,t],y0);
                y1=y1(end,:);
                if ~isreal(y1)
                    %disp('Not Real')
                    disp(FileExists(resu))
                end
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
            
            DistanceSimData = sum(SimulationDataDistanceNormalized(DataInt,[y(1:t,1)+y(1:t,2),y(1:t,3)+y(1:t,4),y(1:t,5)+y(1:t,6)]));
            
            ft = fittype( 'a-a*exp(-x/b)', 'independent', 'x', 'dependent', 'y' );
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.Display = 'Off';
            opts.Robust = 'LAR';
            
            if recoveryTime
                recoverytime=(find(FRAPcurveF(1:recoveryTime+1)-FRAPcurveF(recoveryTime+1)/2>0,1,'first')-1)*10;
                [xData, yData] = prepareCurveData(  0:10:recoveryTime*10,100*FRAPcurveF(1:recoveryTime+1) );
                opts.StartPoint = [FRAPcurveF(recoveryTime+1) -log(0.5)*recoverytime];
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
            if isfinite(coeffvals(1)) && coeffvals(1)<200
                t_half=[t_half -log(0.5)*coeffvals(2)];
                Vmax=[Vmax,coeffvals(1)/100];
                Koff=[Koff,coeffvals(2)];
                DistanceMetric = [DistanceMetric,DistanceSimData];
                RecoveryCurves(ParameterCounter,1:length(FRAPcurveF))=FRAPcurveF;
            else
                t_half=[t_half NaN];
                Vmax=[Vmax,NaN];
                Koff=[Koff,NaN];
                RecoveryCurves(ParameterCounter,1)=NaN;
            end
            
            
            
            
        end
        %avg_RecoveryCurve=[avg_RecoveryCurve;mean(RecoveryCurves(isfinite(RecoveryCurves(:,1))&RecoveryCurves(:,2)>0,:))];
        avg_RecoveryCurve=[avg_RecoveryCurve;RecoveryCurves(isfinite(RecoveryCurves(:,1))&real(RecoveryCurves(:,2))>0,:)];
        avg_t_Half=[avg_t_Half;BleachTime, mean(t_half,'omitnan')];
        std_t_half=[std_t_half,std(t_half,'omitnan')];
        avg_Vmax=[avg_Vmax;BleachTime,mean(Vmax,'omitnan')];
        avg_Koff=[avg_Koff;BleachTime,mean(Koff,'omitnan')];
        std_Vmax=[std_Vmax,std(Vmax,'omitnan')];
        std_Koff=[std_Koff,std(Koff,'omitnan')];
        all_t_half(1:length(t_half),BleachCounter)=t_half';
        all_Vmax(1:length(Vmax),BleachCounter)=Vmax';
        all_Koff(1:length(Koff),BleachCounter)=Koff';
        all_DistanceMetric(1:length(DistanceMetric),BleachCounter)=DistanceMetric';
    end
    figure
    subplot(2,2,1)
    errorbar(avg_t_Half(:,1),avg_t_Half(:,2),std_t_half)
    xlabel('Time of Photobleaching')
    ylabel('Half Recovery Time')
    xlim([BleachTimes(1)-10 BleachTimes(end)+10])
    if LeadBleach
        title('half recovery time of leading compartment')
    else
        title('half recovery time of lagging compartment')
    end
    subplot(2,2,2)
    errorbar(avg_Vmax(:,1),avg_Vmax(:,2),std_Vmax)
    xlabel('Time of Photobleaching')
    ylabel('Vmax from Recovery Curve Fit')
    xlim([BleachTimes(1)-10 BleachTimes(end)+10])
    if LeadBleach
        title('Vmax of leading compartment')
    else
        title('Vmax of lagging compartment')
    end
    subplot(2,2,3)
    errorbar(avg_Koff(:,1),avg_Koff(:,2),std_Koff)
    xlabel('Time of Photobleaching')
    ylabel('K_{off} from Recovery Curve Fit')
    xlim([BleachTimes(1)-10 BleachTimes(end)+10])
    if LeadBleach
        title('K_{off} of leading compartment')
    else
        title('K_{off} of lagging compartment')
    end
    subplot(2,2,4)
    C=winter(length(BleachTimes));
    for i=1:length(BleachTimes)
        if recoveryTime
            plot(BleachTimes(i):10:BleachTimes(i)+recoveryTime*10,avg_RecoveryCurve(i,1:recoveryTime+1),'color',C(i,:))
        else
            plot(BleachTimes(i):10:0,avg_RecoveryCurve(i,1:length(BleachTimes(i):10:0)),'color',C(i,:))
        end
        hold on
    end
    
    xlim([-300 10])
    legend(string(BleachTimes))
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
    all_DistanceMetric
end