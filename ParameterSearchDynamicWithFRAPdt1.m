clear all
Maxruns=1000000;
BleachTime=-170;

UseAverageLeadLagg = 1;
UseSamneParamsLeadLagg = 1;

LeadCutoff = .2;
LaggCutoff = .2;
CytoCutoff = .3;
FRAPCutoff = 10.2;
deltaT=2; %must be a divisor of 10: 10,5,2,1
recoveryTime = 50/(deltaT)+1;
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
DataIntTotl = mean(Totls,2);
lengthDataInt=length(DataIntLagg);

if UseAverageLeadLagg
    DataInt = [(DataIntLead+DataIntLagg)/2 (DataIntLead+DataIntLagg)/2 DataIntCyto];
else
    DataInt = [DataIntLead DataIntLagg  DataIntCyto];
end

%% Set up model
%LeadBleach=1; % 0 = Lagg, 1 == Lead is bleached
t_half=NaN;
Vmax=NaN;
Koff=NaN;


MMHProp=.25;
t_max=lengthDataInt;


C1f = DataIntLead(1); %leading MTOC initial value
C2f = DataIntLagg(2); %Lagging MTOC initial value
Cyf = DataIntCyto(1);
C_0 = C1f+C2f+Cyf; %

C1b=0; %leading MTOC initial value
C2b=0; %Lagging MTOC initial value
Cyb=0;


TotalSamples=1:50;
for mmhStyle = 3
    parfor(Samplings=TotalSamples)
    %for Samplings=1
        %     for Samplings=TotalSamples
        FRAPcurveFLead=[];
        FRAPcurveFLagg=[];
        BLead=0;
        FLead=1-BLead;
        BLagg=0;
        FLagg=1-BLagg;
        runningMMH=1;
        myfilename=['MMH_FRAP_Dynamic_dt1_SymmetricLeadLagg/Run_',num2str(Samplings), '.mat'];
        if isfile(myfilename) % if the file already exists it will not ovberwrite it. It will skip to the next file
            continue
        end
        disp('starting a fit')

        a_1 = .3+.3*rand(1); 	%
        a_2 = rand(1); 	    %
        a_3 = rand(1); 	%
        a_4 = rand(1); 	%
        a_5 = rand(1); 	%
        a_6 = 0.02+.01*rand(1); 	%
        a_7 = a_6; 	%
        n1 = 0.5+rand(1);   	%
        n2 = 0.5+rand(1);

        array=[a_1 a_2 a_3 a_4 a_5 a_6 a_7 n1 n2];
        count=1;
        hits=[];
        hits(1,:)=[array 1 inf inf inf inf inf];
        j=1;

        ParametersLead=[array C_0 FLead BLead];
        ParametersLagg=[array C_0 FLagg BLagg];

        y0Lead=[C1b,C1f,C2b,C2f,Cyb,Cyf];
        y0Lagg=y0Lead;

        ylead=zeros(t_max,6);
        ylead(1,:)=y0Lead;

        ylagg=zeros(t_max,6);
        ylagg(1,:)=y0Lead;

        runningODE=1;
        BadParameters = 0;
        t=1;
        time=0;
        while runningODE
            [~,y1Lead] = ode15s(@(t,y0Lead)ODE_DL_FRAP_Dynamic(t,ParametersLead,y0Lead),[time-0.1*deltaT,time],y0Lead);
            y1Lead=y1Lead(end,:);
            [~,y1Lagg] = ode15s(@(t,y0Lagg)ODE_DL_FRAP_Dynamic(t,ParametersLagg,y0Lagg),[time-0.1*deltaT,time],y0Lagg);
            y1Lagg=y1Lagg(end,:);
            t=t+1;
            time=time+.1*deltaT;
            if time>=t_max-1
                runningODE=0;
            end
            y0Lead=y1Lead;
            y0Lagg=y1Lagg;

            if DataIntT(1)+time*10>=BleachTime && (BLead==0)
                BLead=y0Lead(2);
                BLagg=y0Lagg(4);

                FLead=1-BLead;
                FLagg=1-BLagg;

                ParametersLead(11)=FLead;
                ParametersLead(12)=BLead;

                ParametersLagg(11)=FLagg;
                ParametersLagg(12)=BLagg;

                y0Lead=[1,0,y0Lead(3:end)/FLead];
                y0Lagg=[y0Lagg(1:2)/FLagg,1,0,y0Lagg(5:6)/FLagg];
            end
            ylead(t,:)=y0Lead.*[BLead,FLead,BLead,FLead,BLead,FLead];
            ylagg(t,:)=y0Lagg.*[BLagg,FLagg,BLagg,FLagg,BLagg,FLagg];

            if BLead>0
                FRAPcurveFLead=[FRAPcurveFLead,y0Lead(2)/BLead];
                FRAPcurveFLagg=[FRAPcurveFLagg,y0Lagg(4)/BLagg];
            end

            if y0Lead(1)<0 || y0Lead(2)<0 ||y0Lead(3)<0 || y0Lead(4)<0 || y0Lead(5)<0 ||y0Lead(6)<0 || ~isreal(y0Lead(1)) || ~isreal(y0Lead(2)) || ~isreal(y0Lead(3)) || ~isreal(y0Lead(4)) || ~isreal(y0Lead(5)) || ~isreal(y0Lead(6))
                runningODE=0;
                BadParameters = 1;
            end
            if y0Lagg(1)<0 || y0Lagg(2)<0 ||y0Lagg(3)<0 || y0Lagg(4)<0 || y0Lagg(5)<0 ||y0Lagg(6)<0 || ~isreal(y0Lagg(1)) || ~isreal(y0Lagg(2)) || ~isreal(y0Lagg(3)) || ~isreal(y0Lagg(4)) || ~isreal(y0Lagg(5)) || ~isreal(y0Lagg(6))
                runningODE=0;
                BadParameters = 1;
            end
        end
        ylead=ylead(1:t,:);
        ylagg=ylagg(1:t,:);


        %% Check to see if this run was the best match to the data

        if ~BadParameters
            DistanceSimData = SimulationDataDistanceNormalizedFRAP(DataInt,[ylead(1:10/deltaT:t,1)+ylead(1:10/deltaT:t,2),ylead(1:10/deltaT:t,3)+ylead(1:10/deltaT:t,4),ylead(1:10/deltaT:t,5)+ylead(1:10/deltaT:t,6)],FRAPcurveFLead,FRAPcurveFLagg,recoveryTime,AverageLeadRecovery,AverageLaggRecovery);
            hits(1,:)=[array j DistanceSimData(1) DistanceSimData(2) DistanceSimData(3) DistanceSimData(4) DistanceSimData(5) ];
            parsave(Samplings, array,hits,DistanceSimData,j,mmhStyle)
            if DistanceSimData(1) < LeadCutoff && DistanceSimData(2) < LaggCutoff && DistanceSimData(3) < CytoCutoff && DistanceSimData(4) < FRAPCutoff && DistanceSimData(5) < FRAPCutoff
                runningMMH=0;
            end
        else
            hits(1,:)=[array j 10 10 10 10 10];
            parsave(Samplings, array,hits,[0 0 0 0 0],j,mmhStyle)
        end


        while runningMMH
            FRAPcurveFLead=[];
            FRAPcurveFLagg=[];
            BLead=0;
            FLead=1-BLead;
            BLagg=0;
            FLagg=1-BLagg;
            Candidate=mmh(array,MMHProp,1,mmhStyle);
            %Candidate(7)=max(1.2068*Candidate(6)-0.0026,0);
            NewParametersLead=[Candidate C_0 FLead BLead];
            NewParametersLagg=[Candidate C_0 FLagg BLagg];


            y0Lead=[C1b,C1f,C2b,C2f,Cyb,Cyf];
            y0Lagg=y0Lead;

            ylead=zeros(t_max,6);
            ylead(1,:)=y0Lead;

            ylagg=zeros(t_max,6);
            ylagg(1,:)=y0Lagg;

            runningODE=1;
            BadParameters = 0;
            t=1;
            time=0;
            while runningODE
                [~,y1Lead] = ode15s(@(t,y0Lead)ODE_DL_FRAP_Dynamic(t,NewParametersLead,y0Lead),[time-(0.1*deltaT),time],y0Lead);
                y1Lead=y1Lead(end,:);
                [~,y1Lagg] = ode15s(@(t,y0Lagg)ODE_DL_FRAP_Dynamic(t,NewParametersLagg,y0Lagg),[time-(0.1*deltaT),time],y0Lagg);
                y1Lagg=y1Lagg(end,:);

                t=t+1;
                time=time+.1*deltaT;
                if time>=t_max-1
                    runningODE=0;
                end
                y0Lead=y1Lead;
                y0Lagg=y1Lagg;

                if DataIntT(1)+time*10>=BleachTime && (BLead==0)
                    BLead=y0Lead(2);
                    BLagg=y0Lagg(4);

                    FLead=1-BLead;
                    FLagg=1-BLagg;

                    NewParametersLead(11)=FLead;
                    NewParametersLead(12)=BLead;

                    NewParametersLagg(11)=FLagg;
                    NewParametersLagg(12)=BLagg;

                    y0Lead=[1,0,y0Lead(3:end)/FLead];
                    y0Lagg=[y0Lagg(1:2)/FLagg,1,0,y0Lagg(5:6)/FLagg];
                end
                ylead(t,:)=y0Lead.*[BLead,FLead,BLead,FLead,BLead,FLead];
                ylagg(t,:)=y0Lagg.*[BLagg,FLagg,BLagg,FLagg,BLagg,FLagg];

                if BLead>0
                    FRAPcurveFLead=[FRAPcurveFLead,y0Lead(2)/BLead];
                    FRAPcurveFLagg=[FRAPcurveFLagg,y0Lagg(4)/BLagg];
                end

                if y0Lead(1)<0 || y0Lead(2)<0 ||y0Lead(3)<0 || y0Lead(4)<0 || y0Lead(5)<0 ||y0Lead(6)<0 || ~isreal(y0Lead(1)) || ~isreal(y0Lead(2)) || ~isreal(y0Lead(3)) || ~isreal(y0Lead(4)) || ~isreal(y0Lead(5)) || ~isreal(y0Lead(6))
                    runningODE=0;
                    BadParameters = 1;
                end
                if y0Lagg(1)<0 || y0Lagg(2)<0 ||y0Lagg(3)<0 || y0Lagg(4)<0 || y0Lagg(5)<0 ||y0Lagg(6)<0 || ~isreal(y0Lagg(1)) || ~isreal(y0Lagg(2)) || ~isreal(y0Lagg(3)) || ~isreal(y0Lagg(4)) || ~isreal(y0Lagg(5)) || ~isreal(y0Lagg(6))
                    runningODE=0;
                    BadParameters = 1;
                end
            end
            ylead=ylead(1:t,:);
            ylagg=ylagg(1:t,:);

            %% Check to see if this run was the best match to the data
            if ~BadParameters
                DistanceSimData = SimulationDataDistanceNormalizedFRAP(DataInt,[ylead(1:10/deltaT:t,1)+ylead(1:10/deltaT:t,2),ylead(1:10/deltaT:t,3)+ylead(1:10/deltaT:t,4),ylead(1:10/deltaT:t,5)+ylead(1:10/deltaT:t,6)],FRAPcurveFLead,FRAPcurveFLagg,recoveryTime,AverageLeadRecovery,AverageLaggRecovery);
                %hits(count,:)=[array j DistanceSimData(1) DistanceSimData(2) DistanceSimData(3)];
                %                 sum(DistanceSimData)-sum(hits(count,end-4:end))
                %if sum(DistanceSimData)<sum(hits(count,end-4:end))
                if BetterDistanceSum(DistanceSimData,hits(count,end-4:end),LeadCutoff,LaggCutoff,CytoCutoff,FRAPCutoff)
                    %disp('found a better fit')

                    count=count+1;
                    %                     figure(count)
                    %                     subplot(2,2,1)
                    %                     plot(ylead(:,1)+ylead(:,2))
                    %                     hold on
                    %                     plot(ylagg(:,1)+ylagg(:,1))
                    %                     subplot(2,2,2)
                    %                     plot(ylead(:,3)+ylead(:,4))
                    %                     hold on
                    %                     plot(ylagg(:,3)+ylagg(:,4))
                    %                     subplot(2,2,3)
                    %                     plot(ylead(:,5)+ylead(:,6))
                    %                     hold on
                    %                     plot(ylagg(:,5)+ylagg(:,6))
                    %                     hold off
                    array=Candidate;
                    hits(count,:)=[array j DistanceSimData(1) DistanceSimData(2) DistanceSimData(3) DistanceSimData(4) DistanceSimData(5)];
                    parsave(Samplings, array,hits,DistanceSimData,j,mmhStyle)
                    if DistanceSimData(1) < LeadCutoff && DistanceSimData(2) < LaggCutoff && DistanceSimData(3) < CytoCutoff && DistanceSimData(4) < FRAPCutoff && DistanceSimData(5) < FRAPCutoff
                        runningMMH=0;
                        disp('found fit')
                        Samplings
                    end
                end
            end
            %             FRAPcurveFLead(1:6)
            %             FRAPcurveFLagg(1:6)
            j=j+1;
            if mod(j,10000) == 0
                format long
                [Samplings; array'];
            end
            if j> Maxruns
                disp('maxruns')
                Samplings
                runningMMH=0;
            end
        end

        %% Report Result to compare to other models
        parsave(Samplings, array,hits,hits(end,end-4:end),j,mmhStyle)

    end
end

function parsave(Samplings, array,hits,DistanceSimData,j,mmhStyle)
save(sprintf('MMH_FRAP_Dynamic_Full_V2_dt1_style_%d/Run_%g.mat',mmhStyle, Samplings),'array','hits','DistanceSimData','j')
end

function [Better] = BetterDistanceSum(DistanceSimData,Allhits,LeadCutoff,LaggCutoff,CytoCutoff,FRAPCutoff)
Error=DistanceSimData;
if DistanceSimData(1)<LeadCutoff
    Error(1)=0;
    Allhits(1)=0;
end
if DistanceSimData(2)<LaggCutoff
    Error(2)=0;
    Allhits(2)=0;
end
if DistanceSimData(3)<CytoCutoff
    Error(3)=0;
    Allhits(3)=0;
end
if DistanceSimData(4)<FRAPCutoff
    Error(4)=0;
    Allhits(4)=0;
end
if DistanceSimData(5)<FRAPCutoff
    Error(5)=0;
    Allhits(5)=0;
end
Better = sum(Error)<=sum(Allhits);
end