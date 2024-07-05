clear all
Maxruns = 1500000;
LeadCutoff = 0.2;
LaggCutoff = 0.15;
CytoCutoff = 0.55;


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


%% Set up model


MMHProp=2;
t_max=lengthDataInt;


C1=DataIntLead(1); %leading MTOC initial value
C2=DataIntLagg(1); %Lagging MTOC initial value
c0min = .75;
c0max = 1.25;

TotalSamples=1:50;
for mmhStyle = 1:3
    parfor (Samplings=TotalSamples,4)
        
        runningMMH=1;
        myfilename=['MMH_Conserved_style_',num2str(mmhStyle),'/Run_',num2str(Samplings), '.mat'];
        if isfile(myfilename)
            continue
        end
        
        C_0 = c0min + (c0max-c0min).*rand(1);
        Cy = C_0-C1-C2; 
        ConstantCy=0; 
        
        n1 = 4*rand(1);   	%
        n2 = 4*rand(1);
        a_1 = 1*rand(1); 	%
        a_2 = 1*rand(1); 	%
        a_3 = 1*rand(1); 	%
        a_4 = 1*rand(1); 	%
        a_5 = 1*rand(1); 	%
        
        array=[a_1 a_2 a_3 a_4 a_5 n1 n2 C_0];
        count=1;
        hits=[];
        NewParameters=[array,ConstantCy];
        hits(1,:)=[NewParameters 1 inf inf inf];
        j=2;
        
        
        y0=[C1,C2];
        y=zeros(t_max,3);
        y(1,:)=[y0,Cy];
        runningODE=1;
        BadParameters = 0;
        t=1;
        while runningODE
            [~,y1] = ode15s(@(t,y0)ODE_DL_3C(t,NewParameters,y0),[t-1,t],y0);
            y1=y1(end,:);
%             [~,y1] = ode15s(@(t,y1)ODE_DL_3C(t,NewParameters,y1),[t-.5,t],y1);
%             y1=y1(end,:);
            t=t+1;
            y0=y1;
            y(t,:)=[y0,C_0-sum(y0)];
            if t==t_max
                runningODE=0;
            end
            if y0(1)<0 || y0(2)<0 || ~isreal(y0(1)) || ~isreal(y0(2))
                runningODE=0;
                BadParameters = 1;
            end
        end
        y=y(1:t,:);
       
        
        %% Check to see if this run was the best match to the data
        if ~BadParameters
            DistanceSimData = SimulationDataDistanceNormalized(DataInt,y);
            hits(1,:)=[NewParameters j DistanceSimData(1) DistanceSimData(2) DistanceSimData(3)];
            parsave(Samplings, array,hits,DistanceSimData,j,mmhStyle)
            if DistanceSimData(1) < LeadCutoff && DistanceSimData(2) < LaggCutoff && DistanceSimData(3) < CytoCutoff
                runningMMH=0;
            end
        else
            DistanceSimData = [];
            hits(1,:)=[NewParameters j 10 10 10];
            parsave(Samplings, array,hits,DistanceSimData,j,mmhStyle)
        end
        
        while runningMMH
            Candidate=mmh(array,MMHProp,1,mmhStyle);
            if Candidate(end)<c0min
                Candidate(end)=c0min;
            elseif Candidate(end)>c0max
                Candidate(end)=c0max;
            end
            C_0=Candidate(end);
            NewParameters=[Candidate,ConstantCy];
            y0=[C1,C2];
            y=zeros(t_max,3);
            y(1,:)=[y0,C_0-C1-C2];
            runningODE=1;
            BadParameters = 0;
            t=1;
            while runningODE
                [~,y1] = ode15s(@(t,y0)ODE_DL_3C(t,NewParameters,y0),[t-1,t],y0);
                y1=y1(end,:);
%                 [~,y1] = ode15s(@(t,y1)ODE_DL_3C(t,NewParameters,y1),[t-.5,t],y1);
%                 y1=y1(end,:);
                t=t+1;
                deltay=y1-y0;
                y0=y1;
                y(t,:)=[y0,C_0-sum(y0)];
                if t==t_max
                    runningODE=0;
                end
                if y0(1)<0 || y0(2)<0 || ~isreal(y0(1)) || ~isreal(y0(2))
                    runningODE=0;
                    BadParameters = 1;
                end
            end
            y=y(1:t,:);
            
            %% Check to see if this run was the best match to the data
            if ~BadParameters
                DistanceSimData = SimulationDataDistanceNormalized(DataInt,y);
                %hits(count,:)=[array j DistanceSimData(1) DistanceSimData(2) DistanceSimData(3)];
                if DistanceSimData(1)+DistanceSimData(2)+DistanceSimData(3)<hits(count,end-2)+hits(count,end-1)+hits(count,end)
                    count=count+1;
                    array=Candidate;
                    hits(count,:)=[NewParameters j DistanceSimData(1) DistanceSimData(2) DistanceSimData(3)];
                    parsave(Samplings, array,hits,DistanceSimData,j,mmhStyle)
                    if DistanceSimData(1)<LeadCutoff && DistanceSimData(2)<LaggCutoff && DistanceSimData(3) < CytoCutoff
                        runningMMH=0;
                    end
                end
            end
            j=j+1;
            if j> Maxruns
                disp('maxruns')
                runningMMH=0;
            end
        end
        
        %% Report Result to compare to other models
        parsave(Samplings, array,hits,hits(end,end-3:end),j,mmhStyle)
        
    end
end

function parsave(Samplings, array,hits,DistanceSimData,j,mmhStyle)
save(sprintf('MMH_Conserved_style_%d/Run_%g.mat',mmhStyle, Samplings),'array','hits','DistanceSimData','j')
end