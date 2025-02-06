
clear all
close all
format shortg

mmhStyle = 1;
toDelete=0;

%% Load all MMH sets ran

running=1;

LeadErrorCutoff=.2;
LaggErrorCutoff=.2;
CytoErrorCutoff=.3;
count=0;
FileExists=[];
LeadError=[0,100];
LagError=[0,100];
CytoError=[0,100];
WholeError=[0,100];
LeadLagError=[0,100];
for run=1:2000
    
    myfilename=['MMH_FRAP_Dynamic_dt1_SymmetricLeadLagg_Average',num2str(run),'/Run_',num2str(run), '.mat'];
    if isfile(myfilename)
        count=count+1;
        load(myfilename)
        results{count}=hits;

        FileExists=[FileExists run];
    end
    
end

AllParameters=[];
FinalParameters=zeros(count,10);
passed=[];
for resu=1:count
    %% Plotting results from multiple MMH runs
    
    hits=results{resu};
    
    Leaderrors=hits(end,end-4);
    Laggerrors=hits(end,end-3);
    CytoErrors=hits(end,end-2);
    
    
    if Leaderrors<LeadErrorCutoff && Laggerrors<LaggErrorCutoff && CytoErrors < CytoErrorCutoff
        passed=[passed resu];
        %FileExists(resu)=-1;
        errors=hits(end,end-4:end-2);
        if errors(3)<CytoError(2)
            CytoError=[resu,errors(3)];
        end
        if errors(2)<LagError(2)
            LagError=[resu,errors(2)];
        end
        if errors(1)<LeadError(2)
            LeadError=[resu,errors(1)];
        end
        if sum(errors) < WholeError(2)
            WholeError = [resu,sum(errors)];
        end
        if errors(1)+errors(2) < LeadLagError(2)
            LeadLagError = [resu, errors(1)+errors(2)];
        end
    end
    
end
% deleteList = FileExists(FileExists>0);
% 
% if ~isempty(deleteList) && toDelete
%     for i = 1:length(deleteList)
%         todeletefilename = ['MMH_Dynamic_ExpNorm_300to0_style_',num2str(mmhStyle),'/Run_',num2str(deleteList(i)), '.mat'];
%         delete(todeletefilename)
%     end
% end
SimIntLead = [];
SimIntLagg = [];
SimIntCyto = [];
for resu=passed
    %% Load data to campare model to
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
    t_max=lengthDataInt;
    
    C1 = DataIntLead(1); %leading MTOC initial value
    C2 = DataIntLagg(1); %Lagging MTOC initial value
    Cy = DataIntCyto(1);
    C_0 = C1+C2+Cy; %
    
    %% Plotting results from multiple MMH runs
    
    hits=results{resu};
    array=hits(end,1:9);
    NewParameters=[array,C_0];
    AllParameters=[AllParameters  ;  array];
    y0=[C1,C2,Cy];
    y=zeros(t_max,3);
    y(1,:)=y0;
    running=1;
    t=1;
    while running
        [~,y1] = ode15s(@(t,y0)ODE_DL_3C_Dynamic(t,NewParameters,y0),[t-1,t],y0);
%         y1=y1(end,:);
%         [~,y1] = ode15s(@(t,y1)ODE_DL_3C_Dynamic(t,NewParameters,y1),[t-.5,t],y1);
        y1=y1(end,:);
        t=t+1;
        if t==t_max
            running=0;
        end
        y0=y1;
        y(t,:)=y0;
    end
    
    DistanceSimData = SimulationDataDistanceNormalized(DataInt,y);
    SimIntLead = [SimIntLead y(:,1)];
    SimIntLagg = [SimIntLagg y(:,2)];
    SimIntCyto = [SimIntCyto y(:,3)];
    if resu==WholeError(1)
        NewParameters
        figure
        subplot(2,2,1)
        plot(DataIntT,DataIntLead,'r-','LineWidth',4)
        hold on
        plot(DataIntT,y(:,1),'black','LineWidth',3)
        ylim([0 .2])
        %title(sprintf('Leading; %0.4f%',DistanceSimData(1)))
        title('Leading Compartment/Centrosome ')
        legend('Normalized AIR-1::GFP signal','Simulation')
        
        subplot(2,2,2)
        
        plot(DataIntT,DataIntLagg,'r','LineWidth',4)
        hold on
        plot(DataIntT,y(:,2),'black','LineWidth',3)
        ylim([0 .2])
        title('Lagging Compartment/Centrosome ')
        %title(sprintf('Lagging; %0.4f%',DistanceSimData(2)))
        
        subplot(2,2,3)
        
        plot(DataIntT,DataIntCyto,'r','LineWidth',4)
        hold on
        plot(DataIntT,y(:,3),'black','LineWidth',3)
        ylim([.75 .95])
        title('Cytoplasm')
        xlabel('time (sec) from NEBD')
        %title(sprintf('Cyto; %0.4f%',DistanceSimData(3)))
        
        subplot(2,2,4)
        plot(DataIntT,DataIntCyto+DataIntLagg+DataIntLead,'r','LineWidth',4)
        hold on
        plot(DataIntT,y(:,1)+y(:,2)+y(:,3),'black','LineWidth',3)
        ylim([.9 1.1])
        title('Total')
        xlabel('time (sec) from NEBD')
        %title(sprintf('Total; %0.4f%',sum(DistanceSimData)))
        sgtitle('Dynamic Cytoplasm Model')
        
%         figure
%         subplot(1,2,1)
%         plot(hits(end-25:end,10),hits(end-25:end,1:9),'-*')
% %         hold on
% %         plot(hits(80:end,10),hits(80:end,9),'-*')
%         xlabel('step')
%         ylabel('Parameter Value')
% %         legend('N1','N2')
%         subplot(1,2,2)
%         plot(hits(end-25:end,10),sum(hits(end-25:end,end-2:end),2),'-*')
%         xlabel('step')
%         ylabel('error')
%         title('Total error over annealing steps')
    end
    for j=1:14
        if j==14
            FinalParameters(resu,j)=sum(hits(end,end-2:end));
        else
            FinalParameters(resu,j)=hits(end,j);
        end
    end
    
    
end
figure
subplot(2,2,1)
plot(DataIntT,SimIntLead,'black')
hold on
plot(DataIntT,DataIntLead,'r','LineWidth',2)
ylim([0 .2])
title('Leading')

subplot(2,2,2)
plot(DataIntT,SimIntLagg,'black')
hold on
plot(DataIntT,DataIntLagg,'r','LineWidth',2)
ylim([0 .2])
title('Lagging')

subplot(2,2,3)
plot(DataIntT,SimIntCyto,'black')
hold on
plot(DataIntT,DataIntCyto,'r','LineWidth',2)
ylim([.75 .95])
title('Cyto')

subplot(2,2,4)
plot(DataIntT,SimIntCyto+SimIntLead+SimIntLagg,'black')
hold on
plot(DataIntT,DataIntCyto+DataIntLagg+DataIntLead,'r','LineWidth',2)
ylim([.9 1.1])
title('Total')

figure
FinalParameters=FinalParameters(find(FinalParameters(:,10)),:);
titles={'a1','a2','a3','a4','a5', 'a6','a7','n1','n2','t','Lead Error','Lagger Error','Cytoplasm Error','Total Error'};
for i=1:14
    subplot(4,4,i)
    if i==10
        histogram(FinalParameters(:,i),15)
    else
        toplot=FinalParameters(:,i);
        toplot(toplot>5)=ones(length(find(toplot>5)),1)*-1;
        histogram(toplot,15)
    end
    title(titles(i))
end
sgtitle('Dynamic Cytoplasm')
% filename=strcat('corr_300to0/Style_',num2str(mmhStyle),'Hist.fig');
% savefig(filename)
% filename=strcat('corr_300to0/Style_',num2str(mmhStyle),'Hist.png');
% saveas(gcf,filename)

MinParamaters=min(FinalParameters);
MaxParamaters=max(FinalParameters);
length(FinalParameters);
%titles=['a1','a2','a3','a4','a5','n1','n2','t','Lead Error','Lagger Error','Total Error'];
% save('LHSRange_300to0.mat','MinParamaters','MaxParamaters','FinalParameters')
% for i=1:8
%     for j=i+1:9
%         figure;
%         x=FinalParameters(:,i);
%         y=FinalParameters(:,j);
%         omitx=x>50;
%         omity=y>50;
%         omit=omitx|omity;
%         x=x(~omit);
%         y=y(~omit);
%         sum(x>y)/length(y);
%         b1=x\y;
%         yCalc1=b1*x;
%         scatter(x,y)
%         hold on
%         %plot(x,yCalc1)
%         X = [ones(length(x),1) x];
%         b = X\y;
%         yCalc2 = X*b;
%         Rsq = 1 - sum((y - yCalc2).^2)/sum((y - mean(y)).^2);
%         plot(x,yCalc2,'--','LineWidth',1)
%         legend('Data',sprintf('R^2= %g \n slope = %g \n intercept = %g', Rsq,b(2),b(1)),'Location','best');
%         xlabel(titles(i))
%         ylabel(titles(j))
%         title('Linear Correlation Test');
%     end
% end
for i=8
    for j=9
        figure;
        x=FinalParameters(:,i);
        y=FinalParameters(:,j);
        scatter(x,y)
        hold on
        if i==8 && j==9
            plot(0:.1:5,0:.1:5)
            legend('Data',sprintf('n_1>n_2; %0.2f%%', 100*sum(x>y)/length(y)),'Location','best');
        end
        xlabel(titles(i))
        ylabel(titles(j))
        title('Constant');
    end
end


[h,p,ci,stats]=ttest(FinalParameters(:,8),FinalParameters(:,9))