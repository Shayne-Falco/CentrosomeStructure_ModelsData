
clear all
%close all
format shortg
ConstantCy=0;
mmhStyles = 1:3;
toDelete=0;
if ConstantCy
    preLabel = 'Paper_Data/Hists/FitRecruit_Constant/';
else
    preLabel = 'Paper_Data/Hists/FitRecruit_Conserved/';
end

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


toPlotAll=0;
%% Load all MMH sets ran
count=1;
running=1;

LeadErrorCutoff=.25;
LaggErrorCutoff=.25;
CytoErrorCutoff=1;
count=0;
FileExists=[];
LeadError=[0,100];
LagError=[0,100];
CytoError=[0,100];
WholeError=[0,100];
for mmhStyle = mmhStyles
    for run=1:2000
        if ConstantCy>0
            myfilename=['MMH_Constant_style_',num2str(mmhStyle),'/Run_',num2str(run), '.mat'];
        else
            myfilename=['MMH_Conserved_style_',num2str(mmhStyle),'/Run_',num2str(run), '.mat'];
        end
        if isfile(myfilename)
            count=count+1;
            load(myfilename)
            results{count}=hits;
            FileExists=[FileExists run];
        end

    end
end

AllParameters=[];
FinalParameters=zeros(count,10);
passed=[];
for resu=1:count
    %% Plotting results from multiple MMH runs
    
    hits=results{resu};
    
    Leaderrors=hits(end,end-2);
    Laggerrors=hits(end,end-1);
    CytoErrors=hits(end,end);
    
    
    if Leaderrors<LeadErrorCutoff && Laggerrors<LaggErrorCutoff && CytoErrors < CytoErrorCutoff
        passed=[passed resu];
        FileExists(resu)=-1;
        errors=hits(end,end-2:end);
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
    end
    
end
deleteList = FileExists(FileExists>0);

SimIntLead = [];
SimIntLagg = [];
SimIntCyto = [];
for resu=passed
    
    
    
    %% Plotting results from multiple MMH runs
    
    hits=results{resu};
    array=hits(end,1:8);
    C_0 = array(8); %
    Cy = C_0-C1-C2;
    
    ConstantCy=hits(end,9); %0 for conservation, else that is Cy value
    NewParameters=[array,ConstantCy];
    AllParameters=[AllParameters  ;  array];
    
    
    y0=[C1,C2];
    y=zeros(t_max,3);
    y(1,:)=[y0,Cy];
    
    
    running=1;
    t=1;
    %     rates(resu,t) = RateOfChange(NewParameters,y(1,:));
    %     PropRates(resu,t) = RateOfChange(NewParameters,y0)/sum(y(1,:));
    while running
        [~,y1] = ode15s(@(t,y0)ODE_DL_3C(t,NewParameters,y0),[t-1,t-.5],y0);
        y1=y1(end,:);
        [~,y1] = ode15s(@(t,y1)ODE_DL_3C(t,NewParameters,y1),[t-.5,t],y1);
        y1=y1(end,:);
        t=t+1;
        if t==t_max
            running=0;
        elseif norm(y1-y0)<0.00001 && t>=length(DataIntLagg)
            running=0;
        end
        y0=y1;
        
        if ConstantCy
            y(t,:)=[y0,Cy];
        else
            y(t,:)=[y0,C_0-sum(y0)];
        end
        %         rates(resu,t) = RateOfChange(NewParameters,y(t,:));
        %         PropRates(resu,t) = RateOfChange(NewParameters,y0)/sum(y(t,:));
    end
    
    DistanceSimData = SimulationDataDistanceNormalized(DataInt,y);
    SimIntLead = [SimIntLead, y(:,1)];
    SimIntLagg = [SimIntLagg y(:,2)];
    SimIntCyto = [SimIntCyto y(:,3)];
    
    
    for j=1:14
        if j==14
            FinalParameters(resu,j)=sum(hits(end,end-2:end));
        else
            FinalParameters(resu,j)=hits(end,j);
        end
    end
    if resu==WholeError(1)
        hits(end,:)
        figure
        subplot(2,2,1)
        plot(DataIntT,y(:,1),'black')
        hold on
        plot(DataIntT,DataIntLead,'r','LineWidth',2)
        ylim([0 .2])
        title(sprintf('Leading; %0.4f%',DistanceSimData(1)))
        
        subplot(2,2,2)
        plot(DataIntT,y(:,2),'black')
        hold on
        plot(DataIntT,DataIntLagg,'r','LineWidth',2)
        ylim([0 .2])
        title(sprintf('Lagging; %0.4f%',DistanceSimData(2)))
        
        subplot(2,2,3)
        plot(DataIntT,y(:,3),'black')
        hold on
        plot(DataIntT,DataIntCyto,'r','LineWidth',2)
        ylim([.75 .95])
        title(sprintf('Cyto; %0.4f%',DistanceSimData(3)))
        
        subplot(2,2,4)
        plot(DataIntT,y(:,1)+y(:,2)+y(:,3),'black')
        hold on
        plot(DataIntT,DataIntCyto+DataIntLagg+DataIntLead,'r','LineWidth',2)
        ylim([.9 1.1])
        title(sprintf('Total; %0.4f%',sum(DistanceSimData)))
        sgtitle('Lowest total Error')
        figure
        subplot(1,2,1)
        plot(hits(end-50:end,10),hits(end-50:end,6:7),'-*')
%         hold on
%         plot(hits(end-50:end,10),hits(end-50:end,7),'-*')
        xlabel('step')
        ylabel('Parameter Value')
%         legend('N1','N2')
        subplot(1,2,2)
        plot(hits(end-50:end,10),sum(hits(end-50:end,end-2:end),2),'-*')
        xlabel('step')
        ylabel('error')
        title('Total error over annealing steps')
        
    end
    
end
figure
plot(DataIntT,SimIntLead,'black')
hold on
plot(DataIntT,DataIntLead,'r','LineWidth',2)
ylim([0 .2])
title('Leading')

figure
plot(DataIntT,SimIntLagg,'black')
hold on
plot(DataIntT,DataIntLagg,'r','LineWidth',2)
ylim([0 .2])
title('Lagging')

figure
plot(DataIntT,SimIntCyto,'black')
hold on
plot(DataIntT,DataIntCyto,'r','LineWidth',2)
ylim([.75 .95])
title('Cyto')

figure
plot(DataIntT,SimIntCyto+SimIntLead+SimIntLagg,'black')
hold on
plot(DataIntT,DataIntCyto+DataIntLagg+DataIntLead,'r','LineWidth',2)
ylim([.9 1.1])
title('Total')

SimIntTotl = [DataIntT SimIntLagg+SimIntLead+SimIntCyto];
SimIntLead = [DataIntT SimIntLead];
SimIntLagg = [DataIntT SimIntLagg];
SimIntCyto = [DataIntT SimIntCyto];

figure
FinalParameters=FinalParameters(find(FinalParameters(:,11)),:);
titles={'a1','a2','a3','a4','a5','n1','n2','C_0','Cy','t','Lead Error','Lagger Error','Cytoplasm Error','Total Error'};
BinNumber = 16;
for i=[1:8 10:14]
    if i<9
        subplot(5,3,i)
    else
        subplot(5,3,i-1)
    end
    ToHistPlot = FinalParameters(:,i);
    %ToHistPlot(ToHistPlot>90)=90;
    if i==6 || i==7
        hh = histogram(ToHistPlot,'NumBins', BinNumber,'BinLimits',[0.5,2.5]);
    elseif i==4 || i==5
        hh = histogram(ToHistPlot,'NumBins', BinNumber);
    elseif i==3
        hh = histogram(ToHistPlot,'NumBins', BinNumber,'BinLimits',[0,5]);
    elseif i ==11 || i == 12
        hh = histogram(ToHistPlot,'NumBins', BinNumber*2,'BinLimits',[0,.25]);
    elseif i>12
        hh = histogram(ToHistPlot,'NumBins', BinNumber*2,'BinLimits',[0,1.25]);
    else
        hh = histogram(ToHistPlot,'NumBins', BinNumber);
    end
    title(titles(i))
    HistBinEdges = hh.BinEdges;
    HistBinEdges = HistBinEdges(1:end-1)+hh.BinWidth;
    HistData = [HistBinEdges' hh.Values'];
    filename = append(preLabel,'Histdata_',titles(i),'.dat');
    writematrix(HistData,string(filename),'filetype', 'text','Delimiter','tab')
end

N1mN2 = FinalParameters(:,6)-FinalParameters(:,7);
figure
hh =histogram(N1mN2,26,'BinLimits',[-1.75,1.75]);
xlabel('n_1 - n_2')
sum(N1mN2>0)/(sum(N1mN2>0)+sum(N1mN2<0))
figure
scatter(FinalParameters(:,6),FinalParameters(:,7))
xlim([1 3])
ylim([1 3])
[RHO,PVAL] = corr(FinalParameters(:,6),FinalParameters(:,7),'Type','Spearman');

HistBinEdges = hh.BinEdges;
HistBinEdges = HistBinEdges(1:end-1)+hh.BinWidth;
HistData = [HistBinEdges' hh.Values'];
filename = append(preLabel,'Histdata_n1mn2.dat');
writematrix(HistData,string(filename),'filetype', 'text','Delimiter','tab')

sgtitle('Constant Cytoplasm')



MinParamaters=min(FinalParameters);
MaxParamaters=max(FinalParameters);
length(FinalParameters);

%save('LHSRange_300to0.mat','MinParamaters','MaxParamaters','FinalParameters')
for i=6
    for j=7
        figure;
        x=FinalParameters(:,i);
        y=FinalParameters(:,j);
        %         omitx=x>50;
        %         omity=y>50;
        %         omit=omitx|omity;
        %         x=x(~omit);
        %         y=y(~omit);
        %         sum(x>y)/length(y);
        %         b1=x\y;
        %         yCalc1=b1*x;
        scatter(x,y)
        hold on
        %plot(x,yCalc1)
        %         X = [ones(length(x),1) x];
        %         b = X\y;
        %         yCalc2 = X*b;
        %         Rsq = 1 - sum((y - yCalc2).^2)/sum((y - mean(y)).^2);
        %         plot(x,yCalc2,'--','LineWidth',2)
        if i==6 && j==7
            plot(0:.1:5,0:.1:5)
            legend('Data',sprintf('n_1>n_2; %0.2f%%', 100*sum(x>y)/length(y)),'Location','best');
            %         else
            %             legend('Data',sprintf('R^2= %g', Rsq),'Location','best');
        end
        %plot(MatrixOfDeltas(:,i),yCalc1)
        xlabel(titles(i))
        ylabel(titles(j))
        title('Constant');
    end
end

[h,p,ci,stats]=ttest(FinalParameters(:,6),FinalParameters(:,7))