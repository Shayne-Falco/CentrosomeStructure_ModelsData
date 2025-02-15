clear all
%close all
ToPlot = 0;
mmhstyle = 1;

titlesDyn={'a1','a2','a3','log(a4)','log(a5)','a6','a7','n1','n2','Thalf','Vmax'}; %Dynamic
titlesCon={'a1','a2','a3','a4','a5','n1','n2','C_0','Cy','t','Lead Error','Lagger Error','Cytoplasm Error','Total Error'}; %Constant and Conserved
titles={'a1','a2','a3','a4','a5','a6','a7','n1','n2','Thalf','Vmax'};
mymap = [
    1 0 0
    0 1 0
    0 0 1
    ];


S=load(['FRAPFitToBoth_mmh',num2str(mmhstyle),'_Both.mat']);
FinalParametersDynamic_Both = S.AnalysisParameters(:,1:12);

S=load(['FRAPFitToBoth_mmh',num2str(mmhstyle),'_Lagg.mat']);
FinalParametersDynamic_Lagg = S.AnalysisParameters(:,1:12);

S=load(['FRAPFitToBoth_mmh',num2str(mmhstyle),'_Lead.mat']);
FinalParametersDynamic_Lead = S.AnalysisParameters(:,1:12);




a1 = FinalParametersDynamic_Both(:,1);
a2 = FinalParametersDynamic_Both(:,2);
a3 = FinalParametersDynamic_Both(:,3);
a4 = FinalParametersDynamic_Both(:,4);
a5 = FinalParametersDynamic_Both(:,5);
a6 = FinalParametersDynamic_Both(:,6);
a7 = FinalParametersDynamic_Both(:,7);
n1 = FinalParametersDynamic_Both(:,8);
n2 = FinalParametersDynamic_Both(:,9);

N1mN2 = n1-n2;
figure
hh =histogram(N1mN2,25,'BinLimits',[-1,.25])
xlabel('n_1 - n_2')
sum(N1mN2>0)/(sum(N1mN2>0)+sum(N1mN2<0))
figure
scatter(n1,n2)
xlim([1 2])
ylim([1 2])
[RHO,PVAL] = corr(n1,n2,'Type','Spearman');

HistBinEdges = hh.BinEdges;
HistBinEdges = HistBinEdges(1:end-1)+hh.BinWidth;
HistData = [HistBinEdges' hh.Values'];
filename = append('FRAPcorFRAPFit_Full/mmh_',num2str(mmhstyle),'/Histdata_n1mn2.dat');
writematrix(HistData,string(filename),'filetype', 'text','Delimiter','tab')

maximumParam = max(max(FinalParametersDynamic_Lead),max(max(FinalParametersDynamic_Both),max(FinalParametersDynamic_Lagg)));
BinNumber = 16;
for LeadBleach = 2
    if LeadBleach ==0
        FinalParametersDynamic=FinalParametersDynamic_Lagg;
    elseif LeadBleach == 1
        FinalParametersDynamic=FinalParametersDynamic_Lead;
    else
        FinalParametersDynamic=FinalParametersDynamic_Both;
    end
    f=figure(LeadBleach+1);
    for i = 1:9
        HistData = [];
        %FinalParametersDynamic(FinalParametersDynamic>2000)=2000;
        subplot(3,3,i)
            ToHistPlot = FinalParametersDynamic(:,i);
            %ToHistPlot(ToHistPlot>90)=90;
        if i==8 || i==9
            hh = histogram(ToHistPlot,'NumBins', BinNumber,'BinLimits',[0.5,2.5]);
        elseif i==4 || i==5
            hh = histogram(log(ToHistPlot),'NumBins', BinNumber);
        elseif i==3
            hh = histogram(ToHistPlot,'NumBins', BinNumber,'BinLimits',[0,5]);
        else
            hh = histogram(ToHistPlot,'NumBins', BinNumber);
        end
        title(titlesDyn(i))
        %ylim([0 16])
        HistBinEdges = hh.BinEdges;
        HistBinEdges = HistBinEdges(1:end-1)+hh.BinWidth;
        HistData = [HistBinEdges' hh.Values'];
        filename = append('Paper_Data/Hists/FitFRAP_mmh',num2str(mmhstyle),'/Histdata_',titlesDyn(i),'.dat');
        writematrix(HistData,string(filename),'filetype', 'text','Delimiter','tab')
    end
    if LeadBleach ==0
        sgtitle('Parameter Sets That Match to Lagging Recovery Curve')
        saveas(f,append('FRAPcorFRAPFit_Full/mmh_',num2str(mmhstyle),'/Lagg/BarPlot.jpg'));
    elseif LeadBleach == 1
        sgtitle('Parameter Sets That Match to Leading Recovery Curve')
        saveas(f,append('FRAPcorFRAPFit_Full/mmh_',num2str(mmhstyle),'/Lead/BarPlot.jpg'));
    else
        sgtitle('Parameter Sets That Match to both Recovery Curves')
        saveas(f,append('FRAPcorFRAPFit_Full/mmh_',num2str(mmhstyle),'/Both/BarPlot.jpg'));
        
    end
    
    if ToPlot
        figureCount=4;
        for i=1:8
            for j=i+1:9
                f = figure(figureCount);
                f.Visible = 'off';
                colormap(mymap)
                x=FinalParametersDynamic(:,i);
                y=FinalParametersDynamic(:,j);
                scatter(x,y,[],FinalParametersDynamic(:,12)*(LeadBleach+1))
                hold on
                if LeadBleach ==2
                    xlabel(titles(i))
                    ylabel(titles(j))
                    %if i==8 && j==9
                        %plot(0:.1:2.5,0:.1:2.5)
                        legend('Lagging','Leading','Both','Location','best')
                    %else
                    %    legend('Lagging','Leading','Both','Location','best')
                    %end
                    if i == 8 || i == 9
                        xlim([0 2.5])
                    else
                        xlim([0 min(maximumParam(i)+.01,41)])
                    end
                    if j == 8 || j == 9
                        ylim([0 2.5])
                    else
                        ylim([0 min(maximumParam(j)+.01,41)])
                    end

                    toSave = [x y];
                    writematrix(toSave,append('FRAPcorFRAPFit_Full/mmh_',num2str(mmhstyle),'/Both/',string(titles(i)),string(titles(j)),'.dat'),'Delimiter','tab')
                    saveas(f,append('FRAPcorFRAPFit_Full/mmh_',num2str(mmhstyle),'/All/',string(titles(i)),string(titles(j))),'png');
                end

                figureCount = figureCount+1;
            end
        end



        for i=1:8
            for j=i+1:9
                f = figure('visible','off');
                colormap flag
                x=FinalParametersDynamic(:,i);
                y=FinalParametersDynamic(:,j);
                scatter(x,y,[],FinalParametersDynamic(:,12))
                xlabel(titles(i))
                ylabel(titles(j))
                if LeadBleach == 1
                    title('Leading recovery similar to data')
                elseif LeadBleach == 0
                    title('Lagging recovery similar to data')
                else
                    title('Both recoveries similar to data')
                end
                if i==8 && j==9
                    hold on
                    %plot(0:.1:3,0:.1:3)
                    legend(sprintf('n_1>n_2; %0.2f%%', 100*sum(x>y)/length(y)),'Location','best');
                end
                if i == 8 || i == 9
                    xlim([0 3])
                else
                    xlim([0 min(maximumParam(i)+.01,41)])
                end
                if j == 8 || j == 9
                    ylim([0 3])
                else
                    ylim([0 min(maximumParam(j)+.01,41)])
                end
                if LeadBleach ==1
                    saveas(f,append('FRAPcorFRAPFit_Full/mmh_',num2str(mmhstyle),'/Lead/',string(titles(i)),string(titles(j))),'jpg');
                    %savefig(f,append('FRAPcorFRAPFit_Full/mmh_',num2str(mmhstyle),'/Lead/',string(titles(i)),string(titles(j))));
                elseif LeadBleach == 0
                    saveas(f,append('FRAPcorFRAPFit_Full/mmh_',num2str(mmhstyle),'/Lagg/',string(titles(i)),string(titles(j))),'jpg');
                    %savefig(f,append('FRAPcorFRAPFit_Full/mmh_',num2str(mmhstyle),'/Lagg/',string(titles(i)),string(titles(j))));
                else
                    saveas(f,append('FRAPcorFRAPFit_Full/mmh_',num2str(mmhstyle),'/Both/',string(titles(i)),string(titles(j))),'jpg');
                    %saveas(f,append('FRAPcorFRAPFit_Full/mmh_',num2str(mmhstyle),'/Both/',string(titles(i)),string(titles(j))));
                end

            end
        end
    end
end