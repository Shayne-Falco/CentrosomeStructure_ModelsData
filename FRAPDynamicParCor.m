
maximumsPlot = [1,.25,10,10,10,.12,.12,3,3,60,2];

titles={'a1','a2','a3','a4','a5','a6','a7','n1','n2','Thalf','Vmax'};
for i=1:10
    for j=i+1:11
        f = figure('visible','off');
        colormap flag
        x=AnalysisParametersLead(:,i);
        y=AnalysisParametersLead(:,j);
        scatter(x,y,[],(AnalysisParametersLead(:,12)==1) & (AnalysisParametersLagg(:,12)==1))
        xlabel(titles(i))
        ylabel(titles(j))
        title('Black = Both recovery similar to data')
        if i==8 && j==9
            hold on
            plot(0:.1:3,0:.1:3)
            %legend('Data',sprintf('n_1>n_2; %0.2f%%', 100*sum(x>y)/length(y)),'Location','best');
        end
        if i == 8 || i == 9
            xlim([0 3])
        else
            xlim([0 min(maximumsPlot(i),max(x))+.1])
        end
        if j == 8 || j == 9
            ylim([0 3])
        else
            ylim([0 min(maximumsPlot(j),max(y))+.1])
        end
        
        saveas(f,append('FRAPcor/',string(titles(i)),string(titles(j)),'Both'),'png');
    end
end

