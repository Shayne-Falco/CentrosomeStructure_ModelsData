close all
S=load('PhotobleachCorrection.mat');
DataDisplay = S.CorrectforPhotobleach;
t = DataDisplay(:,1);
TimeFrame = t>-331 & t<1;
t = DataDisplay(TimeFrame,1);
count = DataDisplay(TimeFrame,2);

Cytos = DataDisplay(TimeFrame,10:5:end);
Leads = DataDisplay(TimeFrame,8:5:end);
Laggs = DataDisplay(TimeFrame,9:5:end);
Totls = DataDisplay(TimeFrame,11:5:end);

% for i = 1:size(Totls,2)
%
%     AvgTot = mean(Totls(:,i));
%     Totls(:,i) = Totls(:,i)./AvgTot;
%     Cytos(:,i) = Cytos(:,i)./AvgTot;
%     Leads(:,i) = Leads(:,i)./AvgTot;
%     Laggs(:,i) = Laggs(:,i)./AvgTot;
% end
AvgTotl = mean(Totls,2);
AvgCyto = mean(Cytos,2);
AvgLead = mean(Leads,2);
AvgLagg = mean(Laggs,2);

plot(t,AvgCyto,'g',t,AvgLead,'r',t,AvgLagg,'b',t,AvgTotl,'m')
legend('Average Cytoplam Level','Average Leading Level','Average Lagging Level','Average Total Level')

maxLeadLagg = max(max(max(Leads)),max(max(Laggs)));
figure

plot(t,Leads,'black-')
hold on
plot(t,AvgLead,'r-','linewidth',4)
ylim([0 maxLeadLagg])
title(["Average Leading with individual timecourses"])
xlim([-330 0])
%subtitle()
figure

plot(t,Laggs,'black-')
hold on
plot(t,AvgLagg,'b-','linewidth',4)
ylim([0 maxLeadLagg])
title(["Average Lagging with individual timecourses"])
xlim([-330 0])
figure

plot(t,Leads,'r--','linewidth',.5,'HandleVisibility','Off')
hold on
plot(t,Laggs,'b--','linewidth',.5,'HandleVisibility','Off')
plot(t,AvgLead,'r-','linewidth',5,'DisplayName','Leading')
plot(t,AvgLagg,'b-','linewidth',5,'DisplayName','Lagging')
legend('Location','northwest')
ylim([0 maxLeadLagg])
ylabel('Normalized AIR-1::GFP signal')
xlabel('time (sec) 0=NEBD')
xlim([-330 0])
title("Average Leading and Lagging with individual timecourses")
figure

plot(t,Totls,'black--')
hold on
plot(t,AvgTotl,'b-','linewidth',3)
plot(t,ones(length(t),1),'r','linewidth',2)
ylim([0.6 1.4])
title(["Average Totals with individual timecourses"])
figure

plot(t,Leads'./Laggs','black--')
hold on
errorbar(t,mean(Leads'./Laggs','omitnan'),std(Leads'./Laggs','omitnan'),'r','linewidth',1)
plot(t,mean(Leads'./Laggs','omitnan'),'r','linewidth',4)
%plot(t,ones(length(t),1),'y--','linewidth',2)
ylim([0 2])
xlim([-330 0])
ylabel('Ratio')
xlabel('time (sec) 0=NEBD')
title("Average Ratio of Leading/Lagging with individual ratios")

figure
i=1;
for Bt = 4:2:26
    subplot(3,4,i)
    i=i+1;
    BleachLagg=[AvgLagg(Bt); AvgLagg(Bt+1:Bt+6)-AvgLagg(Bt+1)]./AvgLagg(Bt);
    BleachLead=[AvgLead(Bt); AvgLead(Bt+1:Bt+6)-AvgLead(Bt+1)]./AvgLead(Bt);
    plot(t(Bt:Bt+6),BleachLead,t(Bt:Bt+6),BleachLagg)
    if Bt==4
        legend('Lead','Lagg')
    end
end
sgtitle('FRAP with no Recovery')