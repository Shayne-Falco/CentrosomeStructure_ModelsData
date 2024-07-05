clear all
close all

titlesDyn={'a1','a2','a3','a4','a5', 'a6','a7','n1','n2','t','Lead Error','Lagger Error','Cytoplasm Error','Total Error'}; %Dynamic
titlesCon={'a1','a2','a3','a4','a5','n1','n2','C_0','Cy','t','Lead Error','Lagger Error','Cytoplasm Error','Total Error'}; %Constant and Conserved

S=load('FinalParametersDynamic1.mat');
FinalParameters1 = S.FinalParameters;
S=load('FinalParametersDynamic2.mat');
FinalParameters2 = S.FinalParameters;
S=load('FinalParametersDynamic3.mat');
FinalParameters3 = S.FinalParameters;

n1loc = 8; %6 - Con; 8 - Dyn
n2loc = n1loc+1;

BinNumber = 20;

figure

subplot(3,3,1)
histogram(FinalParameters1(:,n1loc),'NumBins', BinNumber,'BinLimits',[0,4],'FaceColor','r');
ylabel('n1')
title('mmh original objective')
%ylim([0 35])
xticks([0 1 2 3 4])

subplot(3,3,2)
histogram(FinalParameters2(:,n1loc),'NumBins', BinNumber,'BinLimits',[0,4],'FaceColor','r');
title('mmh Constant objective')
%ylim([0 35])
xticks([0 1 2 3 4])

subplot(3,3,3)
histogram(FinalParameters3(:,n1loc),'NumBins', BinNumber,'BinLimits',[0,4],'FaceColor','r');
title('mmh no objective')
%ylim([0 35])
xticks([0 1 2 3 4])

subplot(3,3,4)
histogram(FinalParameters1(:,n2loc),'NumBins', BinNumber,'BinLimits',[0,4]);
%ylim([0 35])
ylabel('n2')
xticks([0 1 2 3 4])

subplot(3,3,5)
histogram(FinalParameters2(:,n2loc),'NumBins', BinNumber,'BinLimits',[0,4]);
%ylim([0 35])
xticks([0 1 2 3 4])

subplot(3,3,6)
histogram(FinalParameters3(:,n2loc),'NumBins', BinNumber,'BinLimits',[0,4]);
%ylim([0 35])
xticks([0 1 2 3 4])

subplot(3,3,7)
x = FinalParameters1(:,n1loc);
y = FinalParameters1(:,n2loc);
scatter(x,y)
hold on
plot(0:.1:3,0:.1:3)
xlabel('n1')
ylabel('n2')
legend('Data',sprintf('n_1>n_2; %0.2f%%', 100*sum(x>y)/length(y)),'Location','best');

subplot(3,3,8)
x = FinalParameters2(:,n1loc);
y = FinalParameters2(:,n2loc);
scatter(x,y)
hold on
plot(0:.1:3,0:.1:3)
xlabel('n1')
ylabel('n2')
legend('Data',sprintf('n_1>n_2; %0.2f%%', 100*sum(x>y)/length(y)),'Location','best');

subplot(3,3,9)
x = FinalParameters3(:,n1loc);
y = FinalParameters3(:,n2loc);
scatter(x,y)
hold on
plot(0:.1:3,0:.1:3)
xlabel('n1')
ylabel('n2')
legend('Data',sprintf('n_1>n_2; %0.2f%%', 100*sum(x>y)/length(y)),'Location','best');

BinLimitSetting=[1 1 1 .5 .5 4 4];

figure
%%
for i =1:7
    subplot(3,3,i)
    h1 = histogram(FinalParameters1(:,i),'NumBins', BinNumber,'BinLimits',[0,BinLimitSetting(i)]);
    hold on
    h2 = histogram(FinalParameters2(:,i),'NumBins', BinNumber,'BinLimits',[0,BinLimitSetting(i)]);
    h3 = histogram(FinalParameters3(:,i),'NumBins', BinNumber,'BinLimits',[0,BinLimitSetting(i)]);
    title(titlesCon(i))
    %ylim([0 10])
end
legend('Original','Uniform','None')
% subplot(2,3,2)
% histogram(FinalParametersDynamic1(:,9),'NumBins', BinNumber,'BinLimits',[0,4]);
% hold on
% histogram(FinalParametersDynamic2(:,9),'NumBins', BinNumber,'BinLimits',[0,4]);
% histogram(FinalParametersDynamic3(:,9),'NumBins', BinNumber,'BinLimits',[0,4]);
% legend('mmh original objective','mmh uniform objective','mmh no objective')
% title('N2')
% ylim([0 10])
% 
% subplot(2,3,3)
% histogram(FinalParametersDynamic1(:,14),'NumBins', BinNumber*2,'BinLimits',[0,1]);
% hold on
% histogram(FinalParametersDynamic2(:,14),'NumBins', BinNumber*2,'BinLimits',[0,1]);
% histogram(FinalParametersDynamic3(:,14),'NumBins', BinNumber*2,'BinLimits',[0,1]);
% title('Total Error')
% ylim([0 20])
% 
% subplot(2,3,4)
% histogram(FinalParametersDynamic1(:,11),'NumBins', BinNumber,'BinLimits',[0,.25]);
% hold on
% histogram(FinalParametersDynamic2(:,11),'NumBins', BinNumber,'BinLimits',[0,.25]);
% histogram(FinalParametersDynamic3(:,11),'NumBins', BinNumber,'BinLimits',[0,.25]);
% title('Leading Error')
% ylim([0 20])
% 
% subplot(2,3,5)
% histogram(FinalParametersDynamic1(:,12),'NumBins', BinNumber,'BinLimits',[0,.25]);
% hold on
% histogram(FinalParametersDynamic2(:,12),'NumBins', BinNumber,'BinLimits',[0,.25]);
% histogram(FinalParametersDynamic3(:,12),'NumBins', BinNumber,'BinLimits',[0,.25]);
% title('Lagging Error')
% ylim([0 20])
% 
% subplot(2,3,6)
% histogram(FinalParametersDynamic1(:,13),'NumBins', BinNumber,'BinLimits',[0,1]);
% hold on
% histogram(FinalParametersDynamic2(:,13),'NumBins', BinNumber,'BinLimits',[0,1]);
% histogram(FinalParametersDynamic3(:,13),'NumBins', BinNumber,'BinLimits',[0,1]);
% title('Cytoplasm Error')
% ylim([0 20])
% for j=1:7
%     subplot(4,3,j);
%     plot(hits(7:end,8),hits(7:end,j))
%     title(titles(j))
%     hold on
% end
% for j=9:11
%     subplot(4,3,j-1);
%     plot(hits(7:end,8),hits(7:end,j))
%     title(titles(j))
%     hold on
% end