clear all
close all

titlesDyn={'a1','a2','a3','a4','a5', 'a6','a7','n1','n2','t','Lead Error','Lagger Error','Cytoplasm Error','Total Error'}; %Dynamic
titlesCon={'a1','a2','a3','a4','a5','n1','n2','C_0','Cy','t','Lead Error','Lagger Error','Cytoplasm Error','Total Error'}; %Constant and Conserved

S=load('FinalParametersDynamic.mat');
FinalParametersDynamic = S.FinalParameters;
S=load('FinalParametersConstant.mat');
FinalParametersConstant = S.FinalParameters;
S=load('FinalParametersConserved.mat');
FinalParametersConserved = S.FinalParameters;

BinNumber = 16;

figure

subplot(2,3,1)
histogram(FinalParametersDynamic(:,8),'NumBins', BinNumber,'BinLimits',[0,4],'FaceColor','r');
ylabel('n1')
title('Dynamic')
ylim([0 10])
xticks([0 1 2 3 4])

subplot(2,3,2)
histogram(FinalParametersConstant(:,6),'NumBins', BinNumber,'BinLimits',[0,4],'FaceColor','r');
title('Constant')
ylim([0 10])
xticks([0 1 2 3 4])

subplot(2,3,3)
histogram(FinalParametersConserved(:,6),'NumBins', BinNumber,'BinLimits',[0,4],'FaceColor','r');
title('Conserved')
ylim([0 10])
xticks([0 1 2 3 4])

subplot(2,3,4)
histogram(FinalParametersDynamic(:,9),'NumBins', BinNumber,'BinLimits',[0,4]);
ylim([0 10])
ylabel('n2')
xticks([0 1 2 3 4])

subplot(2,3,5)
histogram(FinalParametersConstant(:,7),'NumBins', BinNumber,'BinLimits',[0,4]);
ylim([0 10])
xticks([0 1 2 3 4])

subplot(2,3,6)
histogram(FinalParametersConserved(:,7),'NumBins', BinNumber,'BinLimits',[0,4]);
ylim([0 10])
xticks([0 1 2 3 4])





figure
%%
subplot(2,3,1)
histogram(FinalParametersDynamic(:,8),'NumBins', BinNumber,'BinLimits',[0,4]);
hold on
histogram(FinalParametersConstant(:,6),'NumBins', BinNumber,'BinLimits',[0,4]);
histogram(FinalParametersConserved(:,6),'NumBins', BinNumber,'BinLimits',[0,4]);
title('N1')
ylim([0 10])

subplot(2,3,2)
histogram(FinalParametersDynamic(:,9),'NumBins', BinNumber,'BinLimits',[0,4]);
hold on
histogram(FinalParametersConstant(:,7),'NumBins', BinNumber,'BinLimits',[0,4]);
histogram(FinalParametersConserved(:,7),'NumBins', BinNumber,'BinLimits',[0,4]);
legend('Dynamic','Constant','Conservation')
title('N2')
ylim([0 10])

subplot(2,3,3)
histogram(FinalParametersDynamic(:,14),'NumBins', BinNumber*2,'BinLimits',[0,1]);
hold on
histogram(FinalParametersConstant(:,14),'NumBins', BinNumber*2,'BinLimits',[0,1]);
histogram(FinalParametersConserved(:,14),'NumBins', BinNumber*2,'BinLimits',[0,1]);
title('Total Error')
ylim([0 20])

subplot(2,3,4)
histogram(FinalParametersDynamic(:,11),'NumBins', BinNumber,'BinLimits',[0,.25]);
hold on
histogram(FinalParametersConstant(:,11),'NumBins', BinNumber,'BinLimits',[0,.25]);
histogram(FinalParametersConserved(:,11),'NumBins', BinNumber,'BinLimits',[0,.25]);
title('Leading Error')
ylim([0 20])

subplot(2,3,5)
histogram(FinalParametersDynamic(:,12),'NumBins', BinNumber,'BinLimits',[0,.25]);
hold on
histogram(FinalParametersConstant(:,12),'NumBins', BinNumber,'BinLimits',[0,.25]);
histogram(FinalParametersConserved(:,12),'NumBins', BinNumber,'BinLimits',[0,.25]);
title('Lagging Error')
ylim([0 20])

subplot(2,3,6)
histogram(FinalParametersDynamic(:,13),'NumBins', BinNumber,'BinLimits',[0,1]);
hold on
histogram(FinalParametersConstant(:,13),'NumBins', BinNumber,'BinLimits',[0,1]);
histogram(FinalParametersConserved(:,13),'NumBins', BinNumber,'BinLimits',[0,1]);
title('Cytoplasm Error')
ylim([0 20])
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