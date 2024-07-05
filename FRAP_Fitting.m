close all
clear all
t_halfsLagg = [];
v_maxsLagg = [];
r2Lagg = [];
t_halfsLead = [];
v_maxsLead = [];
r2Lead = [];

S=load('FRAP_Data_Lead.mat');
FRAPDataLead = S.FRAPOD142AIR1;
S=load('FRAP_Data_Lagg.mat');
FRAPDataLagg = S.FRAPOD142AIR1S1;

NEBDoffsetLead = [139,145,139,170,175,163,192,209,215,211];
NEBDoffsetLagg = [124,244,186,203,177,242,182,187,196,209];
AverageLeadRec = zeros(10,51);
AverageLaggRec = zeros(10,51);

for worm = 1:10

    
    timeLead = FRAPDataLead(1:end,(worm-1)*3+1);
    FiLead = FRAPDataLead(1:end,(worm-1)*3+2);
    timeRecLead = FRAPDataLead(2:end,(worm-1)*3+1)-FRAPDataLead(2,(worm-1)*3+1);
    FiRecLead = FRAPDataLead(2:end,(worm-1)*3+2);

    AverageLeadRec(worm,:) = interp1(timeRecLead,FiRecLead,0:50);
    
    figure(1)
    subplot(1,2,1)
    plot(timeLead-NEBDoffsetLead(worm),FiLead)
    %xlim([0 60])
    ylim([0 2.5])
    hold on
    
    %% Fit: 'Vmax-Vmax exp(-t/t_off)'.
    [xData, yData] = prepareCurveData( timeRecLead', FiRecLead );
    
    % Set up fittype and options.
    ft = fittype( 'a-a*exp(-x/b)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Robust = 'LAR';
    opts.StartPoint = [15 0.141886338627215];
    
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    coeffvals = coeffvalues(fitresult);
    % Plot fit with data.
%     figure(worm+1)
%     
%     plot(timeLagg,FiLagg)
%     hold on
%     y=coeffvals(1)-coeffvals(1)*exp(-timeLagg/coeffvals(2));
%     plot(timeLagg,y)
    t_half=-coeffvals(2)*log(.5);
%     legend({'Fi',['Vmax-Vmax*exp(-t/toff)' newline ...
%         sprintf('Vmax = %g', coeffvals(1)) newline ...
%         sprintf('koff = %g', coeffvals(2)) newline ...
%         sprintf('t half = %g', t_half) newline ...
%         ]})
    t_halfsLead(worm) = t_half;
    v_maxsLead(worm) = coeffvals(1);
    r2Lead(worm) = gof.rsquare;
    
end

for worm = 1:10
    timeLagg = FRAPDataLagg(1:end,(worm-1)*3+1);
    FiLagg = FRAPDataLagg(1:end,(worm-1)*3+2);
    timeRecLagg = FRAPDataLagg(2:end,(worm-1)*3+1)-FRAPDataLagg(2,(worm-1)*3+1);
    FiRecLagg = FRAPDataLagg(2:end,(worm-1)*3+2);
    
    AverageLaggRec(worm,:) = interp1(timeRecLagg,FiRecLagg,0:50);
    
    figure(1)
    subplot(1,2,2)
    plot(timeLagg-NEBDoffsetLagg(worm),FiLagg)
    %xlim([0 60])
    ylim([0 2.5])
    hold on
    
    %% Fit: 'Vmax-Vmax exp(-t/t_off)'.
    [xData, yData] = prepareCurveData( timeRecLagg', FiRecLagg );
    
    % Set up fittype and options.
    ft = fittype( 'a-a*exp(-x/b)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Robust = 'LAR';
    opts.StartPoint = [15 0.141886338627215];
    
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    coeffvals = coeffvalues(fitresult);
    % Plot fit with data.
%     figure(worm+1)
%     
%     plot(timeLagg,FiLagg)
%     hold on
%     y=coeffvals(1)-coeffvals(1)*exp(-timeLagg/coeffvals(2));
%     plot(timeLagg,y)
    t_half=-coeffvals(2)*log(.5);
%     legend({'Fi',['Vmax-Vmax*exp(-t/toff)' newline ...
%         sprintf('Vmax = %g', coeffvals(1)) newline ...
%         sprintf('koff = %g', coeffvals(2)) newline ...
%         sprintf('t half = %g', t_half) newline ...
%         ]})
    t_halfsLagg(worm) = t_half;
    v_maxsLagg(worm) = coeffvals(1);
    r2Lagg(worm) = gof.rsquare;
end

%% Fit: 'Vmax-Vmax exp(-t/t_off)'.
[xData, yData] = prepareCurveData( 0:50, mean(AverageLeadRec,'omitnan') );

% Set up fittype and options.
ft = fittype( 'a-a*exp(-x/b)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Robust = 'LAR';
opts.StartPoint = [15 0.141886338627215];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
coeffvals = coeffvalues(fitresult);
t_halfAvgLead = -coeffvals(2)*log(.5)
v_maxsAvgLead = coeffvals(1)

%% Fit: 'Vmax-Vmax exp(-t/t_off)'.
[xData, yData] = prepareCurveData( 0:50, mean(AverageLaggRec,'omitnan') );

% Set up fittype and options.
ft = fittype( 'a-a*exp(-x/b)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Robust = 'LAR';
opts.StartPoint = [15 0.141886338627215];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
coeffvals = coeffvalues(fitresult);
t_halfAvgLagg = -coeffvals(2)*log(.5)
v_maxsAvgLagg = coeffvals(1)

figure
plot(0:50,mean(AverageLeadRec,'omitnan'),'r','LineWidth',4)
hold on
plot(0:50,mean(AverageLaggRec,'omitnan'),'blue','LineWidth',4)
plot(0:50,AverageLeadRec,'--r','LineWidth',.5)

plot(0:50,AverageLaggRec,'--blue','LineWidth',.5)
xlabel('Time Post Photobleaching')
ylabel('Flourescence Recovery Normalized to pre-bleach and post-bleach')
legend('Leading','Lagging')

figure
plot(0:50,mean(AverageLeadRec,'omitnan'),'r','LineWidth',4)
hold on
plot(0:50,AverageLeadRec,'--r','LineWidth',.5)
xlabel('Time Post Photobleaching')
ylabel('Flourescence Recovery Normalized to pre-bleach and post-bleach')
legend('Leading')
ylim([0 2.5])

figure
plot(0:50,mean(AverageLaggRec,'omitnan'),'blue','LineWidth',4)
hold on
plot(0:50,AverageLaggRec,'--blue','LineWidth',.5)
xlabel('Time Post Photobleaching')
ylabel('Flourescence Recovery Normalized to pre-bleach and post-bleach')
legend('Lagging')
ylim([0 2.5])

figure

g1 = repmat({'Leading'},10,1);
g2 = repmat({'Lagging'},10,1);
g = [g1; g2];
boxplot([t_halfsLead'; t_halfsLagg'],g)
ylabel('t_{1/2}')

figure
boxplot([v_maxsLead'; v_maxsLagg'],g)
ylabel('V_{max}')
% ylim([0 30])
%subplot(1,2,2)
%boxplot()
%ylim([0 30])

figure
subplot(1,2,1)
scatter(-NEBDoffsetLead,t_halfsLead,'red')
hold on
scatter(-NEBDoffsetLagg,t_halfsLagg,'blue')
%legend('Leading','Lagging')
xlabel('photobleach timing relative to NEBD')
ylabel('t_{half}')

subplot(1,2,2)
scatter(-NEBDoffsetLead,v_maxsLead,'red')
hold on
scatter(-NEBDoffsetLagg,v_maxsLagg,'blue')
legend('Leading','Lagging')
xlabel('photobleach timing relative to NEBD')
ylabel('V_{max}')

figure
scatter(v_maxsLead,t_halfsLead,'red')
hold on
scatter(v_maxsLagg,t_halfsLagg,'blue')
legend('Leading','Lagging')
xlabel('V_{max}')
ylabel('t_{half}')


AverageLeadRecovery = mean(AverageLeadRec,'omitnan');
AverageLaggRecovery = mean(AverageLaggRec,'omitnan');

% mean(t_halfsLead)
% std(t_halfsLead)
% mean(v_maxsLead)
% std(v_maxsLead)
% mean(t_halfsLagg)
% std(t_halfsLagg)
% mean(v_maxsLagg)
% std(v_maxsLagg)

[h,p,ci,stats] = ttest(t_halfsLead,t_halfsLagg)

[h,p,ci,stats] = ttest(v_maxsLead,v_maxsLagg)