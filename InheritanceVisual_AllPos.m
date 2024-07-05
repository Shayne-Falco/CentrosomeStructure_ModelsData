clear all
close all
dU = 0.005;
unknownls = dU:dU:1-dU;
length_unknownls = length(unknownls);
PossibleL1s = zeros(1,length_unknownls);
PossibleL2s = zeros(1,length_unknownls);
PossibleL3s = zeros(1,length_unknownls);
PossibleL4s = zeros(1,length_unknownls);
PossibleL5s = zeros(1,length_unknownls);

% Proba1Target = .32;
% Proba2Target = .26;
% Proba3Target = .2;
% Proba4Target = .22;
Proba1Target = .25;
Proba2Target = .25;
Proba3Target = .25;
Proba4Target = .25;
error = .00001;
Passed = [];
%for l1i = 1:length_unknownls
    %l1 = unknownls(l1i);
    %l1c = 1-l1;
    l1 = 0.5;
    l1c = 0.5;
    for l2i = 1:length_unknownls
        l2 = unknownls(l2i);
        l2c = 1-l2;
        for l3i = 1:length_unknownls
            l3 = unknownls(l3i);
            l3c = 1-l3;
            for l4i = 1:length_unknownls
                l4 = unknownls(l4i);
                l4c = 1-l4;
                for l5i = 1:length_unknownls
                    l5 = unknownls(l5i);
                    l5c = 1-l5;

                    if ~(abs(l1*l2*l3+l1c*l4*l5 - Proba1Target)<error)
                        continue
                    end
                    if ~(abs(l1*l2*l3c+l1c*l4*l5c - Proba2Target)<error)
                        continue
                    end
                    if ~(abs(l1*l2c*l3+l1c*l4c*l5 - Proba3Target)<error)
                        continue
                    end
                    if ~(abs(l1*l2c*l3c+l1c*l4c*l5c - Proba4Target)<error)
                        continue
                    end
                    Passed = [Passed;l1 l2 l3 l4 l5];
                    %PossibleL1s(l1i)=PossibleL1s(l1i)+1;
                    PossibleL2s(l2i)=PossibleL2s(l2i)+1;
                    PossibleL3s(l3i)=PossibleL3s(l3i)+1;
                    PossibleL4s(l4i)=PossibleL4s(l4i)+1;
                    PossibleL5s(l5i)=PossibleL5s(l5i)+1;
                end
            end
        end
    end
%end

subplot(3,2,1)
plot(unknownls,PossibleL1s)
ylim([0 100])
xlim([0 1])
title('lamda_1')

subplot(3,2,2)
plot(unknownls,PossibleL2s)
ylim([0 100])
xlim([0 1])
title('lamda_2')

subplot(3,2,3)
plot(unknownls,PossibleL3s)
ylim([0 100])
xlim([0 1])
title('lamda_3')

subplot(3,2,4)
plot(unknownls,PossibleL4s)
ylim([0 100])
xlim([0 1])
title('lamda_4')

subplot(3,2,5)
plot(unknownls,PossibleL5s)
ylim([0 100])
xlim([0 1])
title('lamda_5')
