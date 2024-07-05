function [tempDistanceSimData,tempBestFitData] = BestDataSimMatchTotal(DataIntLead,DataIntLag,SimData,lengthDataInt,lengthSimInt)
%BestDataSImMatch Returns the best fit of the sim data to the intensity data for any length data set
%   Detailed explanation goes here
tempBestFitData=[];
tempDistanceSimData=[inf,inf,inf]; %Leader,Lagger,Total
SimIntLead=SimData(1:lengthSimInt,1);
SimIntLag =SimData(1:lengthSimInt,2);

DistanceTemp=SimulationDataDistanceNormalized(DataIntLead,DataIntLag,SimIntLead,SimIntLag);
for i=1:3
    if DistanceTemp(i)<tempDistanceSimData(i)
        tempDistanceSimData(i)=DistanceTemp(i);
        tempBestFitData(:,4*i-3:4*i)=[DataIntLead,DataIntLag,SimIntLead,SimIntLag];
    end
end

end

