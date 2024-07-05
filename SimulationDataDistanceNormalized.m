function Distance = SimulationDataDistanceNormalized(DataInt,SimInt)
%SimulationDataDistanceNormalized Takes two equal length vectors and computes their distance from each other
%   Detailed explanation goes here
DataIntLead = DataInt(:,1);
DataIntLagg = DataInt(:,2);
DataIntCyto = DataInt(:,3);

SimIntLead = SimInt(:,1);
SimIntLagg = SimInt(:,2);
SimIntCyto = SimInt(:,3);


LeadDis=sum(abs(DataIntLead-SimIntLead));
LaggDis=sum(abs(DataIntLagg-SimIntLagg));
CytoDis=sum(abs(DataIntCyto-SimIntCyto));
Distance=[LeadDis,LaggDis,CytoDis];


end

