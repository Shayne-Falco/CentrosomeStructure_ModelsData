function [dydt] = ODE_DL_FRAP_Dynamic_Symmetric(~,Parameters,Concentrations)
%ODE_Threecompartment Global morphogen at fixed level, no limiting factor to either component

a1 = Parameters(1);
a4 = Parameters(2);
a6 = Parameters(3);
a7 = Parameters(4);
n1 = Parameters(5);
C_0 = Parameters(6);
F=Parameters(7);
B=Parameters(8);

dydt=zeros(6,1);
% C1b C1f C2b C2f Cyb Cyf
C1b=Concentrations(1);
C1f=Concentrations(2);
C2b=Concentrations(3);
C2f=Concentrations(4);
Cyb=Concentrations(5);
Cyf=Concentrations(6);

% C1b C1f C2b C2f Cyb Cyf
dydt(1) =    C_0^n1*(1-(C1b*B+C1f*F)/a4)*Cyb*(C1b*B+C1f*F)^n1-a1*C1b;
dydt(2) =    C_0^n1*(1-(C1b*B+C1f*F)/a4)*Cyf*(C1b*B+C1f*F)^n1-a1*C1f;
dydt(3) =    C_0^n1*(1-(C2b*B+C2f*F)/a4)*Cyb*(C2b*B+C2f*F)^n1-a1*C2b;
dydt(4) =    C_0^n1*(1-(C2b*B+C2f*F)/a4)*Cyf*(C2b*B+C2f*F)^n1-a1*C2f;
dydt(5) = -a7*Cyb-dydt(1)-dydt(3);
dydt(6) = a6/(C_0*F)-a7*Cyf-dydt(2)-dydt(4);


end

