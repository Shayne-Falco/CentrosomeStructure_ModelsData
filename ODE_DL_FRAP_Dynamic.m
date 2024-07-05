function [dydt] = ODE_DL_FRAP_Dynamic(~,Parameters,Concentrations)
%ODE_Threecompartment Global morphogen at fixed level, no limiting factor to either component

a1 = Parameters(1);
a2 = Parameters(2);
a3 = Parameters(3);
a4 = Parameters(4);
a5 = Parameters(5);
a6 = Parameters(6);
a7 = Parameters(7);
n1 = Parameters(8);
n2 = Parameters(9);
C_0 = Parameters(10);
F=Parameters(11);
B=Parameters(12);

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
dydt(3) = a3*C_0^n2*(1-(C2b*B+C2f*F)/a5)*Cyb*(C2b*B+C2f*F)^n2-a2*C2b;
dydt(4) = a3*C_0^n2*(1-(C2b*B+C2f*F)/a5)*Cyf*(C2b*B+C2f*F)^n2-a2*C2f;
dydt(5) = -a7*Cyb-dydt(1)-dydt(3);
dydt(6) = a6/(C_0*F)-a7*Cyf-dydt(2)-dydt(4);


end

