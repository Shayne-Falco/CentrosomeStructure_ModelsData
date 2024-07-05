function [dydt] = ODE_DL_3C_Dynamic(~,Parameters,Concentrations)
%ODE_Threecompartment Global morphogen at fixed level, no limiting factor to either component
%  array=[a_1 a_2 a_3 a_4 a_5 a_6 a_7 n1 n2 C_0];
dydt=zeros(3,1);
C1 = Concentrations(1);
C2 = Concentrations(2);
Cy = Concentrations(3);

a1 = Parameters(1);
a2 = Parameters(2);
a3 = Parameters(3);
a4 = Parameters(4);
a5 = Parameters(5);
a6 = Parameters(6);
a7 = Parameters(7);
n1 = Parameters(8);
n2 = Parameters(9);
C0 = Parameters(10);


dydt(1) =    C0^n1*(1-C1/a4)*Cy*C1^n1-a1*C1;
dydt(2) = a3*C0^n2*(1-C2/a5)*Cy*C2^n2-a2*C2;
dydt(3) = a6/C0 - a7*Cy - dydt(1) - dydt(2);


end

