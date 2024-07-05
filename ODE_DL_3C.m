function [dydt] = ODE_DL_3C(~,Parameters,Concentrations)
%ODE_Threecompartment Global morphogen at fixed level, no limiting factor to either component
dydt=zeros(2,1);
C1 = Concentrations(1);
C2 = Concentrations(2);


a1 = Parameters(1);
a2 = Parameters(2);
a3 = Parameters(3);
a4 = Parameters(4);
a5 = Parameters(5);

n1 = Parameters(6);
n2 = Parameters(7);

C0 = Parameters(8);

if Parameters(9)
    C_y=Parameters(9);
else
    C_y=C0-C1-C2;
end

dydt(1) = C0^n1*(1-C1/a4)*C_y*C1^n1-a1*C1;
dydt(2) = a3*C0^n2*(1-C2/a5)*C_y*C2^n2-a2*C2;

end

