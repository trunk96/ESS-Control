function [pess, ptl] = optimization(x0, essref, ptlref, newpower, ess_state, Tc)
%Optimize function for ESS Power Flow (Phase I)
%x0 starting point, essref reference charge of ESS, ptlref reference active
%power from/to trasmission line, newpower new transmission active power found with
%Newton-Raphson, ess_state current energy state of ESS, Tc sample time (in
%hours, usually 1/4 of hour)
    f=@(x)((newpower-x(2)-ptlref)^(2)+(ess_state+x(2)*Tc-essref)^2);
    A=[0 1;
       0 Tc];
    b=[2000000; 2000000-ess_state];
    Aeq=[1 1];
    beq=newpower;
    x=fmincon(f, x0, A, b, Aeq, beq);
    pess=x(2);
    ptl=x(1);
end

