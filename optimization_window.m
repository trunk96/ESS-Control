function [pess, ptl]=optimization_window(x0, essref, ptlref, newpowers, ess_state, Tc)
     f=@(x)(100*(newpowers(1)-x(1)-ptlref)^2+50*(newpowers(2)-x(2)-ptlref)^2+25*(newpowers(3)-x(3)-ptlref)^2+15*(newpowers(4)-x(4)-ptlref)^2+5*(newpowers(5)-x(5)-ptlref)^2+...
         (ess_state+x(1)*Tc-essref)^2+(ess_state+(x(1)+x(2))*Tc-essref)^2+(ess_state+(x(1)+x(2)+x(3))*Tc-essref)^2+(ess_state+(x(1)+x(2)+x(3)+x(4))*Tc-essref)^2+(ess_state+(x(1)+x(2)+x(3)+x(4)+x(5))*Tc-essref)^2);
     %{
     A=[1 0 0 0 0;
        -1 0 0 0 0;
        1 0 0 0 0;
        -1 0 0 0 0;
        1 1 0 0 0;
        -1 -1 0 0 0;
        0 1 0 0 0;
        0 -1 0 0 0;
        1 1 1 0 0;
        -1 -1 -1 0 0;
        0 0 1 0 0;
        0 0 -1 0 0;
        1 1 1 1 0;
        -1 -1 -1 -1 0;
        0 0 0 1 0;
        0 0 0 -1 0;
        1 1 1 1 1;
        -1 -1 -1 -1 -1;
        0 0 0 0 1;
        0 0 0 0 -1];
    b=[ess_state/Tc;
       (2000000-ess_state)/Tc;
       2000000;
       2000000;
       ess_state/Tc;
       (2000000-ess_state)/Tc;
       2000000;
       2000000;
       ess_state/Tc;
       (2000000-ess_state)/Tc;
       2000000;
       2000000;
       ess_state/Tc;
       (2000000-ess_state)/Tc;
       2000000;
       2000000;
       ess_state/Tc;
       (2000000-ess_state)/Tc;
       2000000;
       2000000];
   %}
         A=[];
         b=[];
   x=fmincon(f, x0, A, b);
   pess=x(1);
   ptl=newpowers(1)-x(1);
end