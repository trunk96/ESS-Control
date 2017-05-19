%Variables definition
%B =1.0e-03 *[-0.1948   +0.1948         0         0;...
%             +0.1948    -0.37698   +0.1822         0;...
%             0   +0.1822    -0.3644   +0.1822;...
%            0         0   +0.1822    -0.1822];

%G=1.0e-07 *[0.0900   -0.0900         0         0;...
%            -0.0900    0.2309   -0.1400         0;...
%            0   -0.1400    0.2736   -0.1400;...
%            0         0   -0.1400    0.1400];
 
G=[2.2566    -2.2566    0    0;
   -2.2566    4.2208    -1.9643    0;
   0    -1.9643    3.9286    -1.9643;
   0    0    -1.9643    1.9643];

B=[-2    2    0    0;
    2    -2.9529    0.9535    0;
    0    0.9535    -1.907    0.9535;
    0    0    0.9535    -0.9535];

Guess=[20000   20000    20000    20000         0           0           0          0];
P=[0 600000 620000 620000];
Q=[0 190000 50000 50000];
%First Newton_Raphson iteration
[state_variables, p1, q1]=newton_raphson(Guess, G, B, P, Q, 0.00001);
%Now simulation of power demand change on node 2
P_new=P;
ess_state=1000000; % 50% of charge initially
Iter=1:1:97;
ChargeESS=zeros(0,97);
ChargeESS(1)=ess_state;
PTL=zeros(0,97);
PTL(1)=p1;
P_else=zeros(0,97);
P_else(1)=p1;
PESS=zeros(0,97);
P_2=zeros(0, 97);
P_2(1)=P(2);
for i=2:92
    %P_new(2)=rand()*630000;
    %P_new(2)=P(2)+300000*sin(i);
    P_new(2)=P(2)-9000000*((1/sqrt(pi*2))*exp(-(((i-50)/2)^2))); %Dirac delta function
    P_2(i)=P_new(2);
    fprintf('Iteration %d, new P is %f\n',i, P_new(2));
    [state_variables, p1_new, q1_new]=newton_raphson(Guess, G, B, P_new, Q, 0.00001);
    P_else(i)=p1_new;
    [pess, ptl]=optimization([0 0], 1000000, p1, p1_new, ess_state, 1/4);
    ess_state=ess_state+pess*(1/4);
    ChargeESS(i)=ess_state;
    PTL(i)=ptl;
    PESS(i)=pess;
end

%figure=plot(Iter, PTL, Iter, PESS,'--', Iter, ChargeESS, ':');
subplot(4,1,1);
plot(Iter, PTL, Iter, P_else);
title('Power from Transmission Line');
ylabel('W');
xlabel('Time (15 minutes)');
subplot(4,1,2);
plot(Iter, PESS);
title('Power from ESS');
ylabel('W');
xlabel('Time (15 minutes)');
subplot(4,1,3);
plot(Iter, ChargeESS);
title('ESS Charge');
ylabel('Wh');
xlabel('Time (15 minutes)');
subplot(4, 1, 4);
plot(Iter, P_2);
title('Power change on Node 2');
ylabel('W');
xlabel('Time (15 minutes)');