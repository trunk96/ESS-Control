clear all 
close all
clc
%Variables definition
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
P_new=zeros(97, 4);
P_new(1, :)=P;
ess_state=1000000; % 50% of charge initially
Iter=1:1:97;
ChargeESS=zeros(0,97);
ChargeESS(1)=ess_state;
PTL=zeros(0,97);
PTL(1)=p1;
P_else=zeros(0,97);
P_else(1)=p1;
PESS=zeros(0,97);
%P_2=zeros(0, 97);
%P_2(1)=P(2);
P_1=zeros(0,97);
Q_1=zeros(0,97);
for i=2:97
    %P_new(i, 2)=P(2)+300000*sin(i);
    P_new(i, 2)=P(2)-300000*((1/sqrt(pi*2))*exp(-(((i-50)/2)^2))); %Dirac delta function
    fprintf('Iteration %d, new P is %f\n',i, P_new(i, 2));
    [state_variables, P_1(i), Q_1(i)]=newton_raphson(Guess, G, B, P_new(i,:), Q, 0.000001);
end
for i=2:92
    [pess, ptl]=optimization_window([0 0 0 0 0], 1000000, p1, P_1(i:i+4), ess_state, 1/4);
    ess_state=ess_state+pess*(1/4);
    ChargeESS(i)=ess_state;
    PTL(i)=ptl;
    PESS(i)=pess;
end
for i=93:97
    PTL(i)=PTL(92)
    PESS(i)=PESS(92);
    ChargeESS(i)=ChargeESS(92);
end
P_2=P_new(:, 2).';
subplot(4,1,1);
plot(Iter, PTL);
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