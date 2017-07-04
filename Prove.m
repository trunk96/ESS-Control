clear all
close all
clc

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

V_base=20000;
S_base=630000; 
Y_base=S_base/(V_base^2);

[state_variables, p1, q1]=newton_raphson(Guess./V_base, G./Y_base, B./Y_base, P./S_base, Q./S_base, 0.000001);

fprintf('P1: %f\nQ1: %f\nV1: %f\nV2: %f\nV3: %f\nV4: %f\nD1: %f\nD2: %f\nD3: %f\nD4: %f\n', p1*S_base, q1*S_base, state_variables(1)*V_base, state_variables(2)*V_base, state_variables(3)*V_base, state_variables(4)*V_base, state_variables(5), state_variables(6), state_variables(7), state_variables(8));

P(1)=p1*S_base;
Q(1)=q1*S_base;
%No changes on Node 2
ess_state=1000000;
essref=1000000;
ptlref=p1*S_base;
x0=[state_variables(2) state_variables(3) state_variables(4) state_variables(6) state_variables(7) state_variables(8) 0 0];
newpower=[P Q];
Tc=1/4;


num=96;

Iter=1:1:num;
ChargeESS=zeros(1,num);
ChargeESS(1)=ess_state;
PTL=zeros(1,num);
PTL(1)=p1;
PESS=zeros(1,num);
P_2=zeros(1, num);
P_2(1)=P(2);






for i=1:num
fprintf('Iteration %d\n\n', i);
%newpower(2)=P(2)-5000*sin(i);
newpower(2)=P(2)-5000*((1/sqrt(pi*2))*exp(-(((i-50)/2)^2)));
%newpower(2)=P(2)-5000*rand();
P_2(i)=newpower(2);
[state_variables, pess, p1, p2, p3, p4, q1, q2, q3, q4]=optimization_secondary(x0, essref/S_base, ptlref/S_base, newpower./S_base, ess_state/S_base, state_variables, Tc, G./Y_base, B./Y_base);
fprintf('\n\n\nNEW\nV1: %f\nV2: %f\nV3: %f\nV4: %f\nD1: %f\nD2: %f\nD3: %f\nD4: %f\n', state_variables(1)*V_base, state_variables(2)*V_base, state_variables(3)*V_base, state_variables(4)*V_base, state_variables(5), state_variables(6), state_variables(7), state_variables(8));
fprintf('PESS: %f\nP1: %f\nP2: %f\nP3: %f\nP4: %f\nQ1: %f\nQ2: %f\nQ3: %f\nQ4: %f\n',pess*S_base, p1*S_base, p2*S_base, p3*S_base, p4*S_base, q1*S_base, q2*S_base, q3*S_base, q4*S_base);

ess_state=ess_state+pess*S_base*Tc;
ChargeESS(i)=ess_state;
PTL(i)=p1*S_base;
PESS(i)=pess*S_base;
end;




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

