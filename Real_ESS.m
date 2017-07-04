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

A=[27  27  27  27  25  25  25  25  24  24  24  24 24  24  24  24  24  24  24  24  25  25  25  25 27  27  27  27  34  34  34  34  38  38  38  38 40  40  40  40  41  41  41  41  41  41  41  41 40  40  40  40  38  38  38  38  39  39  39  39 39  39  39  39  40  40  40  40  41  41  41  41 40  40  40  40  43  43  43  43  44  44  44  44 42  42  42  42  38  38  38  38  34  34  34  34];


P2=A.*(600000/44);
P3=A.*(620000/44);
P4=A.*(620000/44);
Q2=A.*(190000/44);
Q3=A.*(50000/44);
Q4=A.*(50000/44);

num=96;

P1=zeros(1, num);
Q1=zeros(1, num);
  
V_base=20000;
S_base=630000; 
Y_base=S_base/(V_base^2);

ess_state=1000000;
essref=1000000;
Tc=1/4;

Iter=1:1:num;
ChargeESS=zeros(1,num);
ChargeESS(1)=ess_state;
PTL=zeros(1,num);

PESS=zeros(1,num);
P_2=zeros(1, num);




for i=1:num
P=[P1(i) P2(i) P3(i) P4(i)];
Q=[Q1(i) Q2(i) Q3(i) Q4(i)];
[state_variables, p1, q1]=newton_raphson(Guess./V_base, G./Y_base, B./Y_base, P./S_base, Q./S_base, 0.000001);

P1(i)=p1*S_base;
Q1(i)=q1*S_base;


ptlref=P1(i);

x0=[state_variables(2) state_variables(3) state_variables(4) state_variables(6) state_variables(7) state_variables(8) 0 0];
newpower=[P Q];


fprintf('Iteration %d\n\n', i);
%newpower(2)=P(2)-50000*sin(i);
%newpower(2)=P(2)-5000*((1/sqrt(pi*2))*exp(-(((i-50)/2)^2)));
newpower(2)=P(2)-50000*rand();
%newpower(2)=P(2);
%newpower(2)=A(i)*(500000/44);
P_2(i)=newpower(2);
[state_variables, pess, p1, p2, p3, p4, q1, q2, q3, q4]=optimization_secondary(x0, essref/S_base, ptlref/S_base, newpower./S_base, ess_state/S_base, state_variables, Tc, G./Y_base, B./Y_base);
fprintf('\n\n\nNEW\nV1: %f\nV2: %f\nV3: %f\nV4: %f\nD1: %f\nD2: %f\nD3: %f\nD4: %f\n', state_variables(1)*V_base, state_variables(2)*V_base, state_variables(3)*V_base, state_variables(4)*V_base, state_variables(5), state_variables(6), state_variables(7), state_variables(8));

fprintf('PESS: %f\nP1: %f\nP2: %f\nP3: %f\nP4: %f\nQ1: %f\nQ2: %f\nQ3: %f\nQ4: %f\n',pess*S_base, p1*S_base, p2*S_base, p3*S_base, p4*S_base, q1*S_base, q2*S_base, q3*S_base, q4*S_base);

ess_state=ess_state+pess*S_base*Tc;
ChargeESS(i)=ess_state;
PTL(i)=p1*S_base;
PESS(i)=pess*S_base;
end;

figure

subplot(4,1,1)
plot(Iter, P_2, Iter, P2, 'black', 'LineWidth',2)
xlim([0 num])
ylim([min(P2) max(P2)])
xlabel('Time (15 minutes)', 'FontSize', 11 , 'fontweight', 'b' )
title('Power change on Node 2', 'FontSize', 11 , 'fontweight', 'b' )
ylabel('W', 'FontSize', 11 , 'fontweight', 'b' )
set(gca, 'LineWidth', 2, 'FontSize', 11, 'fontweight', 'b')
set(gcf, 'Color', 'w')
grid on

subplot(4,1,2)
plot(Iter,PESS, 'black', 'LineWidth',2)
xlim([0 num])
ylim([min(PESS) max(PESS)])
xlabel('Time (15 minutes)', 'FontSize', 11 , 'fontweight', 'b' )
title('Power from ESS', 'FontSize', 11 , 'fontweight', 'b' )
ylabel('W', 'FontSize', 11 , 'fontweight', 'b' )
set(gca, 'LineWidth', 2, 'FontSize', 11, 'fontweight', 'b')
set(gcf, 'Color', 'w')
grid on

subplot(4,1,3)
plot(Iter,ChargeESS, 'black', 'LineWidth',2)
xlim([0 num])
ylim([min(ChargeESS) max(ChargeESS)])
xlabel('Time (15 minutes)', 'FontSize', 11 , 'fontweight', 'b' )
title('ESS Charge', 'FontSize', 11 , 'fontweight', 'b' )
ylabel('Wh', 'FontSize', 11 , 'fontweight', 'b' )
set(gca, 'LineWidth', 2, 'FontSize', 11, 'fontweight', 'b')
set(gcf, 'Color', 'w')
grid on

subplot(4,1,4)
plot(Iter,PTL, 'black', 'LineWidth',2)
xlim([0 num])
ylim([min(PTL) max(PTL)])
xlabel('Time (15 minutes)', 'FontSize', 11 , 'fontweight', 'b' )
title('Power from transmission line', 'FontSize', 11 , 'fontweight', 'b' )
ylabel('W', 'FontSize', 11 , 'fontweight', 'b' )
set(gca, 'LineWidth', 2, 'FontSize', 11, 'fontweight', 'b')
set(gcf, 'Color', 'w')
grid on





  

  
  


  
  
 

 

 
 