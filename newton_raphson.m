function [Res, p1, q1]=newton_raphson(Guess, G, B, P, Q, tol)
format long;
syms v1 v2 v3 v4 d1 d2 d3 d4;
temp=[Guess(1) Guess(2) Guess(3) Guess(4) Guess(5) Guess(6) Guess(7) Guess(8)];
%Functions definition
P2(v1, v2, v3, v4, d1, d2, d3, d4)= v2^2*G(2,2)+v2*v1*(B(2,1)*sin(d2-d1)+G(2,1)*cos(d2-d1))+...
    v2*v3*(B(2,3)*sin(d2-d3)+G(2,3)*cos(d2-d3))+ v2*v4*(B(2,4)*sin(d2-d4)+G(2,4)*cos(d2-d4));
Q2(v1, v2, v3, v4, d1, d2, d3, d4)=-v2^2*B(2,2)+v2*v1*(G(2,1)*sin(d2-d1)-B(2,1)*cos(d2-d1))+...
    v2*v3*(G(2,3)*sin(d2-d3)-B(2,3)*cos(d2-d3))+ v2*v4*(G(2,4)*sin(d2-d4)-B(2,4)*cos(d2-d4));
P3(v1, v2, v3, v4, d1, d2, d3, d4)= v3^2*G(3,3)+v3*v1*(B(3,1)*sin(d3-d1)+G(3,1)*cos(d3-d1))+...
    v3*v2*(B(3,2)*sin(d3-d2)+G(3,2)*cos(d3-d2))+ v3*v4*(B(3,4)*sin(d3-d4)+G(3,4)*cos(d3-d4));
Q3(v1, v2, v3, v4, d1, d2, d3, d4)=-v3^2*B(3,3)+v3*v1*(G(3,1)*sin(d3-d1)-B(3,1)*cos(d3-d1))+...
    v3*v2*(G(3,2)*sin(d3-d2)-B(3,2)*cos(d3-d2))+ v3*v4*(G(3,4)*sin(d3-d4)-B(3,4)*cos(d3-d4));
P4(v1, v2, v3, v4, d1, d2, d3, d4)= v4^2*G(4,4)+v4*v1*(B(4,1)*sin(d4-d1)+G(4,1)*cos(d4-d1))+...
    v4*v2*(B(4,2)*sin(d4-d2)+G(4,2)*cos(d4-d2))+ v4*v3*(B(4,3)*sin(d4-d3)+G(4,3)*cos(d4-d3));
Q4(v1, v2, v3, v4, d1, d2, d3, d4)=-v4^2*B(4,4)+v4*v1*(G(4,1)*sin(d4-d1)-B(4,1)*cos(d4-d1))+...
    v4*v2*(G(4,2)*sin(d4-d2)-B(4,2)*cos(d4-d2))+ v4*v3*(G(4,3)*sin(d4-d3)-B(4,3)*cos(d4-d3));

%Jacobian Matrix
J=jacobian([P2, P3, P4, Q2, Q3, Q4], [v2 v3 v4 d2 d3 d4]);


for i=1:1000
   fprintf('Iteration %d\n', i);
   P2g=double(P2(temp(1), temp(2), temp(3), temp(4), temp(5), temp(6), temp(7), temp(8)));
   Q2g=double(Q2(temp(1), temp(2), temp(3), temp(4), temp(5), temp(6), temp(7), temp(8)));
   P3g=double(P3(temp(1), temp(2), temp(3), temp(4), temp(5), temp(6), temp(7), temp(8)));
   Q3g=double(Q3(temp(1), temp(2), temp(3), temp(4), temp(5), temp(6), temp(7), temp(8)));
   P4g=double(P4(temp(1), temp(2), temp(3), temp(4), temp(5), temp(6), temp(7), temp(8)));
   Q4g=double(Q4(temp(1), temp(2), temp(3), temp(4), temp(5), temp(6), temp(7), temp(8)));
   b=[P(2)-P2g; P(3)-P3g; P(4)-P4g; Q(2)-Q2g; Q(3)-Q3g; Q(4)-Q4g];
   Y=double(J(temp(1), temp(2), temp(3), temp(4), temp(5), temp(6), temp(7), temp(8)));
   %x is the increment found at the current iteration
   x=Y\b;
   
   %update temp variable with increment
   temp(2)=temp(2)+x(1);
   temp(3)=temp(3)+x(2);
   temp(4)=temp(4)+x(3);
   temp(6)=temp(6)+x(4);
   temp(7)=temp(7)+x(5);
   temp(8)=temp(8)+x(6);
   if abs(x(1))<tol && abs(x(2))<tol && abs(x(3))<tol && abs(x(4))<tol && abs(x(5))<tol && abs(x(6))<tol
       break;
   else
       bias=[abs(x(1)) abs(x(2)) abs(x(3)) abs(x(4)) abs(x(5)) abs(x(6))];
       fprintf('Maximum bias = %f\n', max(bias));
   end
end
Res=temp;
P1(v1, v2, v3, v4, d1, d2, d3, d4)=v1^2*G(1,1)+v1*v2*(B(1,2)*sin(d1-d2)+G(1,2)*cos(d1-d2))+...
    v1*v3*(B(1,3)*sin(d1-d3)+G(1,3)*cos(d1-d3))+ v1*v4*(B(1,4)*sin(d1-d4)+G(1,4)*cos(d1-d4));
Q1(v1, v2, v3, v4, d1, d2, d3, d4)=-v1^2*B(1,1)+v1*v2*(G(1,2)*sin(d1-d2)-B(1,2)*cos(d1-d2))+...
    v1*v3*(G(1,3)*sin(d1-d3)-B(1,3)*cos(d1-d3))+ v1*v4*(G(1,4)*sin(d1-d4)-B(1,4)*cos(d1-d4));

p1=double(P1(temp(1), temp(2), temp(3), temp(4), temp(5), temp(6), temp(7), temp(8)));
q1=double(Q1(temp(1), temp(2), temp(3), temp(4), temp(5), temp(6), temp(7), temp(8)));

end


