function [state_variables, pess, p1, p2, p3, p4, q1, q2, q3, q4] = optimization_secondary( x0, essref, ptlref, newpower, ess_state, prev_state, Tc, G, B)  
    
    syms v1 v2 v3 v4 d1 d2 d3 d4;
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
    Y=jacobian([P2, P3, P4, Q2, Q3, Q4], [v2 v3 v4 d2 d3 d4]);
    
    P1(v1, v2, v3, v4, d1, d2, d3, d4)=v1^2*G(1,1)+v1*v2*(B(1,2)*sin(d1-d2)+G(1,2)*cos(d1-d2))+...
    v1*v3*(B(1,3)*sin(d1-d3)+G(1,3)*cos(d1-d3))+ v1*v4*(B(1,4)*sin(d1-d4)+G(1,4)*cos(d1-d4));
    Q1(v1, v2, v3, v4, d1, d2, d3, d4)=-v1^2*B(1,1)+v1*v2*(G(1,2)*sin(d1-d2)-B(1,2)*cos(d1-d2))+...
    v1*v3*(G(1,3)*sin(d1-d3)-B(1,3)*cos(d1-d3))+ v1*v4*(G(1,4)*sin(d1-d4)-B(1,4)*cos(d1-d4));

    Y1=jacobian([P1, Q1], [v2 v3 v4 d2 d3 d4]);
    
    p2=double(P2(prev_state(1), prev_state(2), prev_state(3), prev_state(4), prev_state(5), prev_state(6), prev_state(7), prev_state(8)));
    J=double(Y(prev_state(1), prev_state(2), prev_state(3), prev_state(4), prev_state(5), prev_state(6), prev_state(7), prev_state(8)));
    J1=double(Y1(prev_state(1), prev_state(2), prev_state(3), prev_state(4), prev_state(5), prev_state(6), prev_state(7), prev_state(8)));
    
    fprintf('Power on node 2: %f\n', p2);
    fprintf('New Power on node 2: %f\n', newpower(2));
    %%f=@(x)((x(7)-ptlref)^2+((ess_state+(p2-newpower(2)+(J(1,:)*([x(1); x(2); x(3); x(4); x(5); x(6)])))*Tc)-essref)^2);
    %%%%+((ess_state-(p2-newpower(2)+(J(1,:)*([x(1); x(2); x(3); x(4); x(5); x(6)])))*Tc)-essref).^(2)
    %%%%+((ess_state+(x(8)*Tc))-essref)^2
    f=@(x)(10000000*(x(7)-ptlref)^2+((ess_state+(x(8)*Tc))-essref)^2);
    %10000000 in double-weight
    % x(1)=V2-prev_state(2) x(2)=V3-prev_state(3) x(3)=V4-prev_state(4)
    % x(4)=d2-prev_state(6) x(5)=d3-prev_state(7) x(6)=d4-prev_state(8)
    %prev_state=[v1, v2, v3, v4, d1, d2, d3, d4] is prev_state model
    Aeq=[J(1,1)  J(1,2) J(1,3)   J(1,4) J(1,5)   J(1,6) 0  -1 ;
        J(2,1)   J(2,2) J(2,3)   J(2,4) J(2,5)   J(2,6) 0  0  ;
        J(3,1)   J(3,2) J(3,3)   J(3,4) J(3,5)   J(3,6) 0  0  ;
        J(4,1)   J(4,2) J(4,3)   J(4,4) J(4,5)   J(4,6) 0  0  ;
        J(5,1)   J(5,2) J(5,3)   J(5,4) J(5,5)   J(5,6) 0  0  ;
        J(6,1)   J(6,2) J(6,3)   J(6,4) J(6,5)   J(6,6) 0  0  ;
        J1(1,1) J1(1,2) J1(1,3) J1(1,4) J1(1,5) J1(1,6) -1 0  ;
        J1(2,1) J1(2,2) J1(2,3) J1(2,4) J1(2,5) J1(2,6) 0  0  ];
        
       
    beq=[newpower(2)-p2;
        newpower(3)-double(P3(prev_state(1), prev_state(2), prev_state(3), prev_state(4), prev_state(5), prev_state(6), prev_state(7), prev_state(8)));
        newpower(4)-double(P4(prev_state(1), prev_state(2), prev_state(3), prev_state(4), prev_state(5), prev_state(6), prev_state(7), prev_state(8)));
        newpower(6)-double(Q2(prev_state(1), prev_state(2), prev_state(3), prev_state(4), prev_state(5), prev_state(6), prev_state(7), prev_state(8)));
        newpower(7)-double(Q3(prev_state(1), prev_state(2), prev_state(3), prev_state(4), prev_state(5), prev_state(6), prev_state(7), prev_state(8)));
        newpower(8)-double(Q4(prev_state(1), prev_state(2), prev_state(3), prev_state(4), prev_state(5), prev_state(6), prev_state(7), prev_state(8)));
        -double(P1(prev_state(1), prev_state(2), prev_state(3), prev_state(4), prev_state(5), prev_state(6), prev_state(7), prev_state(8)));
        newpower(5)-double(Q1(prev_state(1), prev_state(2), prev_state(3), prev_state(4), prev_state(5), prev_state(6), prev_state(7), prev_state(8)))];
    
    x=fmincon(f, x0, Aeq, beq);
    %x is the new state_variables vector (but x(7) is P1 and x(8) is Pess)
    
    state_variables=[prev_state(1) x(1)+prev_state(2) x(2)+prev_state(3) x(3)+prev_state(4) prev_state(5) x(4)+prev_state(6) x(5)+prev_state(7) x(6)+prev_state(8)];
    pess=x(8);
    p1=x(7);
    p2=double(P2(state_variables(1), state_variables(2), state_variables(3), state_variables(4), state_variables(5), state_variables(6), state_variables(7), state_variables(8)));
    p3=double(P3(state_variables(1), state_variables(2), state_variables(3), state_variables(4), state_variables(5), state_variables(6), state_variables(7), state_variables(8)));
    p4=double(P4(state_variables(1), state_variables(2), state_variables(3), state_variables(4), state_variables(5), state_variables(6), state_variables(7), state_variables(8)));
    q1=double(Q1(state_variables(1), state_variables(2), state_variables(3), state_variables(4), state_variables(5), state_variables(6), state_variables(7), state_variables(8)));
    q2=double(Q2(state_variables(1), state_variables(2), state_variables(3), state_variables(4), state_variables(5), state_variables(6), state_variables(7), state_variables(8)));
    q3=double(Q3(state_variables(1), state_variables(2), state_variables(3), state_variables(4), state_variables(5), state_variables(6), state_variables(7), state_variables(8)));
    q4=double(Q4(state_variables(1), state_variables(2), state_variables(3), state_variables(4), state_variables(5), state_variables(6), state_variables(7), state_variables(8)));
end

