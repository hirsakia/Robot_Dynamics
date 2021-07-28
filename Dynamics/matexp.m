function matefunction = matexp(A,i)

syms theta1 theta2 theta3 theta4 theta5 theta6;
theta1=sym('theta1','real');
theta2=sym('theta2','real');
theta3=sym('theta3','real');
theta4=sym('theta4','real');
theta5=sym('theta5','real');
theta6=sym('theta6','real');

theta=[theta1;theta2;theta3;theta4;theta5;theta6];

w2=A(1:3,i)';
v2=A(4:end,i)';

ww=[0   -w2(3) w2(2);
    w2(3)   0  -w2(1); 
    -w2(2) w2(1)  0];

matew=eye(3)+sin(theta(i))*ww+(1-cos(theta(i)))*(ww*ww);
rr=eye(3)*theta(i)+(1-cos(theta(i)))*ww+(theta(i)-sin(theta(i)))*(ww^2);
rightcolumn=(rr)*[v2(1);v2(2);v2(3)];
matefunction=[matew , rightcolumn ; 0 0 0 , 1];

end