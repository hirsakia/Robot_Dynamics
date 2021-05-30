function matefunction = mate(w2,v2,h)

syms theta1 theta2 theta3 theta4 theta5 theta6
theta=[theta1;theta2;theta3;theta4;theta5;theta6];

ww=[0   -w2(3) w2(2);
    w2(3)   0  -w2(1); 
    -w2(2) w2(1)  0];

matew=eye(3)+sin(theta(h))*ww+(1-cos(theta(h)))*(ww*ww);
rr=eye(3)*theta(h)+(1-cos(theta(h)))*ww+(theta(h)-sin(theta(h)))*(ww^2);
rightcolumn=(rr)*[v2(1);v2(2);v2(3)];
matefunction=[matew , rightcolumn ; 0 0 0 , 1];

end