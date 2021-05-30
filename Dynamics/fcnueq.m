function Ueq = fcnueq(e,edot,xddd,despar,landa,SS,theta,thetadot,Q,C,g)


p1=despar(1)*((abs(e(1)))^despar(5))*sign(e(1))+despar(2)*e(1)+despar(3)*despar(6)*(abs(e(1)))^(despar(6)-1)*edot(1);
p2=despar(1)*((abs(e(2)))^despar(5))*sign(e(2))+despar(2)*e(2)+despar(3)*despar(6)*(abs(e(2)))^(despar(6)-1)*edot(2);
p3=despar(1)*((abs(e(3)))^despar(5))*sign(e(3))+despar(2)*e(3)+despar(3)*despar(6)*(abs(e(3)))^(despar(6)-1)*edot(3);
p4=despar(1)*((abs(e(4)))^despar(5))*sign(e(4))+despar(2)*e(4)+despar(3)*despar(6)*(abs(e(4)))^(despar(6)-1)*edot(4);
p5=despar(1)*((abs(e(5)))^despar(5))*sign(e(5))+despar(2)*e(5)+despar(3)*despar(6)*(abs(e(5)))^(despar(6)-1)*edot(5);
p6=despar(1)*((abs(e(6)))^despar(5))*sign(e(6))+despar(2)*e(6)+despar(3)*despar(6)*(abs(e(6)))^(despar(6)-1)*edot(6);

part=[p1 p2 p3 p4 p5 p6]';

% theta1=theta(1)
% theta2=theta(2)
% theta3=theta(3)
% theta4=theta(4)
% theta5=theta(5)
% theta6=theta(6)
% 
% thetadot1=thetadot(1)
% thetadot2=thetadot(2)
% thetadot3=thetadot(3)
% thetadot4=thetadot(4)
% thetadot5=thetadot(5)
% thetadot6=thetadot(6)
load('t')
load('thetadotsim')
load('thetasim')
Q=subs(Q,{sym('theta1'), sym('theta2'), sym('theta3'), sym('theta4'), sym('theta5'),....
            sym('theta6')},sym('thetadot1'), sym('thetadot2'), sym('thetadot3'),...
            sym('thetadot4'), sym('thetadot5'),sym('thetadot6'),{theta(1) theta(2)...
            theta(3) theta(4) theta(5) theta(6) thetadot(1) thetadot(2) thetadot(3)...
            thetadot(4) thetadot(5) thetadot(6)})
C=subs(C,{sym('theta1'), sym('theta2'), sym('theta3'), sym('theta4'), sym('theta5'),....
            sym('theta6')},sym('thetadot1'), sym('thetadot2'), sym('thetadot3'),...
            sym('thetadot4'), sym('thetadot5'),sym('thetadot6'),{theta(1) theta(2)...
            theta(3) theta(4) theta(5) theta(6) thetadot(1) thetadot(2) thetadot(3)...
            thetadot(4) thetadot(5) thetadot(6)})
g=subs(g,{sym('theta1'), sym('theta2'), sym('theta3'), sym('theta4'), sym('theta5'),....
            sym('theta6')},sym('thetadot1'), sym('thetadot2'), sym('thetadot3'),...
            sym('thetadot4'), sym('thetadot5'),sym('thetadot6'),{theta(1) theta(2)...
            theta(3) theta(4) theta(5) theta(6) thetadot(1) thetadot(2) thetadot(3)...
            thetadot(4) thetadot(5) thetadot(6)})
        
Ueq=-Q*(C+g)+landa*SS-xddd+part
save('Ueq.mat','Ueq', '-append')

