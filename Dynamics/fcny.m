function y = fcny(thetadotsim,thetasim,t,u,Q,C,g)

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
% 

load('Uas')
load('t')
load('thetadotsim')
load('thetasim')

Q=subs(Q,{sym('theta1'), sym('theta2'), sym('theta3'), sym('theta4'), sym('theta5'),....
            sym('theta6')},sym('thetadot1'), sym('thetadot2'), sym('thetadot3'),...
            sym('thetadot4'), sym('thetadot5'),sym('thetadot6'),{thetasim(1) thetasim(2)...
            thetasim(3) thetasim(4) thetasim(5) thetasim(6) thetadotsim(1) thetadotsim(2) thetadotsim(3)...
            thetadotsim(4) thetadotsim(5) thetadotsim(6)})
C=subs(C,{sym('theta1'), sym('theta2'), sym('theta3'), sym('theta4'), sym('theta5'),....
            sym('theta6')},sym('thetadot1'), sym('thetadot2'), sym('thetadot3'),...
            sym('thetadot4'), sym('thetadot5'),sym('thetadot6'),{thetasim(1) thetasim(2)...
            thetasim(3) thetasim(4) thetasim(5) thetasim(6) thetadotsim(1) thetadotsim(2) thetadotsim(3)...
            thetadotsim(4) thetadotsim(5) thetadotsim(6)})
g=subs(g,{sym('theta1'), sym('theta2'), sym('theta3'), sym('theta4'), sym('theta5'),....
            sym('theta6')},sym('thetadot1'), sym('thetadot2'), sym('thetadot3'),...
            sym('thetadot4'), sym('thetadot5'),sym('thetadot6'),{thetasim(1) thetasim(2)...
            thetasim(3) thetasim(4) thetasim(5) thetasim(6) thetadotsim(1) thetadotsim(2) thetadotsim(3)...
            thetadotsim(4) thetadotsim(5) thetadotsim(6)})
u = fcnu(Ueq,Uas,thetasim,thetadotsim,Q)

dx=zeros(12,1)
dx(1)=x(2)
dx(2)=inv(Q)*(u-C-g)
y=[dx(1);dx(2)]

save('myFile.mat','y', '-append')

