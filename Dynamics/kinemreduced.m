clc
clear all

syms theta1 theta2 theta3 theta4 theta5 theta6 l1 l2 l3 l4

theta=[theta1;theta2;theta3;theta4;theta5;theta6];

L=l1+l2+l3;
d=l4;
M=[1 0 0 0;
    0 1 0 L;
    0 0 1 d;
    0 0 0 1]; %% home configuration


%% joint6
w(:,6)=[0;1;0];
v(:,6)=[-d;0;0];


%% joint5
w(:,5)=[1;0;0];
v(:,5)=[0;d;-L];


%% joint4
w(:,4)=[0;1;0];
v(:,4)=[-d;0;0];


%% joint3
w(:,3)=[1;0;0];
v(:,3)=[0;d;-l1];


%% joint2
w(:,2)=[1;0;0];
v(:,2)=[0;d;0];


%% joint1
w(:,1)=[0;0;1];
v(:,1)=[0;0;0];

for h=1:6
    w2=w(:,h);
    v2=v(:,h);
    mates(:,:,h)=mate(w2,v2,h);
end

Tsb=mates(:,:,1)*mates(:,:,2)*mates(:,:,3)*mates(:,:,4)*mates(:,:,5)*mates(:,:,6)*M;

vpa(subs(Tsb, {sym('theta1'), sym('theta2'), sym('theta3'), sym('theta4'), sym('theta5'),....
    sym('theta6'), sym('l1'), sym('l2'), sym('l3'), sym('l4')}, {0, 0, 0, 0,0,0    0.6,0.7,0.5,0.5}));

%% Space Jacobian

S=[w;v];
Jacobianspace=sym(zeros(6,6));

for h=2:6
    Adts(:,:,h)=Adt(mates,h);
    Jacobianspace(:,1)=S(:,1);
    Jacobianspace(:,h)=Adts(:,:,h)*S(:,h);
end

% JS=vpa(subs(Jacobianspace , {sym('theta1'), sym('theta2'), sym('theta3'), sym('theta4'), sym('theta5'),....
% sym('theta6'), sym('l1'), sym('l2'), sym('l3'), sym('l4')}, {0, pi/2, 0, 0,0,0      0.6,0.7,0.5,0.5}))

%% Inverse Kinematic

epsilon=1e-4;
lastthet(:,:)=zeros(6,1);
    
%% Desired Trajectory
segmentnum=100;
[Tsd,z]=Trajectory(segmentnum);
[Tsbn,ztestn,lastthet]=invkinNR(Tsd,segmentnum,Tsb,epsilon,Jacobianspace);

degrees=rem(real(lastthet),pi);

n=[1:segmentnum];

figure(1)
plot(n,degrees(1,:))
hold on
plot(n,degrees(2,:))
plot(n,degrees(3,:))
plot(n,degrees(4,:))
plot(n,degrees(5,:))
plot(n,degrees(6,:))
legend('theta1','theta2','theta3','theta4','theta5','theta6')
hold off

figure (2)
plot(n,real(ztestn(1,:)))
hold on
plot(n,z(1,:))
legend('xcalculated','x')
hold off

figure (3)
plot(n,real(ztestn(2,:)))
hold on
plot(n,z(2,:))
legend('ycalculated','y')
hold off

figure (4)
plot(n,real(ztestn(3,:)))
hold on
plot(n,z(3,:))
legend('zcalculated','z')
hold off
