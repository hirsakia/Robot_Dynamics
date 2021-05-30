clc
clear all

n=6;
gravity=9.81;
l1=0.7;
l2=0.6;
l3=0.5;
l4=0.5;
L=l1+l2+l3;
d=l4;

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

S=[w;v];

R=sym(zeros(3,3,n));
P=sym(zeros(3,n));

R(:,:,1)=eye(3);
R(:,:,2)=[0 0 1;1 0 0;0 1 0];
R(:,:,3)=R(:,:,2);
R(:,:,4)=[0 1 0;0 0 1;1 0 0];
R(:,:,5)=R(:,:,3);
R(:,:,6)=[0 1 0;1 0 0;0 0 1];

P(:,1)=[0 0 0]';
P(:,2)=[0 0 d]';
P(:,3)=[0 l1 d]';
P(:,4)=[0 l1+l2 d]';
P(:,5)=[0 L d]';
P(:,6)=P(:,5);

M=sym(zeros(4,4,n));
for i=1:n
    M(:,:,i)=[R(:,:,i) P(:,i);zeros(1,3) 1];
end

for i=1:n
    Adm(:,:,i)=AdM(M(:,:,i));
    A(:,i)=Adm(:,:,i)*S(:,i);
end

T01=M(:,:,1)*matexp(A,1);
T12=(inv(M(:,:,1))*M(:,:,2))*matexp(A,2);
T23=(inv(M(:,:,2))*M(:,:,3))*matexp(A,3);
T34=(inv(M(:,:,3))*M(:,:,4))*matexp(A,4);
T45=(inv(M(:,:,4))*M(:,:,5))*matexp(A,5);
T56=(inv(M(:,:,5))*M(:,:,6))*matexp(A,6);
% T67=(inv(M(:,:,6))*M(:,:,7))*matexp(A,7)

% T76=matexp(-A,7)*(inv(M(:,:,7))*M(:,:,6))
T65=matexp(-A,6)*(inv(M(:,:,6))*M(:,:,5));
T54=matexp(-A,5)*(inv(M(:,:,5))*M(:,:,4));
T43=matexp(-A,4)*(inv(M(:,:,4))*M(:,:,3));
T32=matexp(-A,3)*(inv(M(:,:,3))*M(:,:,2));
T21=matexp(-A,2)*(inv(M(:,:,2))*M(:,:,1));
T10=matexp(-A,1)*M(:,:,1);

I = sym(zeros(3,3,n));
% II=sym('II',[n 3])
II=2*ones(6,3)
for i=1:n
    I(:,:,i)=[II(i,1) 0 0; 0 II(i,2) 0;0 0 II(i,3)];
end

% m = sym('m', [n 1]);
m=ones(n,1)
GG=sym(zeros(6*n,6*n));
for i=1:n
    G(:,:,i)=[I(:,:,i) zeros(3,3);zeros(3,3) m(i)*eye(3,3)];
    GG(6*i-5:6*i,6*i-5:6*i)=G(:,:,i);
end

AA=sym(zeros(6*n,n));
for i=1:n
    AA(6*i-5:6*i,i)=A(:,i);
end


%% closed dynamics
v0=zeros(6,1);
vdot0=[zeros(3,1);0;-gravity;0];

adVtotal=sym(zeros(6*n,6*n));
thetadot=sym('thetadot',[n 1],'real');
V = sym('V',[6 n],'real');

V(:,1)=Adt(T10)*v0+A(:,1)*thetadot(1);
V(:,2)=Adt(T21)*V(:,1)+A(:,2)*thetadot(2);
V(:,3)=Adt(T32)*V(:,2)+A(:,3)*thetadot(3);
V(:,4)=Adt(T43)*V(:,3)+A(:,4)*thetadot(4);
V(:,5)=Adt(T54)*V(:,4)+A(:,5)*thetadot(5);
V(:,6)=Adt(T65)*V(:,5)+A(:,6)*thetadot(6);

for i=1:n
    adV(:,:,i)=advi(V,i);
    adVtotal(6*i-5:6*i,6*i-5:6*i)=adV(:,:,i);
    
    
    adthet(:,i)=A(:,i)*thetadot(i);
    adTHET(:,:,i)=advi(adthet,i);
    
end

adTHETtotal=sym(zeros(6*n,6*n));

for i=1:6
    adTHETtotal(6*i-5:6*i,6*i-5:6*i)=adTHET(:,:,i);
end

W=sym(zeros(6*n,6*n));


Pji(:,1)=T10(1:3,4);
AdTji(:,:,1)=[T10(1:3,1:3) zeros(3,3);
        skewop(Pji(:,1))*T10(1:3,1:3) T10(1:3,1:3)];
Pji(:,2)=T21(1:3,4);
AdTji(:,:,2)=[T21(1:3,1:3) zeros(3,3);
        skewop(Pji(:,2))*T21(1:3,1:3) T21(1:3,1:3)];
Pji(:,3)=T32(1:3,4);
AdTji(:,:,3)=[T32(1:3,1:3) zeros(3,3);
        skewop(Pji(:,3))*T32(1:3,1:3) T32(1:3,1:3)];
Pji(:,4)=T43(1:3,4);
AdTji(:,:,4)=[T43(1:3,1:3) zeros(3,3);
        skewop(Pji(:,4))*T43(1:3,1:3) T43(1:3,1:3)];
Pji(:,5)=T54(1:3,4);
AdTji(:,:,5)=[T54(1:3,1:3) zeros(3,3);
        skewop(Pji(:,5))*T54(1:3,1:3) T54(1:3,1:3)];
Pji(:,6)=T65(1:3,4);
AdTji(:,:,6)=[T65(1:3,1:3) zeros(3,3);
        skewop(Pji(:,6))*T65(1:3,1:3) T65(1:3,1:3)];


for i=1:5
    W(6*i+1:6*i+6,6*i-5:6*i)=AdTji(:,:,i+1);
end

Pji=T10(1:3,4);
Pji=T10(1:3,4);
AdT10=[T10(1:3,1:3) zeros(3,3);
    skewop(Pji)*T10(1:3,1:3) T10(1:3,1:3)];

vbase=[AdT10*v0;zeros(6*(n-1),1)];
vdotbase=[AdT10*vdot0;zeros(6*(n-1),1)];
% Ftip=[zeros(6*5,1);AdTji'*F]

LL=inv(eye(size(W))-W);
W=simplify(W);
LL=simplify(LL);
AA=simplify(AA);
adTHETtotal=simplify(adTHETtotal);
adVtotal=simplify(adVtotal);

MM=AA'*LL'*GG*LL*AA;
C=-AA'*LL'*(GG*LL*adTHETtotal*W+adVtotal'*GG)*LL*AA*thetadot;
g=AA'*LL'*GG*LL*vdotbase;

MM=simplify(vpa(MM))
% C=simplify(vpa(C))
g=simplify(vpa(g))
Q=inv(MM);























