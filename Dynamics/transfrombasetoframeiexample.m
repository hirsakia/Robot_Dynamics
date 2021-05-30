function M=transfrombasetoframeiexample(l1,l2)
L=l1+l2

R=sym(zeros(3,3,3));
P=sym(zeros(3,3));

R(:,:,1)=eye(3);
P(:,1)=zeros(3,1);

R(:,:,2)=eye(3);
P(:,2)=[0;l1;0];

R(:,:,3)=eye(3);
P(:,3)=[0;L;0];

for i=1:3
   M(:,:,i)=[R(:,:,i) P(:,i);zeros(1,3) 1] 
end
end