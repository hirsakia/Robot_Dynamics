function [Tsbn,ztestn,lastthet]=invkinNR(Tsd,segmentnum,Tsb,epsilon,Jacobianspace)

for n=1:segmentnum
    n
    clear theta Tsbv Tbsv skew AdTbs Js Jb Vb vb skew  argu
    if (n==1)
        theta(:,1)=[0 0 0 0 0 0]';
    else
        theta(:,1)=lastthet(:,n-1);
    end
    
    for i=1:100
i
        %%vpa(subs(Tsb, {sym('theta1'), sym('theta2'), sym('theta3'), sym('theta4'), sym('theta5'),....
        %%sym('theta6'), sym('l1'), sym('l2'), sym('l3'),...
        %%sym('l4')}, {rem(real(lastthet(1,2)),2*pi) rem(real(lastthet(2,2)),2*pi) rem(real(lastthet(3,2)),...
        %%2*pi) rem(real(lastthet(4,2)),2*pi) rem(real(lastthet(5,2)),2*pi) rem(real(lastthet(6,2)),2*pi)    0.6,0.7,0.5,0.5}))
        
        Tsbv=vpa(subs(Tsb, {sym('theta1'), sym('theta2'), sym('theta3'), sym('theta4'), sym('theta5'),....
            sym('theta6'), sym('l1'), sym('l2'), sym('l3'), sym('l4')},...
            {theta(1,i) theta(2,i) theta(3,i) theta(4,i) theta(5,i) theta(6,i)    0.6,0.7,0.5,0.5}));
        
        Tbsv=pinv(Tsbv);
        
        skew=[0   -Tbsv(3,4) Tbsv(2,4);
            Tbsv(3,4)   0  -Tbsv(1,4);
            -Tbsv(2,4) Tbsv(1,4)  0];
        
        AdTbs=[Tbsv(1:3,1:3) zeros(3,3)
            skew*Tbsv(1:3,1:3) Tbsv(1:3,1:3)];
        
        Js=vpa(subs(Jacobianspace , {sym('theta1'), sym('theta2'), sym('theta3'), sym('theta4'), sym('theta5'),....
            sym('theta6'), sym('l1'), sym('l2'), sym('l3'), sym('l4')},...
            {theta(1,i) theta(2,i) theta(3,i) theta(4,i) theta(5,i) theta(6,i)   0.6,0.7,0.5,0.5}));
        
        Jb=(AdTbs)*Js;
        
        argu=Tbsv*Tsd(:,:,n);
        
        Vb=logm(argu);
        
        vb=[-Vb(2,3) Vb(1,3) -Vb(1,2) Vb(1:3,4)']';
        
        
        if (norm(vb(4:end)) > epsilon)%% && (norm(vb(1:3)) > epsilon)
            theta(:,i+1)=theta(:,i)+(pinv(Jb))*vb;
        else
            lastthet(:,n)=theta(:,end);
            break
        end
    end
    
    Tsbn(:,:,n)=vpa(subs(Tsb, {sym('theta1'), sym('theta2'), sym('theta3'), sym('theta4'), sym('theta5'),....
            sym('theta6'), sym('l1'), sym('l2'), sym('l3'), sym('l4')},...
            {lastthet(1,n) lastthet(2,n) lastthet(3,n) lastthet(4,n) lastthet(5,n) lastthet(6,n)    0.6,0.7,0.5,0.5}));
    ztestn(:,n)=Tsbn(1:3,4,n);
end

end