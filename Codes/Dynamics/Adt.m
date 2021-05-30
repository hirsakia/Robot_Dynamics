function Adtcalc=Adt(mates)

%     if h==2
%         T=mates(:,:,1);
%     elseif h==3
%         T=mates(:,:,1)*mates(:,:,2);
%     elseif h==4
%         T=mates(:,:,1)*mates(:,:,2)*mates(:,:,3);
%     elseif h==5
%         T=mates(:,:,1)*mates(:,:,2)*mates(:,:,3)*mates(:,:,4);
%     elseif h==6
%         T=mates(:,:,1)*mates(:,:,2)*mates(:,:,3)*mates(:,:,4)*mates(:,:,5);
%     end
% 
T=mates

skew=[   0   -T(3,4) T(2,4);
      T(3,4)   0  -T(1,4);
     -T(2,4) T(1,4)  0];
Adtcalc=[T(1:3,1:3) zeros(3,3)
     skew*T(1:3,1:3) T(1:3,1:3)];

end