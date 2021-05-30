function [TrajCalc1,TrajCalc2]=Trajectory(segmentnum)
for j=1:segmentnum;
    tt=linspace(0,pi/2,segmentnum+1);
    xx=sin(tt(j));
    yy=cos(tt(j));
    zz=sin(tt(j));
    z(:,j)=[0*xx; 1.6+0.2*yy; 0.5-0.2*zz];
    Tsd(:,:,j)=[eye(3) z(:,j);
        0 0 0 1];
end
TrajCalc1=Tsd;
TrajCalc2=z;
end