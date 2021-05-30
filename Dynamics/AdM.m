function Adtcalc=AdM(Mad,i)
M=inv(Mad);

skew=[   0   -M(3,4) M(2,4);
      M(3,4)   0  -M(1,4);
     -M(2,4) M(1,4)  0];
Adtcalc=[M(1:3,1:3) zeros(3,3)
     skew*M(1:3,1:3) M(1:3,1:3)];

end