function fitnessValue=GA_init(x)

landa=x(1);
nu=x(2);
gamma=x(3);

vu=x(4);
g0=x(5);
g1=x(6);
g2=x(7);
g3=x(8);
epstar=x(9)
S0=x(10);
kappa=x(11);

beta=2*x(3)/(1+x(3));


Simulation_Time=20;
a=sim('SimulationANTSMC.slx','SrcWorkspace','Current');
error = a.get('error');
control_signal=a.get('control_signal');

error_plus=sum(error.^2);
% fitnessValue=0.5*(sum([0.1*error_plus(1:3) 100*error(4:6)]+sum(0.0001*control_signal.^2)));
fitnessValue=log(sum(error_plus));




