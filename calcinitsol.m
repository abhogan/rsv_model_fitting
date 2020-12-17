function calcinitsol(x)
global initcond
% take values for beta0, beta1, phi, nu, and return the solution of the ODEs
betaA0=x(1);
betaB0=x(2);
beta1=x(3);
phi=x(4);
nu=x(5);
gamma=1/1.3;
delta=1/0.57;
year=52; % weeks in a year
tend=year*72;
eta=1/52; % ageing rate - inverse of time spent in each age group (i.e. inverse of 1 year)
mu=431/1658992; % birth rate - update according to location

parameters=[betaA0 betaB0 beta1 phi gamma delta mu nu eta];

I1=0.0001;
E1=I1;
S1=0.005;
I2=0.00001;
E2=I2;
S2=0.003;
R1=0.05-S1-I1-E1;
R2=0.95-S2-I2-E2;
J1=E1;
J2=E2;

initconds=[S1 E1 I1 R1 S2 E2 I2 R2 J1 J2];

tspan=[0:1:tend];
options=[];

[t,y]=ode45(@ODEs,tspan,initconds,options,parameters);

initcond=y(end,:);
end
          