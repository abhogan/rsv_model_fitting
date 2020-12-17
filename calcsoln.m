function [t,y]=calcsoln(x)
global initcond
% take values for beta0, beta1, phi, nu, and return the solution of the ODEs
betaA0=x(1);
betaB0=x(2);
beta1=x(3);
phi=x(4);
nu=x(5);
gamma=1/1.3;
delta=1/0.57;
year=52;
tend=year*72;
eta=1/52;
mu=431/1658992;

parameters=[betaA0 betaB0 beta1 phi gamma delta mu nu eta];

tspan=[0:1:tend];
options=[];
[t,y]=ode45(@ODEs,tspan,initcond,options,parameters);

initcond=y(end,:);

end
