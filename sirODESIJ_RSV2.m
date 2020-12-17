function dy = sirODESIJ_RSV2(t,y,param)
% ode function for the SEIR model - with births and deaths

S1 = y(1); %susceptible
E1 = y(2); %exposed
I1 = y(3); %infected
R1 = y(4); %recovered
S2 = y(5);
E2 = y(6);
I2 = y(7);
R2 = y(8);

gamma = param(1); %1/gamma is av recovery period
delta = param(2); %1/delta is average latency period
nu = param(3); %1/nu is average duration of immunity
mu = param(4); %birth rate
eta = param(5); %ageing rate
beta0 = param(6);
beta1 = param(7);
alpha=param(8);

% beta varies seasonally
beta=beta0*(1+beta1*cos(2*pi*t));
dy=[mu - beta*S1*(I1 + I2) - eta*S1 + nu*R1
    beta*S1*(I1 + I2) - E1*(delta + eta)
    delta*E1 - I1*(eta + gamma)
    gamma*I1 - R1*(eta + nu)
    eta*S1 - alpha*S2*beta*(I1 + I2) - eta*S2 + nu*R2
    eta*E1 + alpha*S2*beta*(I1 + I2) - E2*(delta + eta)
    eta*I1 + delta*E2 - I2*(eta + gamma)
    eta*R1 + gamma*I2 - R2*(eta + nu)];
end