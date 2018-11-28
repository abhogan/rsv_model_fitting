
% Simple SEIR model for RSV transmission - with births and deaths, two age
% classes. One year period.
% Author: Alexandra B Hogan (2016)

function SEIR_RSV_1yearperiod_2ages
format long
tend=100;   % end time of calculations in years
tspan=0:0.001:tend;

%disease parameters
delta=91.4787;
gamma=40.1099;
nu=1.5851;
mu=0.0135;
beta1=0.522;
alpha=0.228;
eta=1;
beta0=3215;

% initial values
I10=0.00001;
S10=0.1;
E10=0.00011;
R10=0;
S20=S10;
E20=E10;
I20=I10;
R20=R10;

% solve ODES using Matlab's ode45 integrator
param=[gamma delta nu mu eta beta0 beta1 alpha];
[t,y1]=ode45(@sirODESIJ_RSV2,tspan,[S10 E10 I10 R10 S20 E20 I20 R20],[],param);

% plot Infectives
figure(1)
plotstart=tend-10;
plot(t-plotstart,y1(:,3)./(y1(:,1)+y1(:,2)+y1(:,3)+y1(:,4)),'b-','LineWidth',1.5)
hold on
plot(t-plotstart,y1(:,7)./(y1(:,5)+y1(:,6)+y1(:,7)+y1(:,8)),'m-','LineWidth',1.5)
xlabel('Years','FontSize',12)
ylabel('Proportion Infected','FontSize',12)
axis([0,10,0,0.1])

figure(2)
plot(t,y1(:,3))

y1(end,:)

end

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

            