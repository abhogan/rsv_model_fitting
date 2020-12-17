
function SEIR_RSV_1yearperiod_2ages

% Simple SEIR model for RSV transmission - with births and deaths, two age
% classes. One year period.
% Author: Alexandra B Hogan (2016)

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


            
