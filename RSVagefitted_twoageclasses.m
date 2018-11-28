% SEIRS model for RSV transmission - with births and deaths, for two age
% classes. This code fits the transmission model to RSV
% laboratory-confirmed cases in Western Australia, for two age classes:
% 0-11 months and 12-23 months. The RSV data is omitted and synthetic data,
% generated from the model with added noise, is instead used (which is why
% the model fits the data so quickly why there is not a unique fitted
% solution - the model will fit for any number of starting parameter
% combinations).

% Code authors: Alexandra B Hogan, Roslyn I Hickson, Geoffry N Mercer,
% Hannah C Moore, Kathryn Glass
% Year: 2016
% Related publication: Hogan AB, Glass K, Moore HC, Anderssen R. Exploring
% the dynamics of respiratory syncytial virus (RSV) transmission in
% children. Theor. Popul. Biol. 2016;110:78–85.
% Available from: http://dx.doi.org/10.1016/j.tpb.2016.04.003

% function descriptions

% ODEs: 2 age class RSV ordinary differential equation model
% calcsoln: solves the ODEs. Note that rates are in weeks.
% calcinitsol: solves the ODEs with a longer burnin time in the
% first instance, and saves the final solution as the initial conditions
% for subsequent ODE solves - makes the fitting routine more efficient
% errorcalc: calculates difference (error) between model solution and data


%% Main routine starts here
format compact

data = csvread("synthetic_data_2ages_raw.csv", 1, 0);
% data age group 1
data1 =data(:,1);
% data age group 2
data2=data(:,2);

lengdata=length(data1);
sumdata1=sum(data1);
sumdata2=sum(data2);
sumtotaldata=sumdata1+sumdata2;

% parameter guesses. Need to start with something in the right ballpark for
% the fitting to work
betaA0=52
betaB0=0.4*betaA0
beta1=0.5630
phi= -0.37
nu= 0.034
paramguess=[betaA0 betaB0 beta1 phi nu];

global initcond
calcinitsol(paramguess);

% either find the optimal solution (findoptimalsoln=1) or just run the
% model once for the paramter guesses (findoptimalsoln=0)
% Note that you may need to experiment with the tolerances - due to the
% fitting algorithm either settling on a solution too quickly, or not at
% all. Currently, the fitting algorithm isn't searching widely enough
% around betaA0 so this needs to be tweaked by manually editing the fminsearch function.
findoptimalsoln=1;
if findoptimalsoln
  options = optimset('TolX',1e-6,'TolFun',1e-6);
  [paramfitted,fval,exitflag,output]=fminsearch(@errorcalc,paramguess,options)
  [t,y]=calcsoln(paramfitted); %finds solution with those final chosen parameters
else
  [t,y]=calcsoln(paramguess);

end
% pull out model data into more sensible names
I1=y(:,3); % proportion infectious age group 1
I2=y(:,7);
J1=y(:,9); % incidence age group 1
J2=y(:,10);   

% plot the whole solution to check it has converged
figure(1)
plot(t,I1,'r-','LineWidth',1.5)
hold on
plot(t,I2,'r--','LineWidth',1.5)
hold off
print -djpeg99 RSV1test

% scale solution to have same total cases and plot solution and data
newcases(1)=0;
under12newcases(1)=0;
under24newcases(1)=0;
for i=2:length(J2)-1
    newcases(i)=J1(i)-J1(i-1)+J2(i)-J2(i-1);
    under12newcases(i)=J1(i)-J1(i-1);
    under24newcases(i)=J2(i)-J2(i-1);
end
year=52;
tend=year*72;
% define the last part of the model output to compare to the data
endmodel=[tend-lengdata+1:tend];
%endmodel=[tend-lengdata+1-52:tend-52];
summodel=sum(newcases(endmodel));
summodel1=sum(under12newcases(endmodel));
summodel2=sum(under24newcases(endmodel));
under12newcases=under12newcases*sumdata1/summodel1;
under24newcases=under24newcases*sumdata1/summodel1;

% plot solutions for specific timeframe
figure(2)
plot(t(endmodel),I1(endmodel),'r-','LineWidth',2)
hold on
plot(t(endmodel),I2(endmodel),'r--','LineWidth',2)
xlabel('Week','FontSize',12)
ylabel('Proportion of total population','FontSize',12)
legend('Group1','Group2')
hold off
print -djpeg99 RSV2test

%plot scaled solution against case data
figure(3)
plot([1:lengdata],data1,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1],...
    'MarkerSize',4,...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1],...
    'DisplayName','Data: Group 1')
hold on
plot([1:lengdata],under12newcases(endmodel),'LineWidth',3,'Color',[0.847058832645416 0.160784319043159 0],...
    'DisplayName','Model: Group 1')
hold on
plot([1:lengdata],data2,'MarkerFaceColor',[1 0 1],...
    'MarkerEdgeColor',[0.749019622802734 0 0.749019622802734],...
    'MarkerSize',4,...
    'Marker','^',...
    'LineStyle','none',...
    'Color',[0.682352960109711 0.466666668653488 0],...
    'DisplayName','Data: Group 2')
hold on
plot([1:lengdata],under24newcases(endmodel),'LineWidth',3,'LineStyle','--',...
    'Color',[0 0.600000023841858 0.200000002980232],...
    'DisplayName','Model: Group 2')
hold on
% Create xlabel
xlabel('Week','FontSize',16);

% Create ylabel
ylabel('Weekly incidence','FontSize',16);

% Create title
title({'Model versus test data'},'FontSize',18);

% Create legend
legend1 = legend('show');
set(legend1,...
    'Position',[0.666426274948517 0.735596774503307 0.191068814055637 0.189381499726327]);
axis([0 lengdata 0 190])
box off
hold off
print -djpeg99 RSV3test

%% Functions

function f=ODEs(t,y,parameters)

betaA0=parameters(1);
betaB0=parameters(2);
beta1=parameters(3);
phi=parameters(4);
gamma=parameters(5); 
delta=parameters(6); 
mu=parameters(7);
nu=parameters(8);
eta=parameters(9);

S1=y(1);  E1=y(2);   I1=y(3);  R1=y(4);
S2=y(5);  E2=y(6);   I2=y(7);  R2=y(8);

betaA=betaA0*(1+beta1*sin(2*pi*t/52+phi));
betaB=betaB0*(1+beta1*sin(2*pi*t/52+phi));
bAsumI=betaA*(I1+I2);
bBsumI=betaB*(I1+I2);

f(1) = mu - S1*(bAsumI)+ nu*R1 - eta*S1 ;
f(2) = S1*(bAsumI) - delta*E1 - eta*E1;
f(3) = delta*E1 - gamma*I1 - eta*I1;
f(4) = gamma*I1 - nu*R1 - eta*R1;
f(5) = eta*S1  + nu*R2 - S2*bBsumI - eta*S2;
f(6) = eta*E1 + S2*bBsumI - delta*E2 - eta*E2;
f(7) = eta*I1 + delta*E2 - gamma*I2 - eta*I2 ;
f(8) = eta*R1 + gamma*I2 - nu*R2 - eta*R2 ;
f(9) = S1*bAsumI;
f(10) = S2*bBsumI;
    
f=f';

end

%-------------------------------------------------------------------

function f=errorcalc(param)
global initcond
% calculate the error between model and data
data1=[];
data2=[];

lengdata=length(data1);
sumdata1=sum(data1);
sumdata2=sum(data2);
sumtotaldata=sumdata1+sumdata2;

year=52;
tend=year*72;
[t,y]=calcsoln(param);  
J1=y(:,9);  J2=y(:,10);   

newcases(1)=0;
under12newcases(1)=0;
under24newcases(1)=0;

for i=2:length(J2)-1
    newcases(i)=J1(i)-J1(i-1)+J2(i)-J2(i-1);
    under12newcases(i)=J1(i)-J1(i-1);
    under24newcases(i)=J2(i)-J2(i-1);
end
% scale so total of new cases is the same
%endmodel=[tend-lengdata+1-52:tend-52];
endmodel=[tend-lengdata+1:tend];
summodel=sum(newcases(endmodel));
summodel1=sum(under12newcases(endmodel));
summodel2=sum(under24newcases(endmodel));
under12newcases=under12newcases*sumdata1/summodel1;
under24newcases=under24newcases*sumdata1/summodel1;

% calculate the error between model and data - sum of least squares
% approach (could instead use Maximum Likelihood Estimation)
%f=sqrt(sum((data1-under12newcases(endmodel)).^2)+sum((data2-under24newcases(endmodel)).^2));

% Maximum Likelihood Estimation option (i.e. minimise negative log
% likelihood)
data=[data1' data2']';
output=[under12newcases(endmodel) under24newcases(endmodel)];
f=-sum(data.*log(output)-output);

disp([f,param])
end

%-----------------------------------------------------------------------

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
