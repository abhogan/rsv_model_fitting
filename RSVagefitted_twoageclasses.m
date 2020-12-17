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
