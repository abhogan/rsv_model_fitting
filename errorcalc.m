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