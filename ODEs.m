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