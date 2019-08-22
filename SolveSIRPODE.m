function SolveSIRPODE

tmax=500;

%Initial individuals & path
P0=40;
S0=960;
I0=40;
R0=0;

inits=[S0,I0,R0,P0];

M=0.25;                    %move
beta=0.00001;             %infect
m=0.00003;                %birth/natural death
mu=0;             %birth/disease death
gamma=0.007;          %recovery
kappa=1;          %shedding
delta=0.01;          %decay

options = odeset('RelTol',1e-12,'Abstol',1e-12);
[t,sol]=ode45(@SIRP,[0 tmax],inits,options);

figure(1)
hold on
plot(t,sol(:,1),'b','linewidth',2)
plot(t,sol(:,2),'r','linewidth',2)
plot(t,sol(:,3),'g','linewidth',2)
xlabel('Time')
ylabel('Populations')
PlotFont

figure(2)
hold on
plot(t,sol(:,4),'m','linewidth',2)
xlabel('Time')
ylabel('Pathogen')
PlotFont

    function dpop=SIRP(t,pop)
        
        S=pop(1);
        I=pop(2);
        R=pop(3);
        P=pop(4);
        
        dS=m*(S+I+R)+mu*I-beta*S*P-m*S;
        dI=beta*S*P-(gamma+m+mu)*I;
        dR=gamma*I-m*R;
        dP=kappa*I-delta*P;
        
        dpop=[dS,dI,dR,dP]';
        
        
    end

end
