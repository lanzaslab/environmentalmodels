function [Scount,Icount,Rcount,Pcount,Tcount]=SIRPSpatialModel
%2-D spatial SIR model with pathogen
%Multiple individuals can occupy a given node

tmax=4000;   %days

%Initial individuals & path
n=1000;
Sn=960;
In=40;
Rn=0;

%generate spatial domain
X=5; Y=X; N=X*Y;
Pn=40;

% Ppos=ceil(N*rand(1,Pn));
% Spos=ceil(N*rand(1,Sn));
% Ipos=ceil(N*rand(1,In));
% Rpos=ceil(N*rand(1,Rn));

Ppos=ones(1,n/N);
Spos=repmat(2:40,1,40);
Ipos=ones(1,n/N);
Rpos=[];


P=zeros(Y,X);
S=zeros(Y,X);
I=zeros(Y,X);
R=zeros(Y,X);

for i=1:Y
    for j=1:X
        P(i,j)=sum(Ppos==(j+X*(i-1)));
        S(i,j)=sum(Spos==(j+X*(i-1)));
        I(i,j)=sum(Ipos==(j+X*(i-1)));
        R(i,j)=sum(Rpos==(j+X*(i-1)));
    end
end

%rates measured per day

M=0.0001;                   %move
beta=0.00001*N;             %infect
m=0.00003;                  %birth/natural death
mu=0;                       %birth/disease death
gamma=0.0008;               %recovery
kappa=0.05;                 %shedding
delta=0.1;                  %decay


t=0;

T=0:tmax;
Ttrack=1;
Tcount=0;
Scount=Sn;
Icount=In;
Rcount=Rn;
Pcount=Pn;

while t<tmax
    t
    
    q=M*n+beta*sum(sum(S.*P))+2*m*n+2*mu*In+gamma*In+kappa*In+delta*Pn;
    
    if q==0
        t=tmax;
        break
    end
    
    %exponentially distributed time step
    dt=exprnd(1/q);
    if t+dt>tmax
        t=tmax;
        break
    else
        t=t+dt;
    end
    
    %choose event type
    
    q1=M*n/q;
    q2=q1+beta*sum(sum(S.*P))/q;
    q3=q2+(m*n+mu*In)/q;
    q4=q3+m*n/q;
    q5=q4+mu*In/q;
    q6=q5+gamma*In/q;
    q7=q6+kappa*In/q;
    q8=q7+delta*Pn/q;
    
    Q=rand;
    
    if Q<q1                                     %move
        box=ceil(n*rand);                       %choose box
        [boxi,boxj]=choosebox(S+I+R,box);
        %choose S,I,R
        pop=rand*(S(boxi,boxj)+I(boxi,boxj)+R(boxi,boxj));
        if pop<S(boxi,boxj)
            S=movement(S,boxi,boxj);            %choose direction'
        elseif pop<(S(boxi,boxj)+I(boxi,boxj))
            I=movement(I,boxi,boxj);
        else
            R=movement(R,boxi,boxj);
        end
        
    elseif Q<q2                                 %infection
        total=sum(sum(S.*P));
        box=rand*total;
        [boxi,boxj]=choosebox(S.*P,box);
        
        Sn=Sn-1;
        In=In+1;
        S(boxi,boxj)=S(boxi,boxj)-1;
        I(boxi,boxj)=I(boxi,boxj)+1;
        
    elseif Q<q3                                  %birth
        boxi=ceil(rand*Y);
        boxj=ceil(rand*X);
        
        n=n+1;
        Sn=Sn+1;
        S(boxi,boxj)=S(boxi,boxj)+1;
        
    elseif Q<q4                                 %natural death
        box=ceil(n*rand);                       %choose box
        [boxi,boxj]=choosebox(S+I+R,box);
        %choose S,I,R
        pop=rand*(S(boxi,boxj)+I(boxi,boxj)+R(boxi,boxj));
        if pop<S(boxi,boxj)
            Sn=Sn-1;                            %choose death
            S(boxi,boxj)=S(boxi,boxj)-1;
        elseif pop<(S(boxi,boxj)+I(boxi,boxj))
            In=In-1;
            I(boxi,boxj)=I(boxi,boxj)-1;
        else
            Rn=Rn-1;
            R(boxi,boxj)=R(boxi,boxj)-1;
        end
        n=n-1;
        
    elseif Q<q5                                  %disease death
        box=ceil(In*rand);                       %choose box
        [boxi,boxj]=choosebox(I,box);
        
        n=n-1;
        In=In-1;
        I(boxi,boxj)=I(boxi,boxj)-1;
        
        
    elseif Q<q6                                  %recovery
        box=ceil(In*rand);                       %choose box
        [boxi,boxj]=choosebox(I,box);
        
        In=In-1;
        Rn=Rn+1;
        I(boxi,boxj)=I(boxi,boxj)-1;
        R(boxi,boxj)=R(boxi,boxj)+1;
        
    elseif Q<q7                                  %shedding
        box=ceil(In*rand);                       %choose box
        [boxi,boxj]=choosebox(I,box);
        
        Pn=Pn+1;
        P(boxi,boxj)=P(boxi,boxj)+1;
        
    elseif Q<q8                                  %decay
        box=ceil(Pn*rand);                   %choose box
        [boxi,boxj]=choosebox(P,box);
        
        Pn=Pn-1;
        P(boxi,boxj)=P(boxi,boxj)-1;
        
    end
    
    
%     if T(Ttrack+1)<t
        Ttrack=Ttrack+1;
        Tcount(Ttrack)=t;
        Scount(Ttrack)=Sn;
        Icount(Ttrack)=In;
        Rcount(Ttrack)=Rn;
        Pcount(Ttrack)=Pn;
%     end
    
end


    function [boxi,boxj]=choosebox(space,number)
        
        count=0;
        for i=1:Y
            for j=1:X
                count=count+space(i,j);
                if count>=number
                    boxi=i;
                    boxj=j;
                    break
                end
            end
            if count>=number
                break
            end
        end
        
    end

    function space=movement(space,boxi,boxj)
        
        dir=rand;
        if dir<0.25         %north
            if boxi~=1
                space(boxi,boxj)=space(boxi,boxj)-1;
                space(boxi-1,boxj)=space(boxi-1,boxj)+1;
            end
        elseif dir<0.5      %east
            if boxj~=X
                space(boxi,boxj)=space(boxi,boxj)-1;
                space(boxi,boxj+1)=space(boxi,boxj+1)+1;
            end
        elseif dir<0.75     %south
            if boxi~=Y
                space(boxi,boxj)=space(boxi,boxj)-1;
                space(boxi+1,boxj)=space(boxi+1,boxj)+1;
            end
        else                %west
            if boxj~=1
                space(boxi,boxj)=space(boxi,boxj)-1;
                space(boxi,boxj-1)=space(boxi,boxj-1)+1;
            end
        end
        
    end
% 
figure(1)
hold on
plot(Tcount,Scount,'b','linewidth',2)
plot(Tcount,Icount,'r','linewidth',2)
plot(Tcount,Rcount,'g','linewidth',2)
title(['M=',num2str(M),', N=',num2str(N)])
PlotFont

figure(2)
hold on
plot(Tcount,Pcount,'m','linewidth',2)
title(['M=',num2str(M),', N=',num2str(N)])
PlotFont

S
I
R
P

end