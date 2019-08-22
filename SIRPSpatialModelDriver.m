

runs=100;

tmax=4000;

tFinal=0:1:tmax;
SFinal=zeros(length(tFinal),runs);
IFinal=zeros(length(tFinal),runs);
RFinal=zeros(length(tFinal),runs);
PFinal=zeros(length(tFinal),runs);

tic
for i=1:runs
    i
    [Scount,Icount,Rcount,Pcount,Tcount]=SIRPSpatialModel;
    
    j=0;
    for t=tFinal
        j=j+1;
        tPos=(Tcount<=t);
        temp=Tcount.*tPos;
        [~,x]=max(temp);
        SFinal(j,i)=Scount(x);
        IFinal(j,i)=Icount(x);
        RFinal(j,i)=Rcount(x);
        PFinal(j,i)=Pcount(x);
    end
    
    
end

toc

% save(['Spatialruns',num2str(runs),'S.mat'],'SFinal')
% save(['Spatialruns',num2str(runs),'I.mat'],'IFinal')
% save(['Spatialruns',num2str(runs),'P.mat'],'PFinal')
% save(['Spatialruns',num2str(runs),'Cavg.mat'],'Cavg')
% save(['Spatialruns',num2str(runs),'Pavg.mat'],'Pavg')


figure(1)
hold on
plot(tFinal,mean(SFinal'),'b','linewidth',2)
plot(tFinal,mean(IFinal'),'r','linewidth',2)
plot(tFinal,mean(RFinal'),'g','linewidth',2)
PlotFont

figure(2)
hold on
plot(tFinal,mean(PFinal'),'m','linewidth',2)
PlotFont
