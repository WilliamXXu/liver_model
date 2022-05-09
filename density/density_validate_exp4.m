function density_validate_exp4()
max_t=180;

atp_decrea5={};
den5={};
for k=1:8
    atp=readtable(append('data2/atp5_',num2str(k),'.csv'));
    basal=atp.Var2(1);
    atp.Var2=atp.Var2/basal;
    atp_decrea5{k}=1-interp1(atp.Var1,atp.Var2,0:120,"spline");
end

for k=1:8
    den=readtable(append('data2/den5_',num2str(k),'.csv'));
    basal=den.Var2(1);
    den.Var2=den.Var2/basal;
    den5{k}=den;
end



atp_decrea6={};
den6={};
for k=1:3
    atp=readtable(append('data2/atp6_',num2str(k),'.csv'));
    basal=atp.Var2(1);
    atp.Var2=atp.Var2/basal;
    atp_decrea6{k}=1-interp1(atp.Var1,atp.Var2,0:180,"spline");
end

for k=1:3
    den=readtable(append('data2/den6_',num2str(k),'.csv'));
    basal=den.Var2(1);
    den.Var2=den.Var2/basal;
    den6{k}=den;
end

%atp_decrea=repmat(0.9,max_t+1,1);
ind=2;
atp_decrea=atp_decrea6{ind};
d=den6{ind};

 
initial=[1 0 0 0];
%r=[ 19.806881436793049   1.285234137013759 21.080026927983923  10.499999951160587];
r=[80.455047279398968   1.073269211540708  88.599050737349799 0  0.004998720928621   0.167137548667601   0.876728114047647];
temp=objectiv(r);


plot(0:max_t,sum(temp(:,1:3).').','-');
%  plot(0:max_t,temp(:,4));
  hold on
  plot(d.Var1,d.Var2,'rx','MarkerSize',20);
  hold off
ax = gca;
ax.FontSize = 16;
  ylabel('% of Basal Hepatocyte Density');
xlabel('Time (minute)');
title('Cell Death by ATP Depletion, Experiment 1');
lgd=legend('Simulation','Experiment','Location','southwest');
lgd.FontSize = 18;
axis tight

    function [simu]=objectiv(para)
        simu=[];
last=initial;
for time=0:max_t
    [t,y]=ode15s(@(t,y)den_exp4(t,y,para,atp_delay(atp_decrea,time,para(4))),[0 1],last);
    last=y(size(y,1),:);
    simu=[simu;last];
    
           % simu=[simu;sum(last(1:3))];
  
end
    end

    function [res]=atp_delay(atp_decrea,time,delay)
        time=time-round(delay);
        time=(time>=0)*time;
        res=atp_decrea(time+1);
    end
% para=[0.1 10 0.01 0.001 0.001 0.1 10 1];
% atp_dec=0.1;
%plot(t,sum(y,2))

%DELAyed response?    ode for delay factor
end