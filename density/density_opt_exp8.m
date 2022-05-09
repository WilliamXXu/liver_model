function density_opt_exp8()   %sqp algo, two data, two stage model no delay


atp_decrea5={};
den5={};
%exclude k=7 temporarily 
temp=readtable(append('data2/atp5_',num2str(5),'.csv'));
for k=1:8
    atp=readtable(append('data2/atp5_',num2str(k),'.csv'));
    basal=atp.Var2(1);
    if k==8
        atp.Var2=atp.Var2/temp.Var2(1);
    else
        atp.Var2=atp.Var2/basal;
    end
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

%plot(1:600,atp_full,'-x')


%  M=eye(4);
%  M(4,3)=1;
% options = odeset('Mass',M);
 

initial=[1 0 0 0];
% r=[ 0.000387803817035  53.723010119627133   6.275002494446798   5.992715304984040   3.790019123820707];
%r=[80.455047279398968   1.073269211540708  88.599050737349799   0.004998720928621   0.167137548667601   0.876728114047647];
r=[12 1 1 1 1 1];
ub=[ 100 5 100  20 10 30 ];
b=zeros(numel(r),1);
%options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective');


rng default % For reproducibility
opts = optimoptions(@fmincon,'Display','iter','Algorithm','sqp','MaxFunctionEvaluations',1500);%,'StepTolerance',1e-6,'OptimalityTolerance',1e-5); %  ,);
problem = createOptimProblem('fmincon','objective',@objectiv,'x0',r,'lb',b.','ub',ub,'options',opts);

ms = MultiStart('StartPointsToRun','bounds-ineqs','UseParallel',1,'Display','iter','MaxTime',1000);
[r,f] = run(ms,problem,30);

%opt = optimoptions('patternsearch','MaxFunctionEvaluations',1200,'FunctionTolerance',1e-4,'Display','iter');
%r=patternsearch(@objectiv,r,[],[],[],[],b.',ub,[],opt);

format long
objectiv(r)       
 
r 



  function [res]=objectiv(para)
res=0;
para=[para(1:3) 0 para(4:6)];
for k=1:8
    
    if k~=7
       flag=1;
       len=numel(den5{k}.Var1);
        simu=[];
        last=initial;
for time=0:120
    
    temp=atp_decrea5{k};
    [t,y]=ode15s(@(t,y)den_exp4(t,y,para,temp(time+1)),[0 1],last);
    last=y(size(y,1),:);

        if abs(time-den5{k}.Var1(flag))<0.5
            flag=(flag+1)*(flag<len)+(flag==len);
            simu=[simu;sum(last(1:3))];
        
        end
end

res=res+norm(den5{k}.Var2-simu)^2;


    end


end

for k=1:3
    if k~=7
       flag=1;
       len=numel(den6{k}.Var1);
        simu=[];
        last=initial;
for time=0:180
     temp=atp_decrea6{k};
    [t,y]=ode15s(@(t,y)den_exp4(t,y,para,temp(time+1)),[0 1],last);
    last=y(size(y,1),:);

        if abs(time-den6{k}.Var1(flag))<0.5
            flag=(flag+1)*(flag<len)+(flag==len);
            simu=[simu;sum(last(1:3))];
        
        end
end

res=res+norm(den6{k}.Var2-simu)^2;


    end

    
end



    end


%DELAyed response?    ode for delay factor
end