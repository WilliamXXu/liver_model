function holder=validate_plot5()  %uncp atp
t = readtable('data/druguncp_atp.csv');
drug=2;
output=4;
res=[];

for i=1:size(t.Var1,1)
    res=[res;validation([drug t.Var1(i) output])];
end


% dec= @(b,x)(b(1)./(b(1)+x));
% inc=@(b,x)(x*b(1)/(x+b(2)));   %MM kinetics
% beta0=[0.2];
% beta = nlinfit(t.Var1,t.Var2/100,dec,beta0)
% 
% beta2=nlinfit(t.Var1,res,dec,beta0)

semilogx(t.Var1,res*100,'-o')
hold on
semilogx(t.Var1,t.Var2/100*100,'rx','MarkerSize',20)
% hold on
% semilogx(t.Var1,dec(beta,t.Var1),'-x')
% hold on
% semilogx(t.Var1,dec(beta2,t.Var1),'-x')
%legend('simulated','original');%,'fittedOriginal','fittedSimulated')
ax = gca;
ax.FontSize = 16;
title('Effects of Uncoupler on Cellular ATP')
xlabel('FCCP Concentration / uM');
ylabel('% of Basal Cellular ATP')
lgd=legend('Simulation','Experiment','Location','southwest');
lgd.FontSize = 18;%,'pp','mid','pv');%,'fittedOriginal','fittedSimulated')
axis tight
end