function holder=validate_plot10()   

Var1=[0.005 0.125 0.25 0.5 1 2];
Var2=[0.97 0.987 0.982 0.963 0.517 0.264];
t = table(Var1.',Var2.');


drug=3;
output=1;
res=[];
% pp=[];
% mid=[];
% pv=[];

for i=1:size(t.Var1,1)
%    out=validation_oxy_mito([drug t.Var1(i)]);
%     pp=[pp out(1,output)];
%     mid=[mid out(2,output)];
%     pv=[pv out(3,output)];
    res=[res;validation([drug t.Var1(i) output])];
end


% dec= @(b,x)(b(1)./(b(1)+x));
% inc=@(b,x)(x*b(1)/(x+b(2)));   %MM kinetics
% beta0=[0.02];
% beta = nlinfit(t.Var1,t.Var2,dec,beta0)
% 
% beta2=nlinfit(t.Var1,res,dec,beta0)
semilogx(t.Var1,res,'-o')
hold on
semilogx(t.Var1,t.Var2,'rx','MarkerSize',20)
hold on
% semilogx(t.Var1,pp,'-x')
% hold on
% semilogx(t.Var1,mid,'-x')
% hold on
% semilogx(t.Var1,pv,'-x')
% 
% % hold on
% % semilogx(t.Var1,dec(beta,t.Var1),'-x')
% % hold on
% % semilogx(t.Var1,dec(beta2,t.Var1),'-x')
% legend('simulated','original','pp','mid','pv');%,'fittedOriginal','fittedSimulated')
ax = gca;
ax.FontSize = 16;
title('Effects of ATPase inhibitor on Oxygen Consumption Rate (OCR)')
xlabel('Oligomycin Concentration / uM');
ylabel('% of Basal OCR')
lgd=legend('Simulation','Experiment','Location','southwest');%,'pp','mid','pv');%,'fittedOriginal','fittedSimulated')
axis tight
lgd.FontSize = 18;
end