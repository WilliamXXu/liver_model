simTime=600;
ind=3;
res1=drug_test(1,1.425,simTime)*100;
res2=drug_test(1,2.849,simTime)*100;
res3=drug_test(1,6.049,simTime)*100;
res4=drug_test(1,23.267,simTime)*100;
res5=drug_test(1,38.462,simTime)*100;
res6=drug_test(1,147.929,simTime)*100;


plot(1:1:simTime,res1(:,ind));
hold on
plot(1:1:simTime,res2(:,ind));
hold on
plot(1:1:simTime,res3(:,ind));
hold on
plot(1:1:simTime,res4(:,ind));
hold on
plot(1:1:simTime,res5(:,ind));
hold on
plot(1:1:simTime,res6(:,ind));
hold on
%plot(1:1:600,res6(:,ind));
yline(50,'--','Toxic');
hold off
ax = gca;
ax.FontSize = 17;
title('Simulated Cellular ATP After Drug Intake (DM)')
xlabel('Minutes') 
ylabel('% of Basal Cellular ATP') 
legend('A','B','C','D','E','F','Orientation','horizontal');




% simTime=600;
% ind=3;
% res1=drug_test(1,1.425,simTime)*100;
% res2=drug_test(1,6.049,simTime)*100;
% res3=drug_test(1,23.267,simTime)*100;
% res4=drug_test(1,38.462,simTime)*100;
% res5=drug_test(1,147.929,simTime)*100;
% 
% 
% 
% plot(1:1:simTime,res1(:,ind));
% hold on
% plot(1:1:simTime,res2(:,ind));
% hold on
% plot(1:1:simTime,res3(:,ind));
% hold on
% plot(1:1:simTime,res4(:,ind));
% hold on
% plot(1:1:simTime,res5(:,ind));
% hold on
% %plot(1:1:600,res6(:,ind));
% yline(50,'--','Criterion');
% hold off
% title('Simulated Cellular ATP After Drug Intake (DM)')
% xlabel('Minutes') 
% ylabel('% of Basal Cellular ATP') 
% legend('Aripiprazole','known\_0,5','known\_1','control\_1A','control\_1B')