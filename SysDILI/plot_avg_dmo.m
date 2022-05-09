typ=5;
pinlv=5;

x=0:pinlv:450;
list=[14 11 10 17 18 15];
%list=[list 20];
for i=list
    avg=[];
    for j=x
        res=load(append('test',num2str(i),'/data/mito',num2str(typ),'_',num2str(j),'.mat')).res;
        avg=[avg mean(mean(res(5:17,:)))];
    end
    
    plot(x,avg*100);
    hold on

end
yline(50,'--','Toxic');
hold off
ax = gca;
ax.FontSize = 17;
title('Simulated Cellular ATP After Drug Intake (DMO)')
xlabel('Minutes') 
ylabel('% of Basal Cellular ATP') 
lgd=legend('A','B','C','D','E','F');%'Orientation','horizontal');
lgd.FontSize = 10;
axis tight

%legend('9','10','11','12','13','14','15','16')


% typ=5;
% pinlv=5;
% 
% x=0:pinlv:450;
% list=10:18;
% %list=[list 20];
% for i=list
%     avg=[];
%     for j=x
%         res=load(append('test',num2str(i),'/data/mito',num2str(typ),'_',num2str(j),'.mat')).res;
%         avg=[avg mean(mean(res(5:17,:)))];
%     end
%     
%     plot(x,avg*100);
%     hold on
% 
% end
% yline(50,'--','Criterion');
% hold off
% 
% title('Simulated Cellular ATP After Drug Intake (DMO)')
% xlabel('Minutes') 
% ylabel('% of Basal Cellular ATP') 
% legend('known\_0.5','unknown\_B','unknown\_D','unknown\_A','Aripiprazole','control\_1A','unknown\_C','known\_1','control\_1B');%,'unknown\_E');
% %legend('9','10','11','12','13','14','15','16')
