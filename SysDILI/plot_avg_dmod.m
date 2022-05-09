typ=6;
pinlv=5;

x=0:pinlv:450;
list=[5 4 2 1 9 8];
%list=[list 19];
for i=list
    avg=[];
    for j=x
        res=load(append('test',num2str(i),'/data/mito',num2str(typ),'_',num2str(j),'.mat')).res;
        avg=[avg mean(mean(res(5:17,:)))];
    end
    
    plot(x,avg*100);
    hold on

end
yline(25,'--','Fatal');
hold off
ax = gca;
ax.FontSize = 17;
title('Simulated Cell Density After Drug Intake (DMOD)')
xlabel('Minutes') 
ylabel('% of Basal Cell Density') 
lgd=legend('A','B','C','D','E','F');%'Orientation','horizontal');
lgd.FontSize = 10;
axis tight

% typ=5;
% pinlv=5;
% 
% x=0:pinlv:450;
% list=1:9;
% %list=[list 19];
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
% yline(25,'--','Criterion');
% hold off
% 
% title('Simulated Hepatocyte ATP After Drug Intake (DMOD)')
% xlabel('Minutes') 
% ylabel('% of Basal Cell Density') 
% legend('known\_1','known\_0.5','unknown\_C','unknown\_B','Aripiprazole','unknown\_A','unknown\_D','control\_1A','control\_1B');%,'unknown\_E');
% 
