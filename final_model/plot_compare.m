function result=plot_compare(in)
%{mito_type [group group_simple] simulation_hour frequency}
%{7 [6 60] 24 15}

%{7 [4 40] 6 10}

% {7 [1 2 9] 4 5}
% {7 [4 5] 6 10}

para=parameterGenerate();
typ=in{1};
list=in{2};
length=in{3};
pinlv=in{4};
x=0:pinlv:(60*length);

%list=[list 19];
    avg=[];
    avg1=[];
    for j=x
        res=load(append('data/',num2str(list(1)),'/mito',num2str(typ),'_',num2str(j),'.mat')).res;
        %res=load(append('experiment/data/',nu
        % m2str(i),'/resp','_',num2str(j),'.mat')).resp;
        avg=[avg circular_avg(res,para)];
        res1=load(append('data/',num2str(list(2)),'/mito',num2str(typ),'_',num2str(j),'.mat')).res;
         avg1=[avg1 mean(res1)];
    end
        avg=avg*100;
        avg1=avg1*100;
        result=avg-avg1;
    plot(x,avg);
    hold on
    plot(x,avg1);



%yline(25,'--','Fatal');
hold off
ax = gca;
ax.FontSize = 17;
title('Simulation After Drug Intake (uncoupler)')
xlabel('Minutes') 
ylabel('% of Basal cell density') 

lgd=legend('normal','simplified','Orientation','horizontal','Location','southwest'); %'tol den','tol','enta den','enta'

%lgd=legend('tolcapone','entacapone','Orientation','horizontal','Location','southwest');
lgd.FontSize = 10;
axis tight
end

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
