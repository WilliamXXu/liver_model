%{
a=eye(5);
b=zeros(5);
save('matrices','a','b')

sa=[];
tic
for i=1:999
    sa=[sa;i i 1];
end
toc
sa=[sa;1000 1000 0.01];
a=sparse(sa(:,1),sa(:,2),sa(:,3),1000,1000);
%}

load('result_modified.mat');
res0=res;
load('result_modified_initial.mat');
sum(sum(res0-res))
