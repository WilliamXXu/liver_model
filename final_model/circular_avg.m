function [res]=circular_avg(matrx,par)
start=par.sinusoid+1;

matrx(1:start,:)=[];
avg=mean(matrx.');

weight=start:(start+size(matrx,1)-1);
res=dot(weight,avg)/sum(weight);
end