function [sum]=q_measure(matri)
sz=size(matri);
N=sz(1)-1;
last=matri(N+1,:);
second=matri(N,:);
sum=0;
   for j=2:(sz(2)-1)
    sum=sum+last(j)-second(j);
   end
if sz(2)==sz(1)
   sum=sum+0.5*(last(1)+last(sz(2))-second(1)-second(sz(2)));
else  
    sum=sum+(last(1)+last(sz(2))-second(1)-second(sz(2)));
    sum=(sum*N)/sz(2); 
end
sum=1-sum;

end