function [sum]=q_measure_test(matri)
sz=size(matri);
N=sz(1)-1;
last=matri(2,:);
second=matri(1,:);

sum=0.5*(last(1)+last(sz(2))-second(1)-second(sz(2)));
for k=2:(sz(2)-1)
    sum=sum+last(k)-second(k);
end

if sz(2)~=sz(1)
   sum=(sum*N)/sz(2); 
end

end