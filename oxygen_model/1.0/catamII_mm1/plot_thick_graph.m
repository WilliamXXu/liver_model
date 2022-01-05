load('q5_1.mat');
y=res2;
sz=size(y);
step=3;
for in=1:sz(3)
   x=y(:,:,in);
   sec=x(201-(in-1)*step-20:201-(in-1)*step,:);
   figure();
   sec(sec==0)=0.5;
   heatmap(sec.','Colormap',copper);
  % xlabel(num2str(in));
end

