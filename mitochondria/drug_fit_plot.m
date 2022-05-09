range=[0.001 0.01 0.015 0.03 0.06 0.12 0.25 0.3 0.5 0.6 0.75 0.8 1 1.2 1.5 2 2.5 3 3.5 4 4.5 5 5.5];
x=arrayfun(@drug_fit,range)
save('nominal.mat',x)
