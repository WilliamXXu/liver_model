names=['drugatpase_mmp.csv' 'drugetc_atp.csv' 'drugetc_mmp' 'druguncp_atp'];


t = readtable('drugatpase_mmp.csv');
ling1=size(t.Var1,1);
drug1=[zeros(2,ling1);t.Var1.'];
target1=[zeros(2,ling1);t.Var2.';zeros(1,ling1)];

t = readtable('drugetc_mmp.csv');
ling2=size(t.Var1,1);
drug2=[t.Var1.';zeros(2,ling2)];
target2=[zeros(2,ling2);t.Var2.';zeros(1,ling2)];

t = readtable('drugetc_atp.csv');
ling3=size(t.Var1,1);
drug3=[t.Var1.';zeros(2,ling3)];
target3=[zeros(3,ling3);t.Var2.'/100];

t = readtable('druguncp_atp.csv');
ling4=size(t.Var1,1);
drug4=[zeros(1,ling4);t.Var1.';zeros(1,ling4)];
target4=[zeros(3,ling4);t.Var2.'/100];

drug=[drug1 drug2 drug3 drug4];
target=[target1 target2 target3 target4];
save('training_graphic.mat','target','drug');