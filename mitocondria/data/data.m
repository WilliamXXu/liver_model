%data group 1: etc drug: rotenone
drug_etc=[0.0005 0.015 0.03 0.06 0.12 0.25 0.5 1];
ocr_etc=[0.9969 0.9514 0.8739 0.8151 0.7082 0.4998 0.318 0.2];
ecar_etc=[1 1 1.1 1.2 1.3 1.5 1.5 1.6];

ling1=size(drug_etc,2);
drug1=[drug_etc; zeros(2,ling1)];
target1=[ocr_etc;ecar_etc;zeros(2,ling1)];

%data group 2: uncp drug: FCCP (normal)
drug_uncp=[0.0005 0.05 0.141 0.25 0.5 1];
ocr_uncp=[1 1 1.4 2.2 3.9 4];
ecar_uncp=[1 1.2 1.39 2 2.6 2.5];
ling2=size(drug_uncp,2);
drug2=[zeros(1,ling2);drug_uncp;zeros(1,ling2)];
target2=[ocr_uncp;ecar_uncp;zeros(2,ling2)];


%data group3: uncp drug: FCCP (potential)
drug_uncp2=[0.005 0.012 0.037 0.111 0.333 1 3 9];
potential=[0.999 0.931 0.921 1.05 1.14 0.949 0.621 0.276];
ling3=size(drug_uncp2,2);
drug3=[zeros(1,ling3);drug_uncp2;zeros(1,ling3)];
target3=[zeros(2,ling3);potential;zeros(1,ling3)];


%data group4: atpase drug: oligomycin
drug_atpase=[0.005 0.125 0.25 0.5 1 2];
ocr_atpase=[0.97 0.987 0.982 0.963 0.517 0.264];
ecar_atpase=[1 1.1 1.1 1.1 1.4 1.7];
ling4=size(drug_atpase,2);
drug4=[zeros(2,ling4);drug_atpase];
target4=[ocr_atpase;ecar_atpase;zeros(2,ling4)];

drug=[drug1 drug2 drug3 drug4];
target=[target1 target2 target3 target4];
save('training.mat','target','drug');