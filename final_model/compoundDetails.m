function details=compoundDetails()

details=dictionary();


c=Compound();
c.name='aripiprazole';

c.para_val=[0.0765];
c.drug_type=1;
c.tmax_hour=4;
c.halfLife_hour=75;
details(c.name)=c;

c=Compound();
c.name='nefazodone_old';

c.para_val=[0.2456];
c.drug_type=1;
c.tmax_hour=1;
c.halfLife_hour=1.4;
details(c.name)=c;

c=Compound();
c.name='nefazodone';

c.para_val=[1.5654];
c.drug_type=1;
c.tmax_hour=1;
c.halfLife_hour=1.4;
details(c.name)=c;





c=Compound();
c.name='tolcapone';

c.para_val=[397.0337 5.0000e+04];
c.drug_type=2;
c.tmax_hour=2;
c.halfLife_hour=2.9;
details(c.name)=c;


c=Compound();
c.name='entacapone';

c.para_val=[56.2668 5.0000e+04];
c.drug_type=2;
c.tmax_hour=1;
c.halfLife_hour=0.8;
details(c.name)=c;


end