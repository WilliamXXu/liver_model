function details=trialDetails()

compoundBook=compoundDetails();





details=dictionary();

s=Trial();
s.id=6;
s.c_max=5;
s.compound=compoundBook('aripiprazole');
s.options=[1 0];
s.simulation_hour=24;
s.monitoring_freq_min=15;

details(s.id)=s;


s=Trial();
s.id=60;
s.c_max=5;
s.compound=compoundBook('aripiprazole');
s.options=[1 1];
s.simulation_hour=24;
s.monitoring_freq_min=15;
details(s.id)=s;



s=Trial();
s.id=9;
s.c_max=1.3;
s.compound=compoundBook('nefazodone_old');
s.options=[1 0];
s.simulation_hour=5;
s.monitoring_freq_min=5;

details(s.id)=s;


s=Trial();
s.id=900;
s.c_max=0.417;
s.compound=compoundBook('nefazodone_old');
s.options=[1 0];
s.simulation_hour=24*3;
s.monitoring_freq_min=20;
s.dosing_freq_hour=12;
details(s.id)=s;

s=Trial();
s.id=9000;
s.c_max=0.417;
s.compound=compoundBook('nefazodone');
s.options=[1 0];
s.simulation_hour=24*3;
s.monitoring_freq_min=20;
s.dosing_freq_hour=12;
details(s.id)=s;

s=Trial();
s.id=90000;
s.c_max=0.417;
s.compound=compoundBook('nefazodone_old');
s.options=[1 1];
s.simulation_hour=24*7;
s.monitoring_freq_min=20;
s.dosing_freq_hour=12;
details(s.id)=s;




s=Trial();
s.id=4;
s.c_max=20.8;
s.compound=compoundBook('tolcapone');
s.options=[1 0];
s.simulation_hour=6;
s.monitoring_freq_min=10;
details(s.id)=s;


s=Trial();
s.id=40;
s.compound=compoundBook('tolcapone');
s.options=[1 1];
s.simulation_hour=6;
s.monitoring_freq_min=10;
details(s.id)=s;


s=Trial();
s.id=5;
s.c_max=4;
s.compound=compoundBook('entacapone');
s.options=[1 0];
s.simulation_hour=6;
s.monitoring_freq_min=10;
details(s.id)=s;


end
