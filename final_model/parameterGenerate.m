function par=parameterGenerate()

par.time=[4.5e-3 60];      %timesteps
par.error=[1e-4 1e-4 1e-5];      %[oxygen cellular ATP drug] mmHg/second  unit/second
par.batch_number=[4 100];

par.unit=2; % miu m, space step
hepatocyte_diamter=27.5;   % miu m,
sinusoid_diameter=12.5;
sinusoid_length=275;

par.sinusoid=round(sinusoid_diameter/2/par.unit);
par.radius=round(par.sinusoid+hepatocyte_diamter/par.unit);
par.length=round(sinusoid_length/par.unit);
par.total=(par.length+2)*(par.radius+2)  ;





par.velo=400;   %miu/sec?

d1=1620;     % diffusion coeff of sinusoid    miu m^2/sec
d2=1600 ;    %diffusionx coeff of hepatocyte     3600(not sure)

par.miu1=par.time(1)*d1/(par.unit^2);
par.miu2=par.time(1)*d2/(par.unit^2);

par.basal=[0.6 199.6*5 337.2 0.3 0.1722 4.0 1 ];
par.ocr=199.6/60/9*100/1.34;

par.oxy_pp=70;     %mmHg 
par.hemo_pp=4*2300/1.34/(1+(26/par.oxy_pp)^2.73); %mmHg

conversion=1.34e-3;    %mmHg to mol/m3, Henry's law
par.v_respire=0.044/conversion;  %mmHg/sec  
par.k_respire=6.24e-3/conversion; %mmHg

par.coeff=[par.ocr/0.6 198.650/0.6 198.650/par.basal(2) 1/90 (par.basal(2)+par.basal(3))/90/4 ];   %%[etc2ocr etc2mmp atp2mmp atp2atp_cell  atp_cell_consump
par.para_den=[80.455047279398968   1.073269211540708  88.599050737349799   0.004998720928621   0.167137548667601   0.876728114047647];
par.para_mito=load('r_9.mat').r;

end