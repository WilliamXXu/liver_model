
function sysDILI(in)   % in={c_max para_val drug_type tmax_hour halfLife_hour 1/0density testId simulation_hour monitoring_freq_min}

%nefazodone {1.5 [0.2456] 1 1.5 2.5 1 9 5 5}
%tolcapone {20.8 [397.0337 5.0000e+04] 2 2 2.9 [1 0] 4 6 10}
% entacapone {4 [56.2668 5.0000e+04] 2 1 0.8 [1 0] 5 6 10}
% ---------------------- input processing -----------------------------------------

details=trialDetails();
s=details(in);


par=parameterGenerate();
method=Method(par);

listHepato=[];
    for i=2:(par.length+1)      
        for j=(par.sinusoid+3):(par.radius+1)
            n=method.up(i,j);
            listHepato=[listHepato n];  
        end
    end



option=s.options;

drug=Drug(s);
oxygen=Oxygen(par,option(1),option(2),listHepato);
mitochondria=Mitochondria(par,option(1),option(2),listHepato);
saver=Saver(par,option(2));

saver.id=in;
saver.frequency=s.monitoring_freq_min;


%----------------------  zoning -------------------------------




%---------------------- iterate ---------------------------------------

saver.makeFolder()
count=0;
tim=0;
%while m_t>error(2)
while tim<=(s.simulation_hour*60)
    tim=count*par.time(2)/60

    tic    
    drug.iterate(tim);
    if option(2)==0
        oxygen.iterate(mitochondria.M,drug);
    end
    mitochondria.iterate(oxygen.C,drug);
    toc

   
if mod(count,saver.frequency)==0
    saver.snapshot(count,oxygen,mitochondria);
end
    count=count+1;
    
end


end