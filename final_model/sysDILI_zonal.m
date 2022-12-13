
function sysDILI_zonal(in)   % in={c_max para_val drug_type tmax_hour halfLife_hour 1/0density testId simulation_hour monitoring_freq_min}
%aripiprazole {5 [0.0765] 1 4 75 [1 1] 60 24 15}
%tolcapone {20.8 [397.0337 5.0000e+04] 2 2 2.9 [1 0] 40 6 10}


%nefazodone {1.5 [0.2456] 1 1.5 2.5 1 9 5 5}
%tolcapone {20.8 [397.0337 5.0000e+04] 2 2 2.9 1 4 6 10}
% entacapone {4 [56.2668 5.0000e+04] 2 1 0.8 1 5 6 10}
% ----------------------  canvas -----------------------------------------


par=parameterGenerate();
option=in{6};


drug=Drug();
oxygen=Oxygen(par,option(1),option(2));
mitochondria=Mitochondria(par,option(1),option(2));
saver=Saver(par,option(2));

saver.id=in{7};
saver.frequency=in{9};
drug.type=in{3};
drug.coeff{in{3}}=in{2};
drug.c_max=in{1};
drug.t_max=in{4};
drug.half_life=in{5};
simulationTime=in{8};




%---------------------- iterate ---------------------------------------

saver.makeFolder()
count=0;
tim=0;
%while m_t>error(2)
while tim<=(simulationTime*60)
    tim=count*par.time(2)/60


    drug.iterate(tim);
    mitochondria.iterate(oxygen.C,drug);
    
    if mod(count,saver.frequency)==0
         saver.snapshot(count,oxygen,mitochondria);
    end
    count=count+1;
    
end


end