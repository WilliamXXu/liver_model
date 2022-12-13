classdef Drug < handle
    properties 
        dose

    end
    properties 
        compound
        coeff
        c_max
        dosing_freq_hour
    end
    methods
        function obj=Drug(trial)
            compound=trial.compound;
            obj.compound=compound;
            obj.coeff={1 [1 1] 1};
            obj.coeff{compound.drug_type}=compound.para_val;
            obj.dose=[0 0 0];
            obj.c_max=trial.c_max;
            obj.dosing_freq_hour=trial.dosing_freq_hour;

        end

        function iterate(obj,t)   %ez drug concentration
            t=t/60;
            effect=0;
            while t>=0                
                effect=effect+iterate_single(obj,t);
                t=t-obj.dosing_freq_hour;
            end
            obj.dose(obj.compound.drug_type)=effect;
        end


function conc=iterate_single(obj,t)   %ez drug concentration
    

    grad=obj.c_max/obj.compound.tmax_hour;
    
    k=log(2)/obj.compound.halfLife_hour;
    
    conc=(t<=obj.compound.tmax_hour)*t*grad+(t>obj.compound.tmax_hour)*obj.c_max*exp(-k*(t-obj.compound.tmax_hour));
    
end


    end

end