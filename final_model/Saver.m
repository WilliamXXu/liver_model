classdef Saver < handle
    properties (Access=private)
    method
    par
    simple
    end
   properties
   id
   frequency
   directory
   end
    methods 
        function obj=Saver(par,simple)
            obj.method=Method(par);
            obj.par=par;
            obj.simple=simple;
        end
        function makeFolder(obj)
            obj.directory=append('data/',num2str(obj.id));
            lastwarn('')
            mkdir(obj.directory);
            if strcmp(lastwarn,'Directory already exists.')
                error('Check ID to prevent override.')
            end
        end

        function snapshot(obj,count,oxygen,mitochondria)
            if obj.simple==0
        matrix_assemble=@(x)(obj.method.matrix_assemble(x));
        
        length=obj.par.length;
        radius=obj.par.radius;
         if mod(count,obj.frequency)==0
                 oxy=matrix_assemble(oxygen.C);
                 hemo=matrix_assemble(oxygen.H);
                 %-oxygen.b_current/obj.par.time(1)
                 resp=matrix_assemble(-oxygen.b_current/obj.par.time(1));
%remove boundry units
oxy=oxy(2:radius+1,2:length+1);
hemo=hemo(2:radius+1,2:length+1);
resp=resp(2:radius+1,2:length+1);

for ind=1:sum(mitochondria.states)
    res=matrix_assemble(mitochondria.M(:,ind));
    res=res(2:radius+1,2:length+1);
    res=res/obj.par.basal(ind);
    save(append(obj.directory,'/mito',num2str(ind),'_',num2str(count),'.mat'),"res");
end

%save('drug.mat','sol');
save(append(obj.directory,'/oxy','_',num2str(count),'.mat'),"oxy");
save(append(obj.directory,'/hemo','_',num2str(count),'.mat'),'hemo');
save(append(obj.directory,'/resp','_',num2str(count),'.mat'),'resp');
        end
        
            else
                for ind=1:sum(mitochondria.states)
                    res=mitochondria.M(:,ind);
                    res=res/obj.par.basal(ind);
                    save(append(obj.directory,'/mito',num2str(ind),'_',num2str(count),'.mat'),"res");
                end
            end

    end
    end
end