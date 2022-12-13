classdef Method
        
    properties (Access=private)
        par
    end

    methods 
        function obj=Method(para)
            obj.par=para;
        end
        function [res]=MM(obj,c,v,k)
            res=c.*v./(c+k);
        end
        %{
            function res=HILL(obj,c,k,n)
                temp=c.^n;
        res=temp./(temp+k.^n);
            end
        %}
            function res=HillNorm(obj,c,basal,k,n)
                temp=c.^n;
        res=(temp+basal.^n)./(temp+k.^n);
        end

function [res]=up(obj,i,j)  %matrix index to vector index
    res=(i-1)*(obj.par.radius+2)+j;
end


function [c]=matrix_assemble(obj,C)
    c=zeros(obj.par.radius+2,obj.par.length+2);
    for n=1:obj.par.total
        c(n)=C(n);
    end
end
    

        function cel=parti(obj,list,batch)
        cel={};
        m=numel(list);
        lef=mod(m,batch);
        quo=(m-lef)/batch;
        st=1;
        for l=1:lef
            ed=st+quo;
            cel{l}=list(st:ed);
            st=ed+1;
        end

        for l=(lef+1):batch
            ed=st+quo-1;
            cel{l}=list(st:ed);
            st=ed+1;
        end
    end

    end
end