classdef Mitochondria < handle
    properties (Access=private)
        par
        method
        option
        listHepato
    end
    properties
        M
        states
    end
    
    methods 
        function obj=Mitochondria(para,option,simple,listHepato) % simple==1: 3-zone model for comparison
            
            obj.par=para;
            obj.method=Method(para);
            obj.option=option;
            if option==1     %with density
                   obj.states=[3 4];  %[mitochondria implicit states, explicit states]
            else
                     obj.states=[3 3];  %[mitochondria implicit states, explicit states]
            end
            if simple==0
                obj.listHepato=listHepato;
            obj.M=zeros(para.total,sum(obj.states)); 
            obj.M(listHepato,:)=repmat(para.basal(1:sum(obj.states)),numel(listHepato),1);
            else
                obj.listHepato=[1 2 3];
                obj.M=repmat(para.basal(1:sum(obj.states)),3,1);
            end
        end

        function out=implicit(obj,y,drug)
                
    para=[cell2mat(drug.coeff) obj.par.para_mito];
    drug=drug.dose;
            basal=obj.par.basal;
            HillNorm=@(x,y,z,w)obj.method.HillNorm(x,y,z,w);
            MM=@(x,y,z)obj.method.MM(x,y,z);

       drug_atpase=MM(para(4),1,drug(3));
       k_atp=basal(5)/(basal(5)*para(7)-1); 
       out=zeros(size(y,1),3);
       out(:,1)=basal(1)*HillNorm(para(5),basal(4),y(:,1),para(6)).*HillNorm(para(8),basal(5),y(:,2),para(9));
       out(:,2)=basal(2)*MM(y(:,2),para(7)*k_atp,k_atp).*MM(basal(6),2,y(:,3))*drug_atpase;
       out(:,3)=basal(3)*HillNorm(para(10),basal(6),y(:,3),para(11));

    end


    function out=mito_sys(obj,t,y,c,drug_struct)  %para=[k_etc   v_max_uncp k_uncp    k_atpase      k_mito  n_mito  v_max_atp ]  drug=[etc uncp atpase]
    HillNorm=@(x,y,z,w)obj.method.HillNorm(x,y,z,w);
    MM=@(x,y,z)obj.method.MM(x,y,z);
    
    
    para=[cell2mat(drug_struct.coeff) obj.par.para_mito];
    drug=drug_struct.dose;
    basal=obj.par.basal;
    coeff=obj.par.coeff;
    para2=obj.par.para_den;

    y=reshape(y,[],obj.states(2));


    drug_etc=MM(para(1),1,drug(1));
    drug_uncp=y(:,2)*MM(drug(2),para(2)*para(3),para(3));
    ox2ocr=(drug_etc/basal(4))*MM(c,obj.par.v_respire,obj.par.k_respire);
    ox2etc=ox2ocr/coeff(1);
    
    atp_dec=1-(y(:,3)/basal(6));
    atp_dec=(atp_dec>=0).*atp_dec;   
    z=obj.implicit(y,drug_struct);
    %y(:,1)*ox2etc
    %z(:,1)
    out=[-y(:,1).*ox2etc+z(:,1) coeff(2)*y(:,1).*ox2etc-coeff(3)*z(:,2)-drug_uncp coeff(4)*(z(:,2)+z(:,3))-coeff(5)*y(:,3)];

    if obj.option==1
        den=-(para2(1)*HillNorm(atp_dec,0,para2(2),para2(3))+para2(4)*HillNorm(atp_dec,0,para2(5),para2(6))).*y(:,4);
        den=den+0.00045*y(:,4).*(1-y(:,4));
        out=[out den];
    end
    out=out(:);


end


function iterate(obj,c,drug)
            colomn_impli=1:obj.states(1);
            colomn_expli=(obj.states(1)+1):sum(obj.states);

            p=obj.method.parti(obj.listHepato,obj.par.batch_number(2));
            for q=1:numel(p)
                activ=p{q};
                if ~isempty(activ)
                temp=obj.M(activ,colomn_expli);
            [t,y]=ode15s(@(t,y)obj.mito_sys(t,y,c(activ),drug),[0 obj.par.time(2)/60],temp(:));
            
            last=reshape(y(size(y,1),:),[],obj.states(2));
             obj.M(activ,colomn_expli)=last;
             obj.M(activ,colomn_impli)=obj.implicit(last,drug);
                end
            end

             %save('y_new.mat',"y");
             %save('new.mat',"mitocon");

end




    end


end