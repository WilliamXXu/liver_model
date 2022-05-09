function out=drug_sys(t,y,para,input)   %para=[c to p, p to c, elimination, metabolism[

    out=[  input(t)-(para(1)+para(3))*y(1)+para(2)*y(2)  %central 
            para(1)*y(1)-(para(2)+para(4))*y(2)        %peripharal
    ];


end