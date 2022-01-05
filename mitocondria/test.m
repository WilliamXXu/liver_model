conversion=1.34e-3;
global v_max
%v_max=0.0154/conversion;  %mmHg/sec  
v_max=0.044/conversion;
global k_max
k_max=6.24e-3/conversion; %mmHg

syms c
    f(c)=c*v_max/(c+k_max);
g=finverse(f)

function [res]=MM(c)
    global v_max;
    global k_max;
    res=c*v_max/(c+k_max);
end