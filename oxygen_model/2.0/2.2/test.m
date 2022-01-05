conversion=1.34e-3;
global v_max
v_max=0.044/conversion;  %mmHg/sec  
global k_max
k_max=6.24e-3/conversion; %mmHg



x=arrayfun(@MM,[0:5:70])
plot(x)
function [res]=MM(c)
    global v_max;
    global k_max;
    res=c*v_max/(c+k_max);
end
