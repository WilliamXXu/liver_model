%----------------------  adjustable parameters  ----------------------
iterations=1.5e4;
unit=1; % miu m, space step
timestep=1e-4; %second, time step


%----------------------   parameters ----------------------

%diameter of hepatocyte rescaled, 25-30 miu m
%diameter of sinusoid rescaled, 10-15 miu m
%length of sinusoid rescaled, 275 miu m

hepatocyte_diamter=27.5;   
sinusoid_diameter=12.5;
sinusoid_length=275;
%{
hepatocyte_diamter=10;   
sinusoid_diameter=10;
sinusoid_length=20;
%}

sinusoid=round(sinusoid_diameter/2/unit);
global radius
radius=round(sinusoid+hepatocyte_diamter/unit)
global length
length=round(sinusoid_length/unit)
global total
total=(length+2)*(radius+2)    %enlarge the entire matrix by 1 unit, for systematic treatment of boundary conditions


%advection and diffusion   
velo=400;   %miu/sec?
f=timestep/unit;   %check courant for upwind method'
if abs(f*velo)>1
    'courant number>1. upwind method diverges. quit'
    return
end

d1=1620;     % diffusion coeff of sinusoid    miu m^2/sec
d2=1400 ;    %diffusion coeff of hepatocyte     3600(not sure)

miu1=timestep*d1/(unit*unit);
miu2=timestep*d2/(unit*unit);

%boundary conditions: oxygen concentration and flux at entry
oxy_pp=65;     %mmHg 
conversion=1.34e-3;    %to mol/m3, Henry's law
flux=oxy_pp*600;     %mmHg*miu/sec, slightly higher than velo*oxy_pp

 %  Michalis-Menton equation, metabolism for oxygen
global v_max
v_max=0.044/conversion;  %mmHg/sec  
global k_max
k_max=6.24e-3/conversion; %mmHg

MM(68)
%Michalis-Menton
function [res]=MM(c)
    global v_max;
    global k_max;
    res=c*v_max/(c+k_max);
end
