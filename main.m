%----------------------  adjustable parameters  ----------------------
iterations=1000;
unit=1; % miu m, space step
timestep=1e-3; %second, time step


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
d2=1400;    %diffusion coeff of hepatocyte     3600(not sure)

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



%----------------------  construction of matrices   ----------------------
% A c^(n+1)= B c^n+ b
C=zeros(total,1);     %vector of each grid
b=zeros(total,1);

sA=[]; %sparse matrix for iteration matrices A and B
sB=[]; 

%construction of A and B
corners=[1 up(1,sinusoid+2) up(1,radius+2) up(length+2,1) up(length+2,sinusoid+2) up(length+2,radius+2)];  %4 corners and 2 ends of the interface set to 0
for k=1:6
    n=corners(k);
    sA=[sA;n n 1];
end

for j=2:(sinusoid+1)   %sinusoid lower/upper boundary
     l=up(1,j); 
     m=up(2,j);
     n=up(3,j);

     p=up(length+2,j);  
     q=up(length+1,j);
     r=up(length,j);

     

     sA=[sA;l m 1;m m f*velo;m n -f*0.5*d1/unit;m l f*0.5*d1/unit;p p 1;q q 1;];
     sB=[sB;q q 1-f*velo;q r f*velo ];
     b(l)=oxy_pp;
     b(m)=f*flux;
end

for j=(sinusoid+3):(radius+1)   %hepatocyte lower/upper boundary, zero flux
     n=up(1,j);   %lower
     m=up(3,j);          
     p=up(length+2,j);  %higher
     q=up(length,j);

     sA=[sA;n n 1;n m -1;p p 1;p q -1];
end


for i=2:(length+1)    %inside
    n=up(i,1);  %sinusoid left boundary, symmetry
    sA=[sA;n n 1;n up(i,3) -1];

    for j=2:(sinusoid+1)  %sinusoid
       if and(i>2,i<length+1)
        n=up(i,j);
        sA=[sA;n up(i+1,j) -0.5*miu1;n up(i-1,j) -0.5*miu1;n up(i,j+1) -0.5*miu1;n up(i,j-1) -0.5*miu1;n n 1+2*miu1];
        sB=[sB;n up(i+1,j) 0.5*miu1;n up(i-1,j) 0.5*miu1+f*velo;n up(i,j+1) 0.5*miu1;n up(i,j-1) 0.5*miu1;n n 1-2*miu1-f*velo];
       end
    end

    j=sinusoid+2;    %interface
    n=up(i,j);
    sA=[sA; n n -miu1-miu2;n up(i,j+1) miu2;n up(i,j-1) miu1];

    
    for j=(sinusoid+3):(radius+1)  %hepatocytes
        n=up(i,j);
        sA=[sA;n up(i+1,j) -0.5*miu2;n up(i-1,j) -0.5*miu2;n up(i,j+1) -0.5*miu2;n up(i,j-1) -0.5*miu2;n n 1+2*miu2];
        sB=[sB;n up(i+1,j) 0.5*miu2;n up(i-1,j) 0.5*miu2;n up(i,j+1) 0.5*miu2;n up(i,j-1) 0.5*miu2;n n 1-2*miu2];
    end 
    n=up(i,radius+2);       %hepatocytes right boundary, zero flux
    sA=[sA;n n 1;n up(i,radius) -1];
end

%sparse matrix defined
A=sparse(sA(:,1),sA(:,2),sA(:,3),total,total);
B=sparse(sB(:,1),sB(:,2),sB(:,3),total,total);


%pre-processing to get the matrices for iterations
A_inv=inv(A);
update=A\B;



% ----------------------iterations ----------------------
C_previous=C; %two-step method for the reaction term
count=0;
while count<iterations
    b_current=b;  
    count    
    for i=2:(length+1)    % update the oxygen consumption term
        for j=(sinusoid+2):(radius+1)
            n=up(i,j);
            b_current(n)=0.5*timestep*(MM(C_previous(n))-3*MM(C(n)));            
        end
    end

    result= update*B*C+A_inv*b_current; %iterate
    C_previous=C;   %update
    C=result;
    count=count+1;
end

%---------------------- visualisation ----------------------

%map back to 2-d
res=matrix_assemble(C);
save('result.mat','res')
%plot
%heatmap(res);




%---------------------- helper functions ----------------------

function [res]=up(i,j)  %matrix index to vector index
    global radius;
    res=(i-1)*(radius+2)+j;
end

function [i,j]=down(k) %vector index to matrix index
    global radius;
    r=radius+2;
    j=mod(k,r);
    if j==0
        j=r;
    end
    i=(k-j)/r+1;

end

%vector-matrix conversion

function [c]=matrix_assemble(C)
    global length;
    global radius;
    global total;
    c=zeros(radius+2,length+2);
    for n=1:total
        c(n)=C(n);
    end
end

%Michalis-Menton
function [res]=MM(c)
    global v_max;
    global k_max;
    res=c*v_max/(c+k_max);
end
