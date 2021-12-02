%----------------------  adjustable parameters  ----------------------
iterations=5;
unit=0.5; % miu m, space step
timestep=1; %second, tine step


%----------------------   parameters ----------------------

%diameter of hepatocyte rescaled, 25-30 miu m
%diameter of sinusoid rescaled, 10-15 miu m
%length of sinusoid rescaled, 275 miu m
hepatocyte_diamter=27.5;   
sinusoid_diameter=12.5;
sinusoid_length=275;


sinusoid=round(sinusoid_diameter/2/unit);
global radius
radius=round(sinusoid+hepatocyte_diamter/unit)
global length
length=round(sinusoid_length/unit)
global total
total=(length+2)*(radius+2)    %enlarge the entire matrix by 1 unit, for systematic treatment of boundary conditions


%boundary conditions: oxygen concentration and flux at entry (TBD)
oxy_pp=65;   
flux=10;  

%advection and diffusion
flow_rate=0.9;   %8-10 mm/sec?
d1=5;     % diffusion coeff of sinusoid
d2=0.4;    %diffusion coeff of hepatocyte

velo=flow_rate*1000/unit;   
miu1=timestep*d1;
miu2=timestep*d2;


 %  Michalis-Menton equation, metabolism for oxygen
global v_max
v_max=10;  
global k_max
k_max=5;



%----------------------  construction of matrices   ----------------------
% A c^(n+1)= B c^n+ b

c=zeros(length+2,radius+2);
c(2,2:sinusoid)=oxy_pp;       %Dirichlet boundary condition for sinusoid
C=vector_assemble(c);     %vector of each grid

b=zeros(total,1);

sA=[]; %sparse matrix for iteration matrices A and B
sB=[]; 

%construction of A and B
corners=[1 up(1,sinusoid+1) up(1,radius+2) up(length+2,1) up(length+2,sinusoid+1) up(length+2,radius+2)];  %4 corners and 2 ends of the interface set to 0
for k=1:6
    n=corners(k);
    sA=[sA;n n 1];
end

for j=2:sinusoid   %sinusoid lower/upper boundary, Neumann boundary conditions
     n=up(1,j);    %lower: flux
     m=up(3,j);
     l=up(length+2,j);  % upper: to be fixed
     sA=[sA;n n 1;n m -1;l l 1];
     b(n)=-2*flux;

end

for j=(sinusoid+2):(radius+1)   %hepatocyte lower/upper boundary, zero flux
     n=up(1,j);   %lower
     m=up(3,j);          
     p=up(length+2,j);  %higher
     q=up(length,j);

     sA=[sA;n n 1;n m -1;p p 1;p q -1];
end


for i=2:(length+1)    %inside
    n=up(i,1);  %sinusoid left boundary, symmetry
    sA=[sA;n n 1;n (n+2) -1];

    for j=2:sinusoid  %sinusoid
        n=up(i,j);
        sA=[sA;n up(i+1,j) -0.5*miu1;n up(i-1,j) -0.5*miu1;n up(i,j+1) -0.5*miu1;n up(i,j-1) -0.5*miu1;n n 1+2*miu1];
        sB=[sB;n up(i+1,j) 0.5*miu1;n up(i-1,j) 0.5*miu1+timestep*velo;n up(i,j+1) 0.5*miu1;n up(i,j-1) 0.5*miu1;n n 1-2*miu1-timestep*velo];

    end

    j=sinusoid+1;    %interface
    n=up(i,j);
    sA=[sA; n n -d1-d2;n up(i,j+1) d2;n up(i,j-1) d1];

    
    for j=(sinusoid+2):(radius+1)  %hepatocytes
        n=up(i,j);
        sA=[sA;n up(i+1,j) -0.5*miu2;n up(i-1,j) -0.5*miu2;n up(i,j+1) -0.5*miu2;n up(i,j-1) -0.5*miu2;n n 1+2*miu2];
        sB=[sB;n up(i+1,j) 0.5*miu2;n up(i-1,j) 0.5*miu2;n up(i,j+1) 0.5*miu2;n up(i,j-1) 0.5*miu2;n n 1-2*miu2];
    end 
    n=up(i,radius+2);       %hepatocytes right boundary, zero flux
    sA=[sA;n n 1;n (n-2) -1];
end

%sparse matrix defined
A=sparse(sA(:,1),sA(:,2),sA(:,3),total,total);
B=sparse(sB(:,1),sB(:,2),sB(:,3),total,total);

%save('matrices.mat')

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
%plot
heatmap(res);




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
    i=(k-j)/r;

end

%vector-matrix conversion

function [C]=vector_assemble(c)
    global total;
    C=zeros(total,1);
    for n=1:total
        C(n)=c(n);
    end
end

function [c]=matrix_assemble(C)
    global length;
    global radius;
    global total;
    c=zeros(length+2,radius+2);
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
