%parameters
oxy_pp=65;
oxy_pv=32;
velo=1;    %advection velocity

v=10;   %mm-equation
k=5;

sinusoid=5;
radius=sinusoid*5;
length=sinusoid*100;
total=length*radius;

timestep=1;
spacestep=1;
d1=1;   
d2=1;    %diffusion coeff of sinusoid and hepatocyte
miu1=timestep*d1/(spacestep^2);
miu2=timestep*d2/(spacestep^2);

%initialisation of matrix    A c^(n+1)= B c^n+ b

c=zeros(length,radius);
c(1,:)=oxy_pp;       %boundary conditions
c(length,:)=oxy_pv;
C=vector_assemble(c);
C_previous=vector_assemble(c);

A=zeros(total);
B=zeros(total);
b=zeros(total,1);




%iteration matrix construction

for i=2:(length-1)   
    n=up(i,1);  %axis condition: C_r=0 
    A(n,n)=1;
    A(n,n+1)=-1;    

    for j=2:(sinusoid-1)  %sinusoid
        n=up(i,j);
        row=-0.5*miu1*CN(i,j)+0.25*velo*timestep/spacestep*advect(i,j); % Crank-Nicolson with advection
        A(n,:)=row;
        A(n,n)=A(n,n)+1;   %identity matrix
        B(n,:)=-row;
        B(n,n)=B(n,n)+1;       
    end

    j=sinusoid;    %interface
    n=up(i,j);
    A(n,n)=-d1-d2;
    A(n,up(i,j+1))=d2;
    A(n,up(i,j-1))=d1;
    
    for j=(sinusoid+1):(radius-1)  %hepatocytes
        n=up(i,j);
        row=-0.5*miu1*CN(i,j); % Crank-Nicolson 
        A(n,:)=row;
        A(n,n)=A(n,n)+1;   %identity matrix
        B(n,:)=-row;
        B(n,n)=B(n,n)+1;   

      
    end 
end

count=0;
while count<10000


    %  b(n)=timestep*(1.5*MM(C(n))-0.5*MM(C_previous(n)));


    C_previous=C;   %update
    count=count+1;
end



%plot
heatmap(x);








%functions

%vector-matrix conversion
function [res]=up(i,j)
    global radius;
    res=(i-1)*radius+j;
end

function [i,j]=down(k)
    global radius;
    j=mod(k,radius);
    if j==0
        j=radius;
    end
    i=(k-j)/radius;

end

function [C]=vector_assemble(c)
    global total;
    C=zeros(total,1);
    for n=1:total
        [i,j]=down(n);
        C(n)=c(i,j);
    end
end

function [c]=matrix_assemble(C)
    global length;
    global radius;
    c=zeros(length,radius);
    for i=1:length
        for j=1:radius
            c(i,j)=C(up(i,j),1);
        end
    end
    

end

%Michalis-Menton
function [res]=MM(c,v,k)
    res=c*v/(c+k);
end

%crank-nicolson component
function [res]=CN(i,j)
    global total;
    res=zeros(1,total);
    res(up(i+1,j))=1;
    res(up(i-1,j))=1;
    res(up(i,j))=-4;
    res(up(i,j+1))=1+0.5/j;
    res(up(i,j-1))=1-0.5/j;
end

%advection component
function [res]=advect(i,j)
    global total;
    res=zeros(1,total);
    res(up(i+1,j))=1;
    res(up(i-1,j))=-1;
end
