
%function dru()


para=[log(2)/500 log(2)/510 log(2)/800 log(2)/400 ];


tspan=[0 1200];
initial=[0 0];

sol=ode15s(@(t,y)drug_sys(t,y,para,@input),tspan,initial);
deval(sol,1)

% x=linspace(0,tspan(2),size(y,1));
% 
% plot(t,y(:,2),'-x')
% hold on
% plot(t,y(:,1),'-x')
% hold off
function out=input(t)
    thresh=120;
    if or(t<0,t>thresh)
        out=0;
    else
        out=0.01/thresh;
    end
end

%end