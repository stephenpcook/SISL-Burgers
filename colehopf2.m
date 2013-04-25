clear all
N=99;
nabla = 2; 
theta = 1/2+0.1; 
epsilon = 1/200;

dx = 1/(N+1);
dt = nabla*(dx^2);

x = (0:(N+1))'./(N+1);
rho0 = exp(1/(4*epsilon*pi)*(cos(2*pi*x) + cos(pi*x) - 2));

A = -2*eye(N+2);
for i = 1:N
    A(i+1,i)=1; 
    A(i,i+1)=1;
end
A(1,2)=2; 
A(N+2,N+1)=2;

rho = rho0;

u = -2*epsilon*([0;1/2*(diff(rho(1:N+1))+diff(rho(2:N+2)));0]./dx)./(rho);
subplot(2,1,1)
plot(x,rho)
title(['t = 0'])
drawnow()
subplot(2,1,2)
plot(x, u)
hold on;
plot(x, sin(2*pi*x) + 1/2*sin(pi*x),'-k')
hold off
pause

for t = 0:dt:1
    RHS_n= (eye(N+2)+epsilon*nabla*(1-theta)*A)*rho;
    M = eye(N+2) - epsilon*nabla*theta*A;
    Ma= [0;diag(M,-1)]; Mb= diag(M); Mc= diag(M,1);
    %rho= M\RHS_n;
    rho = TDMAsolver(Ma,Mb,Mc,RHS_n)';
    
    u = -2*epsilon*([0;diff(rho)]./dx)./(rho);
    subplot(2,1,1)
    plot(x,rho)
    title(['t =',num2str(t)])
    drawnow()
    subplot(2,1,2)
    plot(x, u)
    hold on;
    plot(x, sin(2*pi*x) + 1/2*sin(pi*x))
    hold off
    drawnow()
end
