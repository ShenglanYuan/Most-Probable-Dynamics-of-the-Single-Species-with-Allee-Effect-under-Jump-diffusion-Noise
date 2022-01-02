
%dX = rX(1-X/T)(1-X/K)-lambdaX.,  beta=K/T,
%%0 < lambda/r < ((1+beta)^2/4beta  -1)
clear all;
x1=0; x2=8; dx=1/2^10;
x=x1:dx:x2;
a1=1; a2=0.1; a3=2.67; a4=1/a3;   % parameter values
D=@(x) a1*x-a2*x.^2-a3*(x./(a3*a4*x+1));  %deterministic term D(x)=x(Ax^2+Bx+C), 
A=a2*a3*a4; B=(a1*a3*a4-a2); C=a3;%(a1-a3);
dD=@(x) a1-2*a2*x-a3/(a3*a4*x+1)^2;
% Equilibrium points 
%%% 0 < a3 < (a1+a2)^2/(4*a2) or 0< a3 < 3.025; and 0 < a2 < 1
Bp=(a1*a3*a4+a2)^2-(4*a2*a3);
y1=0;
y2=(B-sqrt(Bp))/(2*A); y3=(B+sqrt(Bp))/(2*A); y4= B/(2*A);
D(y2);
D(y3);
%Bp=(a1*a3*a4-a2)^2-(a2*a3*a4*(a3-a1)) % bifurcation point must be greater to 4
%%%-----------Potential function
U=@ (x)-0.5*a1*x.^2+(a2/3)*(x.^3)+(a3/(a3*a4)^2)*(((a3*a4)*x+1)-log((a3*a4)*x+1)); % potential fun
aC=(U(y1)+0.5*a1*y3^2-(a2/3)*(y3^3))/(((y3+1)-log((y3+1))));  % U(y1)=U(y3) crtical value

M=round((x2-x1)/dx);
X=zeros(1,M+1);
for i=1:M+1
 X(i)=U(x(i));
end

plot(x,X,'-')

%ylabel('dx/dt')
%ylabel('U(x)')
xlabel('x')
hold on
box off
grid off