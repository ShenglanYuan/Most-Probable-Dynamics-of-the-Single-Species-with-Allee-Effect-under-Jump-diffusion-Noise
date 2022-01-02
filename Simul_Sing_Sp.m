 %M?todo de Euler?Maruyama 
 %Generating 100 sample paths of stochastic logistic growth model 
 %dx=mudt+sigmadBW or where mu=r*x*(1-x/k), and sigma=lambda*x
 clear all 
 delta=2^(-10);              %X1=0, X2=2.6159; ,X3=6.3841
 T=20;   n=T/delta ;   X0=4.5; % initial pop  (X1<x0<X2,   XX2<x0<X3,  X3<x0)
 t0=0; %initial time
 num_trajetorys=1; 
 Y=zeros(num_trajetorys ,n) ; 
 Ye=zeros (1,n) ; 
 Yt=X0; 
 Yte=X0; 
 meantime=zeros(n,10) ; 
 extinction=0; 
 N=0.8; 
a1=1; a2=0.1; a3=2.67; a4=1/a3;   % parameter values
%D=@(x) a1*x-a2*x.^2-a3*(x./(a3*a4*x+1));  %deterministic term D(x)=x(Ax^2+Bx+C),
%X1=0, X2=2.6159; ,X3=6.3841
 hold on 
 for i=1:num_trajetorys 
     dW1=sqrt(delta)*randn(1,n) ; 
     %dW2=sqrt(delta)*randn(1,n) ; 
     Yt=X0; 
     for j=1:n
 mu=a1*Yt-a2*Yt.^2-a3*(Yt./(a3*a4*Yt+1)) ;  % deteminant logistic iwth Allee
  %effect
   lam=0.050; %noise intensity-Gaussian
   eps=0.0;  %Levy noise
  lambda=sqrt(abs(lam*Yt))/(2*sqrt(N)) ; %Gaussian
 % epsilon=sqrt(abs(eps*Yt))/(2*sqrt(N)) ; % non-Gaussian
  Yt=Yt+mu*delta+lambda*dW1(j) ;%random solution 
  Y(i , j)=Yt; 
     end
     bool=Y(i ,:) >0; 
     if (sum(bool)==n) 
         plot ([0: delta :T] ,[X0,Y(i ,:) ]) ;
     else
         meantime(: , extinction+1)=1-bool ; 
         extinction=extinction+1; 
         time=find(1-bool , 1 ) ; 
         plot ([0: delta :time*delta ] ,[X0,Y(i ,1: time) ] , '- ') ; 
     end
 end
%  for j=1:n 
%      mu=a1*Yt-a2*Yt.^2-a3*(Yt./(a3*a4*Yt+1)); % f(x)=rx(1-x/k), here Yte=x
%      Yte=Yte+mu*delta ; 
%      Ye(j)=Yte; % deterministic solution 
%  end
% plot ([0: delta :T] ,[X0,Ye],'-') ; 
 hold on 
 xlabel('time')
 ylabel('X_t')
 %legend('x_0=0.2','x_0=0.4', 'x_0=0.6','x_0=0.8','x_0=1.0','x_0=1.2','x_0=1.4','x_0=1.6','x_0=1.8','x_0=2.0','x_0=2.2','x_0=2.4','x_0=2.8','x_0=3.0')
 %legend('\epsilon=0.1','\epsilon=0.3','\epsilon=0.6','\epsilon=0.9')
 %legend('logistic model without Allee effect','logistic model with Allee effect')
 MeanTimeExtinction=find(meantime(: ,1) ,1)*delta ;
 grid off
 box off