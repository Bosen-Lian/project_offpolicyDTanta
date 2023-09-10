clear all;
close all;
clc;


%系统
 A=[   0.9933    0.3645    0.0279   -0.0066;
   -0.0244    0.8295    0.1192   -0.0480;
   -0.2551   -0.0533    0.5017   -0.4950;
    0.0299    0.0056    0.0003    0.9999];


B=[ 0.0066;
    0.0480;
    0.4950;
    0.0001];


D=[ -0.0062;
    0.0001;
    0.0009;
   -0.0001];


% A=[0.9704  0.6629  0.0849 -0.0446;
%   -0.0762  0.6724  0.1584 -0.1462;
%   -0.3954  -0.1663 0.2367 -0.7403;
%    0.0594  0.0212  0.0019 0.9993];
% 
% 
% B=[0.0446;
%     0.1462;
%     0.7403;
%     0.0007];
% 
% D=[-0.7924;
%     0.0230;
%     0.1893;
%     -0.0239]; 



%专家参数
Qe=1*[1 0 0 0;
      0 1 0 0;
      0 0 1 0;
      0 0 0 1];
Re=1;
gammae=5;
Pe=dare(A,[B,D],Qe,[Re 0; 0, -gammae^2]);
Ke=inv(Re+B'*Pe*B+B'*Pe*D*inv(gammae^2-D'*Pe*D)*D'*Pe*B)*(B'*Pe*A+B'*Pe*D*inv(gammae^2-D'*Pe*D)*D'*Pe*A);
Le=inv(-gammae^2+D'*Pe*D-D'*Pe*B*inv(Re+B'*Pe*B)*B'*Pe*D)*(D'*Pe*A-D'*Pe*B*inv(Re+B'*Pe*B)*B'*Pe*A);



%学生参数
Q=0*[0.01 0 0 0;
      0 0.01 0 0;
      0 0 0.01 0;
      0 0 0   0.01];
R=4;
gamma=10;
 K1=[1.8152    0.3844    1.3118    0.4207];
% L1=[-0.0477 -0.0244 0.0001 -0.0123];
L1=[0 0 0 0];
K=K1;
L=L1;



%迭代初值
x(1:4,1)=[2;-2;2;-2];
xx(1:4,1)=[2;4;-4;-3];
k=0;
kk=0;
ii=0;
jj=0;
eK=norm(K1-Ke);
eL=[];
eQ=[];
stopq=0;
alpha=0.1;
stop=0.01;
uq=0;

for i=1:2000
    
    %专家数据
    d(i)=0.0005*(rand(1));
    ue(:,i)=-Ke*x(:,i)+d(i);
%       we(:,i)=-Le*x(:,i)+0.010*(rand(1));
      we(:,i)=0.0010*(rand(1));
    x(:,i+1)=A*x(:,i)+B*ue(:,i)+D*we(:,i);
    
      uue(:,i)=-Ke*xx(:,i);
%         wwe(:,i)=-Le*xx(:,i)+0.0000*(rand(1));
        wwe(:,i)=0.00010*(rand(1));
    xx(:,i+1)=A*xx(:,i)+B*uue(:,i)+D*wwe(:,i);  
    
    
    
    ii=ii+1;
    Hxx(ii,:)=[x(1,i)^2, 2*x(1,i)*x(2,i), 2*x(1,i)*x(3,i), 2*x(1,i)*x(4,i), x(2,i)^2, 2*x(2,i)*x(3,i), 2*x(2,i)*x(4,i),x(3,i)^2, 2*x(3,i)*x(4,i), x(4,i)^2]...
              -[x(1,i+1)^2, 2*x(1,i+1)*x(2,i+1), 2*x(1,i+1)*x(3,i+1), 2*x(1,i+1)*x(4,i+1), x(2,i+1)^2, 2*x(2,i+1)*x(3,i+1), 2*x(2,i+1)*x(4,i+1),x(3,i+1)^2, 2*x(3,i+1)*x(4,i+1),x(4,i+1)^2];%1-10  
    Hxu(ii,:)=2*kron(x(:,i+1)'*K1',x(:,i+1)')+2*kron(ue(:,i)',x(:,i)');%11-14
    Hxw(ii,:)=2*kron(x(:,i+1)'*L1',x(:,i+1)')+2*kron(we(:,i)',x(:,i)');%15-18
    Huu(ii,:)=kron(ue(:,i)',ue(:,i)')-kron(x(:,i+1)'*K1',x(:,i+1)'*K1');%19
    Huw(ii,:)=2*kron(we(:,i)',ue(:,i)')-2*kron(x(:,i+1)'*L1',x(:,i+1)'*K1');%20
    Hww(ii,:)=kron(we(:,i)',we(:,i)')-kron(x(:,i+1)'*L1',x(:,i+1)'*L1');%21   
    r(ii)=x(:,i)'*K1'*R*K1*x(:,i)+x(:,i)'*Q*x(:,i)+alpha*(ue(:,i)+K1*x(:,i)-d(i))'*R*(ue(:,i)+K1*x(:,i)-d(i))...
          -gamma^2*(we(:,i)-L1*x(:,i))'*(we(:,i)+L1*x(:,i))-gamma^2*x(:,i)'*(L1'*L1)*x(:,i)+(ue(:,i)-K1*x(:,i))'*R*(ue(:,i)+K1*x(:,i));
    zbar=[Hxx Hxu Hxw Huu Huw Hww];
    
    if stopq==0
        jj=jj+1;
       X(jj,:)=[x(1,i)^2, 2*x(1,i)*x(2,i), 2*x(1,i)*x(3,i), 2*x(1,i)*x(4,i), x(2,i)^2, 2*x(2,i)*x(3,i), 2*x(2,i)*x(4,i), x(3,i)^2, 2*x(3,i)*x(4,i), x(4,i)^2];
       Hqxx(jj,:)=Hxx(ii,:);%kron(x(:,i)',x(:,i)')-kron(x(:,i+1)',x(:,i+1)');
       Hqxu(jj,:)=Hxu(ii,:);
       Hqxw(jj,:)=Hxw(ii,:);
       Hquu(jj,:)=Huu(ii,:);
       Hquw(jj,:)=Huw(ii,:);
       Hqww(jj,:)=Hww(ii,:);
       rq(jj)=r(ii)-x(:,i)'*Q*x(:,i)-alpha*(ue(:,i)-d(i)+K1*x(:,i))'*R*(ue(:,i)-d(i)+K1*x(:,i));

% Hq(jj)=1*xx(:,i)'*Q*xx(:,i)+alpha*(uue(:,i)+K1*xx(:,i))'*R*(uue(:,i)+K1*xx(:,i));


    end
    if rank(X)==10
        stopq=1;
    end
    
    
%     m=zbar'*zbar; 
%     m=zbar; 
%     q=zbar'*r';
    a=rank(zbar)
        
    if a==21
        
       H=zbar\r';

       K1=inv(H(19)-H(20)*inv(H(21))*H(20))*([H(11) H(12) H(13) H(14)]-H(20)*inv(H(21))*[H(15) H(16) H(17) H(18)]);
       L1=inv(H(21)-H(20)*inv(H(19))*H(20))*([H(15) H(16) H(17) H(18)]-H(20)*inv(H(19))*[H(11) H(12) H(13) H(14)]);
       Hq=Hqxx*[H(1),H(2),H(3),H(4),H(5),H(6) H(7) H(8) H(9) H(10)]'...
            +Hqxu*[H(11) H(12) H(13) H(14)]'+Hqxw*[H(15) H(16) H(17) H(18)]'+Hquu*H(19)+Hquw*H(20)+Hqww*H(21)-rq;

       
       K=[K;K1];
       L=[L,L1];
       
       eK=[eK,norm(K1-Ke)];
       eL=[eL,norm(L1-Le)];    

        Hxx=[];
        Hxu=[];
        Hxw=[];
        Huu=[];
        Huw=[];
        Hww=[];
        zbar=[];
        r=[];
        ii=0;
        m=[];
        
%         Q=Q+0.7*(Ke-K1)'*R*(Ke-K1);
%            eQ=[eQ,norm(Q)];
               % Hq=rq;
       
       
        rank(X)
        HQ=X\Hq;
        if eK(end)>=stop        
        Q=[HQ(1),HQ(2),HQ(3),HQ(4);
           HQ(2),HQ(5),HQ(6),HQ(7);
           HQ(3),HQ(6),HQ(8),HQ(9);
           HQ(4),HQ(7),HQ(9),HQ(10)];
        eQ=[eQ,norm(Q)];
        end
        jj=0;
        stopq=0;
        X=[];
        Hqxx=[];
        Hqxu=[];
        Hqxw=[];
        Hquu=[];
        Hquw=[];
        Hqww=[];
        Hq=[];
        
    end
    
        if eK(end)< stop
        uq=1;
        end


end
 
figure(1)
subplot(2,1,1);plot(eQ,'-*','LineWidth',1)
xlabel('Iteration step $j$');
ylabel('$\Vert Q \Vert$');
set(get(gca,'XLabel'),'FontSize',18);
set(get(gca,'YLabel'),'FontSize',18);
set(get(gca,'XLabel'),'Interpreter','latex');set(get(gca,'YLabel'),'Interpreter','latex');
subplot(2,1,2);plot(eK,'-*','LineWidth',1)
xlabel('Iteration step $j$');ylabel('$\Vert K_j-K_e \Vert$');
set(get(gca,'XLabel'),'FontSize',18);
set(get(gca,'YLabel'),'FontSize',18);
set(get(gca,'XLabel'),'Interpreter','latex');set(get(gca,'YLabel'),'Interpreter','latex');




figure(2)
plot(eL,'-*','LineWidth',1)
legend('L-Le')


%%Trajectory
xt(1:4,1)=[1;-1;-1;1];
xxt(1:4,1)=[1;-1;-1;1];


for j=1:200
    
       ute(:,j)=-Ke*xt(:,j);
      wte(:,j)=0.0010*(rand(1));
    xt(:,j+1)=A*xt(:,j)+B*ute(:,j)+D*wte(:,j); 
    
          uute(:,j)=-K1*xxt(:,j);
      wwte(:,j)=0.0010*(rand(1));
    xxt(:,j+1)=A*xxt(:,j)+B*uute(:,j)+D*wte(:,j);  
    
    
end

j=0:1:200;
subplot(2,1,1);plot(j,xt(1,:),'--',j,xt(2,:),'--',j,xt(3,:),'--',j,xt(4,:),'--','LineWidth',2);
hold on 
plot(j,xxt(1,:),':',j,xxt(2,:),'k:',j,xxt(3,:),':',j,xxt(4,:),':','LineWidth',2);
xlabel('Time step');ylabel('States');
legend(' $ x_{k_1} $','$ x_{k_2}$','$ x_{k_3}$','$ x_{k_4}$',' $ \hat{x}_{k_1} $','$ \hat{x}_{k_2}$',' $ \hat{x}_{k_3} $','$ \hat{x}_{k_4}$');
set(get(gca,'XLabel'),'FontSize',18);
set(get(gca,'YLabel'),'FontSize',18);
set(get(gca,'XLabel'),'Interpreter','latex');set(get(gca,'YLabel'),'Interpreter','latex');
set(get(gca,'Legend'),'FontSize',20);set(legend,'Interpreter','latex');
subplot(2,1,2);plot(ute(1,:),'--','LineWidth',2);
hold on
plot(uute(1,:),':','LineWidth',2);
xlabel('Time step');ylabel('Control inputs');
legend(' $ u^*_{k} $',' $ \hat{u}_{k} $');
set(get(gca,'XLabel'),'FontSize',18);
set(get(gca,'YLabel'),'FontSize',18);
set(get(gca,'XLabel'),'Interpreter','latex');set(get(gca,'YLabel'),'Interpreter','latex');
set(get(gca,'Legend'),'FontSize',20);set(legend,'Interpreter','latex');



