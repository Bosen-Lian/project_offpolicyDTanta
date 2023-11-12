clear all;
close all;
clc;


%
A=[0.995   0.229 ;
   -0.018   0.867];


B=[-0.416 0.004;
    0.004   0.03];

D=[1   0;
    0   1]; 



Qe=1*[1 0;
      0 0];
Re=eye(2);
gammae=5;

%%
Pe=dare(A,[B,D],Qe,[Re zeros(2); zeros(2), -gammae^2*eye(2)]);
Ke=inv(Re+B'*Pe*B+B'*Pe*D*inv(gammae^2*eye(2)-D'*Pe*D)*D'*Pe*B)*(B'*Pe*A+B'*Pe*D*inv(gammae^2*eye(2)-D'*Pe*D)*D'*Pe*A);
Le=inv(-gammae^2*eye(2)+D'*Pe*D-D'*Pe*B*inv(Re+B'*Pe*B)*B'*Pe*D)*(D'*Pe*A-D'*Pe*B*inv(Re+B'*Pe*B)*B'*Pe*A);




Q=0*[0.01 0 ;
      0 0.01];
R=2*eye(2);
gamma=10;

 K1=[-1.0152    -0.844;
     0.3118    0.4207];
 L1=[0.00 0.00;
     0.000 0.0];
 
K=K1;
L=L1;



%%
x(1:2,1)=[1;-1];
xx(1:2,1)=[1;-1];
k=0;
kk=0;
ii=0;
jj=0;
eK=0.1;
eL=[];
eQ=[];
stopq=0;
alpha=0.9;

for i=1:10000
    

    d(:,i)=0.001*[rand(1); rand(1)];
    ue(:,i)=-Ke*x(:,i)+d(:,i);
      we(:,i)=0.00010*[rand(1); rand(1)];
    x(:,i+1)=A*x(:,i)+B*ue(:,i)+D*we(:,i);
    
      uue(:,i)=-Ke*xx(:,i);
       wwe(:,i)=0.0010*[rand(1); rand(1)];
    xx(:,i+1)=A*xx(:,i)+B*uue(:,i)+D*wwe(:,i);  
    
    
    
    ii=ii+1;
    u1(:,i)=K1*x(:,i+1);
    w1(:,i)=L1*x(:,i+1);
    
    Hxx(ii,:)=[x(1,i)^2, 2*x(1,i)*x(2,i), x(2,i)^2]...
              -[x(1,i+1)^2, 2*x(1,i+1)*x(2,i+1), x(2,i+1)^2];%1-3  
    Hxu(ii,:)=2*kron(x(:,i+1)',x(:,i+1)'*K1')+2*kron(x(:,i)',ue(:,i)');%4-7
    Hxw(ii,:)=2*kron(x(:,i+1)',x(:,i+1)'*L1')+2*kron(x(:,i)',we(:,i)');%8-11
    Huu(ii,:)=[ue(1,i)^2 2*ue(1,i)*ue(2,i) ue(2,i)^2]-[u1(1,i)^2 2*u1(1,i)*u1(2,i) u1(2,i)^2];    
    Huw(ii,:)=2*kron(ue(:,i)',we(:,i)')-2*kron(x(:,i+1)'*K1',x(:,i+1)'*L1');%15-18   
    Hww(ii,:)=[we(1,i)^2 2*we(1,i)*we(2,i) we(2,i)^2]-[w1(1,i)^2 2*w1(1,i)*w1(2,i) w1(2,i)^2];
    r(ii)=x(:,i)'*K1'*R*K1*x(:,i)+x(:,i)'*Q*x(:,i)-gamma^2*x(:,i)'*(L1'*L1)*x(:,i)-gamma^2*(we(:,i)-L1*x(:,i))'*(we(:,i)+L1*x(:,i))...
        +alpha*(ue(:,i)+K1*x(:,i)-d(:,i))'*R*(ue(:,i)+K1*x(:,i)-d(:,i))+(ue(:,i)-K1*x(:,i))'*R*(ue(:,i)+K1*x(:,i));
    
    jj=jj+1;
    if stopq==0
       X(jj,:)=[x(1,i)^2, 2*x(1,i)*x(2,i), x(2,i)^2];
       XX(jj,:)=[xx(1,i)^2, 2*xx(1,i)*xx(2,i), xx(2,i)^2];
       Hqxx(jj,:)=Hxx(ii,:);
       Hqxu(jj,:)=Hxu(ii,:);
       Hqxw(jj,:)=Hxw(ii,:);
       Hquu(jj,:)=Huu(ii,:);
       Hquw(jj,:)=Huw(ii,:);
       Hqww(jj,:)=Hww(ii,:);
  rq(jj)=xx(:,i)'*Q*xx(:,i)+alpha*(uue(:,i)+K1*xx(:,i))'*R*(uue(:,i)+K1*xx(:,i));


    end
    if rank(X)==3
        stopq=1;
    end
    
    zbar=[Hxx Hxu Hxw Huu Huw Hww];
    m=zbar; 
    a=rank(m);
        
    if a==21
        
       H=m\r';

       
       K1=inv([H(12) H(13);H(13) H(14)]-[H(15) H(16);H(17) H(18)]*inv([H(19) H(20);H(20) H(21)])*[H(15) H(16);H(17) H(18)]')...
           *([H(4) H(5);H(6) H(7)]'-[H(15) H(16);H(17) H(18)]*inv([H(19) H(20);H(20) H(21)])*[H(8) H(9);H(10) H(11)]');
       L1=inv([H(19) H(20);H(20) H(21)]-[H(15) H(16);H(17) H(18)]'*inv([H(12) H(13);H(13) H(14)])*[H(15) H(16);H(17) H(18)])...
           *([H(8) H(9);H(10) H(11)]'-[H(15) H(16);H(17) H(18)]'*inv([H(12) H(13);H(13) H(14)])*[H(4) H(5);H(6) H(7)]');
      
       
             
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
        

 HQ=XX\rq';

        Q=[HQ(1),HQ(2);
           HQ(2),HQ(3)];

        eQ=[eQ,norm(Q)];
        jj=0;
        stopq=0;
        X=[];
        XX=[];
        Hqxx=[];
        Hqxu=[];
        Hqxw=[];
        Hquu=[];
        Hquw=[];
        Hqww=[];
        Hq=[];
        
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



