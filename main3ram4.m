set(0,'DefaultAxesLinewidth',2,'DefaultLineLineWidth',2);
set(0,'defaultAxesFontSize',14);
set(0,'defaultAxesFontName','arial');
set(0,'defaultTextFontName','arial');

close all
clear

A = [0.7 0.5 -0.1;0 0.7 0.1;-0.3 0 0.9];
Bd = [0.1;0.1;0.2];
Bu = [-1.2;-0.8;1.4];
C = [1 2 -1;0 -5 -0.2];
D = [0.01;0.02];

E = [1 1 1];

A0 = A;
A1 = A^2;
A2 = A^3;
A3 = A^4;
%A4 = A^5;
B0 = sqrt(3)*[zeros(3,1) Bd];
D0 = sqrt(3)*[D zeros(2,1)];
B1 = sqrt(4)*[zeros(3,1) Bd A*Bd];
D1 = sqrt(4)*[D zeros(2,1) zeros(2,1)];
B2 = sqrt(5)*[zeros(3,1) Bd A*Bd A^2*Bd];
D2 = sqrt(5)*[D zeros(2,1) zeros(2,1) zeros(2,1)];
B3 = sqrt(6)*[zeros(3,1) Bd A*Bd A^2*Bd A^3*Bd];
D3 = sqrt(6)*[D zeros(2,1) zeros(2,1) zeros(2,1) zeros(2,1)];


gammaoptimal = 100000;

for i = 1:12
    for j = 1:12
        for k = 1:15
            for l = 1:20

alpha0 = 0.3+0.03*i;
alpha1 = 0.3+0.03*j;
alpha2 = 0.4+0.03*k;
alpha3 = 0.05*l;
setlmis([])

[gamma,n,sgamma] = lmivar(1,[1 1]);
[P,n,sP] = lmivar(1,[3 1]);
[Y0,n,sY0] = lmivar(2,[3 2]);
[Y1,n,sY1] = lmivar(2,[3 2]);
[Y2,n,sY2] = lmivar(2,[3 2]);
[Y3,n,sY3] = lmivar(2,[3 2]);


S0 = newlmi;
lmiterm([-S0 1 1 P],(1-alpha0),1)
lmiterm([-S0 3 1 P],1,A0)
lmiterm([-S0 3 1 Y0],eye(3),C)
lmiterm([-S0 3 2 P],1,B0)
lmiterm([-S0 2 2 0],alpha0)
lmiterm([-S0 3 2 Y0],eye(3),D0)
lmiterm([-S0 3 3 P],1,1)

S1 = newlmi;
lmiterm([-S1 1 1 P],(1-alpha1),1)
lmiterm([-S1 3 1 P],1,A1)
lmiterm([-S1 3 1 Y1],1,C)
lmiterm([-S1 3 2 P],1,B1)
lmiterm([-S1 2 2 0],alpha1)
lmiterm([-S1 3 2 Y1],1,D1)
lmiterm([-S1 3 3 P],1,1)

S2 = newlmi;
lmiterm([-S2 1 1 P],(1-alpha2),1)
lmiterm([-S2 3 1 P],1,A2)
lmiterm([-S2 3 1 Y2],1,C)
lmiterm([-S2 3 2 P],1,B2)
lmiterm([-S2 2 2 0],alpha2)
lmiterm([-S2 3 2 Y2],1,D2)
lmiterm([-S2 3 3 P],1,1)


S3 = newlmi;
lmiterm([-S3 1 1 P],(1-alpha3),1)
lmiterm([-S3 3 1 P],1,A3)
lmiterm([-S3 3 1 Y3],1,C)
lmiterm([-S3 3 2 P],1,B3)
lmiterm([-S3 2 2 0],alpha3)
lmiterm([-S3 3 2 Y3],1,D3)
lmiterm([-S3 3 3 P],1,1)

Sn = newlmi;
lmiterm([-Sn 1 1 P],1,1)
lmiterm([-Sn 2 1 0],E)
lmiterm([-Sn 2 2 gamma],1,1)

LMIs = getlmis;
c = [1 zeros(1,30)];

[copt,xopt] = mincx(LMIs,c);

if isempty(xopt)
else
Popt = dec2mat(LMIs,xopt,P);
Y0opt = dec2mat(LMIs,xopt,Y0);
Y1opt = dec2mat(LMIs,xopt,Y1);
Y2opt = dec2mat(LMIs,xopt,Y2);
Y3opt = dec2mat(LMIs,xopt,Y3);
gammaopt2 = dec2mat(LMIs,xopt,gamma);
L0p = -Popt^(-1)*Y0opt;
L1p = -Popt^(-1)*Y1opt;
L2p = -Popt^(-1)*Y2opt;
L3p = -Popt^(-1)*Y3opt;

al0=abs(max(eig(A0-L0p*C)));
al1=abs(max(eig(A1-L1p*C)));
al2=abs(max(eig(A2-L2p*C)));
al3=abs(max(eig(A3-L3p*C)));

if (gammaopt2 < gammaoptimal)&&(alpha0<1-al0^2)&&(alpha1<1-al1^2)&&(alpha2<1-al2^2)&&(alpha3<1-al3^2)
    gammaoptimal = gammaopt2;
    L0 = L0p;
    L1 = L1p;
    L2 = L2p;
    L3 = L3p;
    iopti = i;
    jopti = j;
    kopti = k;
    lopti = l;
end
end
dellmi(LMIs,S0);
dellmi(LMIs,S1);
dellmi(LMIs,S2);
dellmi(LMIs,S3);
dellmi(LMIs,Sn);
            end
        end
    end
end


x(1,1) = 0;
x(2,1) = 0;
x(3,1) = 0;
xa(1,1) = 0;
xa(2,1) = 0;
xa(3,1) = 0;
x0(1,1) = 0;
x0(2,1) = 0;
x0(3,1) = 0;
x1(1,1) = 0;
x1(2,1) = 0;
x1(3,1) = 0;
x2(1,1) = 0;
x2(2,1) = 0;
x2(3,1) = 0;
x3(1,1) = 0;
x3(2,1) = 0;
x3(3,1) = 0;
xtotal(1,1) = 0;
xtotal(2,1) = 0;
xtotal(3,1) = 0;
Y(1,1) = x1(1,1);
Y(2,1) = x1(2,1);

Ts = 0.01;
t = 0:1:600;
a = 0:0.01:2*pi;
u = sin(a);
F = [];

e_proposed = 0;

co1 = 0;
co2 = 0;
co3 = 0;
co4 = 0;

for k = 1:600;
    x(:,k+1) = A * x(:,k) + Bu * u(:,k);
    y(:,k) = C * x(:,k);
    xa(:,k+1) = A * xa(:,k) + Bu * u(:,k) + Bd*(rand-0.5)*2;
    ya(:,k) = C * xa(:,k) + D*(rand-0.5)*2;

    if k < 5

    x0(:,k+1) = (A-L0*C)* xtotal(:,k) + Bu * u(:,k) + L0 * ya(:,k);
    x1(:,k+1) = (A-L0*C)* xtotal(:,k) + Bu * u(:,k) + L0 * ya(:,k);
    x2(:,k+1) = (A-L0*C)* xtotal(:,k) + Bu * u(:,k) + L0 * ya(:,k);
    x3(:,k+1) = (A-L0*C)* xtotal(:,k) + Bu * u(:,k) + L0 * ya(:,k);
            F(1,1) = norm(x0(:,k+1));
            F(2,1) = norm(x1(:,k+1));
            F(3,1) = norm(x2(:,k+1));
            F(4,1) = norm(x3(:,k+1));
            F(5,1) = norm(x1(:,k+1));
            F(6,1) = norm(x2(:,k+1));
            F(7,1) = norm(x3(:,k+1));
            F(8,1) = norm(x1(:,k+1));
            F(9,1) = norm(x2(:,k+1));

            [I1,I2] = sort(F);



            if I2(2) == 1
                xtotal(:,k+1) = x0(:,k+1);
            elseif I2(2) == 2
                xtotal(:,k+1) = x1(:,k+1);
            elseif I2(2) == 3
                xtotal(:,k+1) = x2(:,k+1);
            else
                xtotal(:,k+1) = x3(:,k+1);
            end
                     F2(k) = 0;
    else
            if rem(k,400) == 0
                ya(1,k) = 100;
                figure(2)
                plot(k,100,'gs',...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5]);
hold on;
            elseif rem(k,200) == 0
                ya(1,k) = -100;
                figure(2)
                plot(k,-100,'gs',...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5]);
hold on;
            elseif rem(k,100) == 0
                ya(1,k) = 0;
                figure(2)
                plot(k,0,'gs',...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5]);
hold on;
            end
            
            if rem(k,60) == 10
                ya(2,k) = 80;
                                figure(2)
                plot(k,80,'gs',...
    'MarkerSize',10,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor',[0.5,0.5,0.5]);
hold on;
            elseif rem(k,30) == 5
                ya(2,k) = -80;
                                figure(2)
                plot(k,-80,'gs',...
    'MarkerSize',10,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor',[0.5,0.5,0.5]);
hold on;
            end

            x0(:,k+1) = (A-L0*C)* xtotal(:,k) + Bu * u(:,k) + L0 * ya(:,k);
            x1(:,k+1) = (A^2-L1*C)*xtotal(:,k-1) + L1*ya(:,k-1) + A*Bu*u(:,k-1) + Bu*u(:,k);
            x2(:,k+1) = (A^3-L2*C) * xtotal(:,k-2) + L2*ya(:,k-2) + A^2*Bu*u(:,k-2) + A*Bu*u(:,k-1) + Bu*u(:,k);
            x3(:,k+1) = (A^4-L3*C) * xtotal(:,k-3) + L3*ya(:,k-3) + A^3*Bu*u(:,k-3) + A^2*Bu*u(:,k-2) + A*Bu*u(:,k-1) + Bu*u(:,k);

            F(1,1) =  [1 1] * ya(:,k) - [1 1] * C * xtotal(:,k-1) + [1 1] * C * Bu * u(:,k-1);
            F(2,1) =  [1 1] * ya(:,k) - [1 1] * C * xtotal(:,k-1) + [1 1] * C * Bu * u(:,k-1);
            F(3,1) = [1 1] * ya(:,k-1) - [1 1] * C * xtotal(:,k-2) + [1 1] * C * Bu * u(:,k-2);
            F(4,1) = [1 1] * ya(:,k-1) - [1 1] * C * xtotal(:,k-2) + [1 1] * C * Bu * u(:,k-2);
            F(5,1) = [1 1] * ya(:,k-2) - [1 1] * C * xtotal(:,k-3) + [1 1] * C * Bu * u(:,k-3);
            F(6,1) = [1 1] * ya(:,k-2) - [1 1] * C * xtotal(:,k-3) + [1 1] * C * Bu * u(:,k-3);
            F(7,1) = [1 1] * ya(:,k-3) - [1 1] * C * xtotal(:,k-4) + [1 1] * C * Bu * u(:,k-4);
            F(8,1) = [1 1] * ya(:,k-3) - [1 1] * C * xtotal(:,k-4) + [1 1] * C * Bu * u(:,k-4);
            [I1,I2] = sort(F);


            if I2(4) <= 2
                xtotal(:,k+1) = x0(:,k+1);
                         F2(k) = 1;
                         co1 = co1+1;
            elseif I2(4) <= 4
                xtotal(:,k+1) = x1(:,k+1);
                          F2(k) = 2;
                          co2 = co2+1;
            elseif I2(4) <= 6
                xtotal(:,k+1) = x2(:,k+1);
                          F2(k) = 3;
                          co3 = co3+1;
            elseif I2(4) <= 8
                xtotal(:,k+1) = x3(:,k+1);
                          F2(k) = 4;
                          co4 = co4+1;
            end

  % xtotal(:,k+1) = (A-L0*C)* xtotal(:,k) + B * u(:,k) + L0 * ya(:,k);

    end
         e_proposed = e_proposed + norm(xtotal(:,k)-xa(:,k));

end

figure(1),stairs(t,xa(1,:),'k','LineStyle','--'),hold on,stairs(t,xa(2,:),'k','LineStyle','--'),hold on,stairs(t,xa(3,:),'k','LineStyle','--')
figure(1),stairs(t,xtotal(1,:),'k'),hold on,stairs(t,xtotal(2,:),'k'),hold on,stairs(t,xtotal(3,:),'k'),axis([0 600 -30 30]);
%figure(2),stairs(t,x2(1,:),'b'), hold on
%figure(3),stairs(t,x3(1,:),'k'), hold on
%figure(4),stairs(t,xtotal(1,:),'k'), hold on

ylabel('x_p,x')
xlabel('k')

ya(1,601) = 0;
ya(2,601) = 0;


figure(2)
stairs(t,ya(1,:),'k')
hold on,axis([0 600 -105 105]);
stairs(t,ya(2,:),'k'),axis([0 600 -105 105]);
hold on

ylabel('y')
xlabel('k')

figure(3) %z,zQ1
stairs(t,E(1)*(xa(1,:)-xtotal(1,:))+E(2)*(xa(2,:)-xtotal(2,:))+E(3)*(xa(3,:)-xtotal(3,:)),'k','LineWidth',0.5);
axis([0 600 -2.5 2.5])
%set(gca,'fontsize',14);
xlabel('k')
ylabel('z_e')


figure(4) %z,zQ1
subplot(3,1,1)
stairs(t,(xa(1,:)-xtotal(1,:)),'k','LineWidth',0.5);
axis([0 600 -.5 .5])
ylabel('e_{1}')
xlabel('k')
%set(gca,'fontsize',14);
subplot(3,1,2)
stairs(t,xa(2,:)-xtotal(2,:),'k','LineWidth',0.5);
axis([0 600 -.5 .5])
%set(gca,'fontsize',14);
ylabel('e_{2}')
xlabel('k')
subplot(3,1,3)
stairs(t,xa(3,:)-xtotal(3,:),'k','LineWidth',0.5);
axis([0 600 -.5 .5])
%set(gca,'fontsize',14);
ylabel('e_{3}')
xlabel('k')

figure(5),plot(F2,'ok'),axis([0 600 0.5 4.5])

ylabel('selected candidate i')
xlabel('k')

gammaopt = sqrt(gammaoptimal)


    iopti
    jopti
    kopti
    L0
    L1
    L2
    L3