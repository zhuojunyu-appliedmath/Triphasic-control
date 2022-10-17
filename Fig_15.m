% Fig. 15: Solutions of piecewise linear Aplysia model with several values of a1

clear; clc;

a2 = 0.01; a3 = 0.01;
rho = 3;
tF = 15; dt = 0.001; T = 0:dt:tF;

%% a1 = 0.01;
a1 = 0.01;
initials1 = [0.500342977762597 0.0172224538514646 0.237135024241301];
P1 = Aplysia(initials1,0,tF,dt,a1,a2,a3,rho);
X1 = P1(:,1); Y1 = P1(:,2); Z1 = P1(:,3);

%% a1 = 0.02;
a1 = 0.02;
initials2 = [0.500250761643061,0.0224291763985092,0.217964182087826];
P2 = Aplysia(initials2,0,tF,dt,a1,a2,a3,rho);
X2 = P2(:,1); Y2 = P2(:,2); Z2 = P2(:,3);

%% a1 = 0.0005;
a1 = 0.0005;
initials3 = [0.500240461019345,0.0150111754691418,0.253098779673146];
P3 = Aplysia(initials3,0,tF,dt,a1,a2,a3,rho);
X3 = P3(:,1); Y3 = P3(:,2); Z3 = P3(:,3);

%% Plot solutions

figure

subplot(3,1,1)
plot(T,X2,'-k','LineWidth',2); hold on
plot(T,Y2,'-b','LineWidth',2); 
plot(T,Z2,'-r','LineWidth',2); hold off
ylabel('x, y, z'); title('a1 = 0.02');
set(gca,'FontSize',12);
subplot(3,1,2)
plot(T,X1,'-k','LineWidth',2); hold on
plot(T,Y1,'-b','LineWidth',2); 
plot(T,Z1,'-r','LineWidth',2); hold off
ylabel('x, y, z'); title('a1 = 0.01');
set(gca,'FontSize',12);
subplot(3,1,3)
plot(T,X3,'-k','LineWidth',2); hold on
plot(T,Y3,'-b','LineWidth',2); 
plot(T,Z3,'-r','LineWidth',2); hold off
xlabel('time'); ylabel('x, y, z'); title('a1 = 0.0005');
set(gca,'FontSize',12);

figure
a1=0.01;
%trajectory
plot3(X1,Y1,Z1,'-k','LineWidth',2); hold on
plot3(X2,Y2,Z2,'-r','LineWidth',2);
plot3(X3,Y3,Z3,'-b','LineWidth',2);
%fixed points
plot3(1-(a1-a2)*rho,-a2,a3,'.m','MarkerSize',28);
plot3(a1,1-(a2-a3)*rho,-a3,'.g','MarkerSize',28);
plot3(-a1,a2,1-(a3-a1)*rho,'.y','MarkerSize',28);
%switching boundaries
%region 1
plot3([0 1-(a1-a2)*rho],[-(a1+a2)/2,-a2],[a3,a3],'-m','LineWidth',2)
plot3([1 1-(a1-a2)*rho],[-a2,-a2],[1,a3],'-m','LineWidth',2)
plot3([1 0],[-a2,-(a1+a2)/2],[1,a3],'-m','LineWidth',2)
plot3([0 1],[-(a1+a2)/2,1-(a1+a2)/2],[a3,1],'-m','LineWidth',2)
plot3([1-(a1-a2)*rho,1],[-a2,1-(a1+a2)/2],[a3,a3],'-m','LineWidth',2)
plot3([0,1],[-(a1+a2)/2,1-(a1+a2)/2],[a3,a3],'-m','LineWidth',2)
plot3([1,1],[1-(a1+a2)/2,1-(a1+a2)/2],[a3,1],'-m','LineWidth',2)
plot3([1 1],[1-(a1+a2)/2,-a2],[1,1],'-m','LineWidth',2)
%region 2
plot3([1,a1],[1 1-(a2-a3)*rho],[-a3,-a3],'-g','LineWidth',2)
plot3([1,a1],[1,0],[-a3,-(a2+a3)/2],'-g','LineWidth',2)
plot3([a1,a1],[0,1-(a2-a3)*rho],[-(a2+a3)/2,-a3],'-g','LineWidth',2)
plot3([1,1],[1,1],[-a3,1-(a2+a3)/2],'-g','LineWidth',2)
plot3([1,a1],[1,1],[1-(a2+a3)/2,1-(a2+a3)/2],'-g','LineWidth',2)
plot3([a1,a1],[1-(a2-a3)*rho,1],[-a3,1-(a2+a3)/2],'-g','LineWidth',2)
plot3([a1,1],[0,1],[-(a2+a3)/2,1-(a2+a3)/2],'-g','LineWidth',2)
%region 3
plot3([-(a1+a3)/2,-a1],[a2,a2],[0,1-(a3-a1)*rho],'-y','LineWidth',2)
plot3([-a1,-a1],[a2,1],[1-(a3-a1)*rho,1],'-y','LineWidth',2)
plot3([-(a2+a3)/2,-a1],[a2,1],[0,1],'-y','LineWidth',2)
plot3([-a1,1-(a1+a3)/2],[a2,a2],[1-(a3-a1)*rho,1],'-y','LineWidth',2)
plot3([1-(a1+a3)/2,-(a1+a3)/2],[a2,a2],[1,0],'-y','LineWidth',2)
plot3([1-(a1+a3)/2,1-(a1+a3)/2],[a2,1],[1,1],'-y','LineWidth',2)
plot3([1-(a1+a3)/2,-(a1+a3)/2],[1,a2],[1,0],'-y','LineWidth',2)
plot3([1-(a1+a3)/2,-a1],[1,1],[1,1],'-y','LineWidth',2); hold off
set(gca,'XDir','reverse'); set(gca,'YDir','reverse'); 
xlabel('x'); ylabel('y'); zlabel('z'); 
grid on; axis([-2*a2 1+2*a2 -2*a2 1+2*a2 -2*a2 1+2*a2])
set(gca,'FontSize',12);