% Fig. 5: HC and LC solutions of piecewise linear Aplysia model

clear; clc;

% Heteroclinic cycle case
a1 = 0.0001; a2 = 0.0001; a3 = 0.0001;
rho = 3;

tF = 40; dt = 0.001; T = 0:dt:tF;
initials1 = [ 0.4462    0.0001    0.3403];

P1 = Aplysia(initials1,0,tF,dt,a1,a2,a3,rho);
X1 = P1(:,1); Y1 = P1(:,2); Z1 = P1(:,3);

% Limit cycle case
a1 = 0.01; a2 = 0.01; a3 = 0.01;
rho = 3;

initials2 = [0.500342977762597 0.0172224538514646 0.237135024241301];

P2 = Aplysia(initials2,0,tF,dt,a1,a2,a3,rho);
X2 = P2(:,1); Y2 = P2(:,2); Z2 = P2(:,3);

figure

subplot(1,2,1)  %hetroclinic cycle
%trajectory
plot3(X1,Y1,Z1,'-k','LineWidth',2); hold on
%fixed points
plot3(1,0,0,'.m','MarkerSize',28);
plot3(0,1,0,'.g','MarkerSize',28);
plot3(0,0,1,'.y','MarkerSize',28);
%switching boundaries
%region 1
plot3([0 1],[0,0],[0,0],'-m','LineWidth',2)
plot3([1 1],[0,0],[1,0],'-m','LineWidth',2)
plot3([1 0],[0,0],[1,0],'-m','LineWidth',2)
plot3([0 1],[0,1],[0,1],'-m','LineWidth',2)
plot3([1,1],[0,1],[0,0],'-m','LineWidth',2)
plot3([0,1],[0,1],[0,0],'-m','LineWidth',2)
plot3([1,1],[1,1],[0,1],'-m','LineWidth',2)
plot3([1 1],[1,0],[1,1],'-m','LineWidth',2)
%region 2
plot3([1,0],[1 1],[0,0],'-g','LineWidth',2)
plot3([1,0],[1,0],[0,0],'-g','LineWidth',2)
plot3([0,0],[0,1],[0,0],'-g','LineWidth',2)
plot3([1,1],[1,1],[0,0],'-g','LineWidth',2)
plot3([1,0],[1,1],[1,1],'-g','LineWidth',2)
plot3([0,0],[1,1],[0,1],'-g','LineWidth',2)
plot3([0,1],[0,1],[0,1],'-g','LineWidth',2)
%region 3
plot3([0,0],[0,0],[0,1],'-y','LineWidth',2)
plot3([0,0],[0,1],[1,1],'-y','LineWidth',2)
plot3([0,0],[0,1],[0,1],'-y','LineWidth',2)
plot3([0,1],[0,0],[1,1],'-y','LineWidth',2)
plot3([1,0],[0,0],[1,0],'-y','LineWidth',2)
plot3([1,1],[0,1],[1,1],'-y','LineWidth',2)
plot3([1,0],[1,0],[1,0],'-y','LineWidth',2)
plot3([1,0],[1,1],[1,1],'-y','LineWidth',2); hold off
set(gca,'XDir','reverse'); set(gca,'YDir','reverse'); 
xlabel('x'); ylabel('y'); zlabel('z'); 
title('(a) a_1 = 0')
grid on; axis([0 1 0 1 0 1])

subplot(1,2,2)  %limit cycle
%trajectory
plot3(X2,Y2,Z2,'-k','LineWidth',2); hold on
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
title('(b) a_1 = 0.01')
grid on; axis([-2*a2 1+2*a2 -2*a2 1+2*a2 -2*a2 1+2*a2])