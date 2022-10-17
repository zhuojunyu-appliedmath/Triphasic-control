% Fig. 6: A solution of competitive threshold-linear model

epsilon = 0.25; 
delta = 0.5;
theta1 = 1; theta2 = 1; theta3 = 1;

W13 = -1+epsilon; W21 = -1+epsilon; W32 = -1+epsilon; 
W12 = -1-delta; W23 = -1-delta; W31 = -1-delta; 

tF = 20; dt = 0.001; tspan = 0:dt:tF;
initials = [0.501193190683074,0.0139591696003762,0.425786237781174];
[T,P] = ode45(@threshold_linear,tspan,initials,[],epsilon,delta,theta1,theta2,theta3);
x1 = P(:,1); x2 = P(:,2); x3 = P(:,3);

theta1_pert = 1.05;
[Tp,Pp] = ode45(@threshold_linear,tspan,initials,[],epsilon,delta,theta1_pert,theta2,theta3);
x1p = Pp(:,1); x2p = Pp(:,2); x3p = Pp(:,3);

% Mark nonsmooth points
term1 = max(W13*x3+W12*x2+theta1,0); 
term2 = max(W21*x1+W23*x3+theta2,0);
term3 = max(W32*x2+W31*x1+theta3,0);

zeropt1 = find(term1==0); ind1 = diff(zeropt1); tempt1 = find(ind1~=1);
zeropt2 = find(term2==0); ind2 = diff(zeropt2); tempt2 = find(ind2~=1);
zeropt3 = find(term3==0); ind3 = diff(zeropt3); tempt3 = find(ind3~=1);

nonsmoothpt1 = [zeropt1(1) zeropt1(tempt1) zeropt1(tempt1+1) zeropt1(end)];
nonsmoothpt2 = [zeropt2(tempt2(1)) zeropt2(tempt2(1)+1) zeropt2(tempt2(2)) zeropt2(tempt2(2)+1)];
nonsmoothpt3 = [zeropt3(1) zeropt3(tempt3) zeropt3(tempt3+1) zeropt3(end)];

figure(1)
p1 = plot(T,x1,'-k','LineWidth',2); hold on
p2 = plot(T,x2,'-b','LineWidth',2);
p3 = plot(T,x3,'-r','LineWidth',2); 
plot(T(nonsmoothpt1),x1(nonsmoothpt1),'*k','MarkerSize',8)
plot(T(nonsmoothpt2),x2(nonsmoothpt2),'*b','MarkerSize',8)
plot(T(nonsmoothpt3),x3(nonsmoothpt3),'*r','MarkerSize',8); hold off
legend([p1 p2 p3],'x1','x2','x3')
xlabel('time'); ylabel('node activity')
set(gca,'FontSize',13);

figure(2)
plot3(x1,x2,x3','-k','LineWidth',2); hold on
plot3(x1(nonsmoothpt1),x2(nonsmoothpt1),x3(nonsmoothpt1),'*k','MarkerSize',8)
plot3(x1(nonsmoothpt2),x2(nonsmoothpt2),x3(nonsmoothpt2),'*b','MarkerSize',8)
plot3(x1(nonsmoothpt3),x2(nonsmoothpt3),x3(nonsmoothpt3),'*r','MarkerSize',8);
%region 1
plot3([1 1],[0 1],[0 0],'-m','LineWidth',2)
plot3([1 1],[0 0],[0 1],'-m','LineWidth',2)
plot3([1 0],[0 0],[0 0],'-m','LineWidth',2)
plot3([1 0],[1 0],[0 0],'-m','LineWidth',2)
plot3([1,1],[1 1],[0 1],'-m','LineWidth',2)
plot3([1,1],[0 1],[1 1],'-m','LineWidth',2)
plot3([1 0],[0,0],[1,0],'-m','LineWidth',2)
plot3([1 0],[1,0],[1,0],'-m','LineWidth',2)
%region 3
plot3([0 0],[0 0],[0,1],'-y','LineWidth',2)
plot3([0 0],[0 1],[1,1],'-y','LineWidth',2)
plot3([0 0],[0 1],[0,1],'-y','LineWidth',2)
plot3([0 1],[0 0],[1,1],'-y','LineWidth',2)
plot3([1 0],[0 0],[1,0],'-y','LineWidth',2)
plot3([1 1],[0 1],[1,1],'-y','LineWidth',2)
plot3([1 0],[1 0],[1,0],'-y','LineWidth',2)
plot3([1 0],[1,1],[1,1],'-y','LineWidth',2); 
%region 2
plot3([1 0],[1 1],[0 0],'-g','LineWidth',2)
plot3([1,0],[1,0],[0 0],'-g','LineWidth',2)
plot3([0 0],[0,1],[0 0],'-g','LineWidth',2)
plot3([1,1],[1,1],[0 1],'-g','LineWidth',2)
plot3([1 0],[1,1],[1 1],'-g','LineWidth',2)
plot3([0 0],[1,1],[0 1],'-g','LineWidth',2)
plot3([0 0],[0,1],[0 1],'-g','LineWidth',2)
plot3([0 1],[0,1],[0 1],'-g','LineWidth',2); hold off
xlabel('x_1'); ylabel('x_2'); zlabel('x_3')
grid on
set(gca,'XDir','reverse'); set(gca,'YDir','reverse'); 
set(gca,'FontSize',13);