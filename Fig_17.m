% Fig. 17: Solutions of the competitive threshold-linear model in the
% perturbed and unperturbed cases

epsilon = 0.25; 
delta = 0.5;
theta1 = 1; theta2 = 1; theta3 = 1;
theta1_pert = 1.05;

W13 = -1+epsilon; W21 = -1+epsilon; W32 = -1+epsilon; 
W12 = -1-delta; W23 = -1-delta; W31 = -1-delta; 

T0=11.25; %unperturbed period
Tp=11.4;  %perturbed period
dt = 0.001;
tspan_u = 0:dt:T0/3+2*Tp;
tspan_p = T0/3:dt:T0/3+2*Tp;  %apply perturbation since time T0/3

initials = [0.0173559447648594,0.466544170059006,0.465889877817247];
initials_p = [0.466275884479895,0.0173266745359748,0.466137590369761];

[~,P] = ode45(@threshold_linear,tspan_u,initials,[],epsilon,delta,theta1,theta2,theta3);
x1 = P(:,1); x2 = P(:,2); x3 = P(:,3);

[~,Pp] = ode45(@threshold_linear,tspan_p,initials_p,[],epsilon,delta,theta1_pert,theta2,theta3);
x1p = Pp(:,1); x2p = Pp(:,2); x3p = Pp(:,3);

figure
plot(tspan_u,x1,'--k','LineWidth',1.5); hold on
plot(tspan_u,x2,'--b','LineWidth',1.5); 
plot(tspan_u,x3,'--r','LineWidth',1.5); 
plot(tspan_p,x1p,'-k','LineWidth',1.5); 
plot(tspan_p,x2p,'-b','LineWidth',1.5); 
plot(tspan_p,x3p,'-r','LineWidth',1.5); 
plot(T0/3,initials_p(1),'.k','MarkerSize',25); hold off
legend('unperturbed x_1','unperturbed x_2','unperturbed x_3','perturbed x_1','perturbed x_2','perturbed x_3')
xlim([0 max(tspan_p)])
xlabel('time'); ylabel('node activity')
set(gca','FontSize',13);