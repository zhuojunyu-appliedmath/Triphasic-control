% Fig. 8: Phase plane for cell 1 in synaptic release, with several values of d1
% Change the following parameters in relaxation.m
% G = struct('NaP', 6.8, 'L', 3, 'I', 0.4, 'E', 0.1); % conductances, nS
% Theta = struct('h', -40, 'mp', -37); % half activations (mV)
% Sigma = struct('h', 6, 'mp', -6); % slopes

clear; clc;

d2 = 1; d3 = 1;
theta_I = -25;  %intrinsic release
tF = 100; dt = 0.01; tspan = 0:dt:tF;

%% d1 = 1;

d1 = 1; 
initials1 = [-10.0000  -62.7983  -63.8956    0.4055    0.7024    0.3903];
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[~,P1] = ode15s(@relaxation,tspan,initials1,options,d1,d2,d3,theta_I);
v1_1 = P1(:,1); h1_1 = P1(:,4); 

% Nullclines of cell 1
E = struct('Na', 50, 'L', -65, 'I', -80, 'E', 0); % reversal potentials, mV
sigma_I = -0.01;
G = struct('NaP', 6.8, 'L', 3, 'I', 0.4, 'E', 0.1); % conductances, nS
Theta = struct('h', -40, 'mp', -37); % half activations (mV)
Sigma = struct('h', 6, 'mp', -6); % slopes
B = struct('b12', 1, 'b13', 1, 'b21', 1, 'b23', 1, 'b31', 1, 'b32', 1); % coupling constants
mp_inf = @(v) 1./(1+exp((v-Theta.mp)/Sigma.mp));  
v = -70:0.1:16; 

h_null = 1./(1+exp((v-Theta.h)/Sigma.h)); 
vfree_null1 = (-G.L*(v-E.L)-G.E*d1*(v-E.E))./(G.NaP*mp_inf(v).*(v-E.Na));
vinhibited_null1 = (-G.L*(v-E.L)-G.I*B.b21*(v-E.I)-G.E*d1*(v-E.E))./(G.NaP*mp_inf(v).*(v-E.Na));

%% d1 = 1.3;

d1 = 1.3; 
initials2 = [-61.8148  -63.5400  -17.7075    0.7356    0.5051    0.3160];
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[~,P2] = ode15s(@relaxation,tspan,initials2,options,d1,d2,d3,theta_I);
v1_2 = P2(:,1); h1_2 = P2(:,4); 

% Nullclines
vfree_null2 = (-G.L*(v-E.L)-G.E*d1*(v-E.E))./(G.NaP*mp_inf(v).*(v-E.Na));
vinhibited_null2 = (-G.L*(v-E.L)-G.I*B.b21*(v-E.I)-G.E*d1*(v-E.E))./(G.NaP*mp_inf(v).*(v-E.Na));

%% d1 = 0.7;

d1 = 0.7; 
initials3 = [-63.4510  -63.7734  -13.0830    0.7363    0.4314    0.3660];
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[~,P3] = ode15s(@relaxation,tspan,initials3,options,d1,d2,d3,theta_I);
v1_3 = P3(:,1); h1_3 = P3(:,4); 

% Nullclines
vfree_null3 = (-G.L*(v-E.L)-G.E*d1*(v-E.E))./(G.NaP*mp_inf(v).*(v-E.Na));
vinhibited_null3 = (-G.L*(v-E.L)-G.I*B.b21*(v-E.I)-G.E*d1*(v-E.E))./(G.NaP*mp_inf(v).*(v-E.Na));


%% Plot
figure  
plot(v1_1,h1_1,'.k','MarkerSize',6); hold on
plot(v1_2,h1_2,'.r','MarkerSize',6);
plot(v1_3,h1_3,'.b','MarkerSize',6);
plot(theta_I*ones(1,length([0:.01:1])),[0:.01:1],':m','LineWidth',2);
plot(v,h_null,'-g','LineWidth',1.5)
plot(v,vfree_null1,'-k','LineWidth',1)
plot(v,vfree_null2,'-r','LineWidth',1)
plot(v,vfree_null3,'-b','LineWidth',1)
plot(v,vinhibited_null1,'--k','LineWidth',1); 
plot(v,vinhibited_null2,'--r','LineWidth',1);
plot(v,vinhibited_null3,'--b','LineWidth',1);hold off
axis([-68, 15, 0, 1])
xlabel('v_1 (mV)'); ylabel('h_1')
set(gca,'FontSize',13)