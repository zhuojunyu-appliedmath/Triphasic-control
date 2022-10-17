% Fig. 12: Phase plane for three cells in synaptic escape, with several values of d1
% Change the following parameters in relaxation.m
% G = struct('NaP', 6.8, 'L', 3, 'I', 0.4, 'E', 0.1); % conductances, nS
% Theta = struct('h', -40, 'mp', -37); % half activations (mV)
% Sigma = struct('h', 5, 'mp', -6); % slopes

clear; clc;

d2 = 1; d3 = 1;
theta_I = -62;  %synaptic escape
tF = 80; dt = 0.01; tspan = 0:dt:tF;

%% d1 = 1;

d1 = 1; 
initials1 = [-10.0000  -62.6063  -63.9030    0.4049    0.7455    0.3885];
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[~,P1] = ode15s(@relaxation,tspan,initials1,options,d1,d2,d3,theta_I);
v1_1 = P1(:,1); v2_1 = P1(:,2); v3_1 = P1(:,3); h1_1 = P1(:,4); h2_1 = P1(:,5); h3_1 = P1(:,6);

E = struct('Na', 50, 'L', -65, 'I', -80, 'E', 0); % reversal potentials, mV
sigma_I = -0.01;
G = struct('NaP', 6.8, 'L', 3, 'I', 0.4, 'E', 0.1); % conductances, nS
Theta = struct('h', -40, 'mp', -37); % half activations (mV)
Sigma = struct('h', 6, 'mp', -6); % slopes
B = struct('b12', 1, 'b13', 1, 'b21', 1, 'b23', 1, 'b31', 1, 'b32', 1); % coupling constants
mp_inf = @(v) 1./(1+exp((v-Theta.mp)/Sigma.mp));  
v = -70:0.1:16; 

h_null = 1./(1+exp((v-Theta.h)/Sigma.h)); 

% V-nullclines of cell 1
vfree_null1_1 = (-G.L*(v-E.L)-G.E*d1*(v-E.E))./(G.NaP*mp_inf(v).*(v-E.Na));
vinhibited_null1_1 = (-G.L*(v-E.L)-G.I*B.b21*(v-E.I)-G.E*d1*(v-E.E))./(G.NaP*mp_inf(v).*(v-E.Na));
% V-nullclines of cell 2
vfree_null2 = (-G.L*(v-E.L)-G.E*d2*(v-E.E))./(G.NaP*mp_inf(v).*(v-E.Na));
vinhibited_null2 = (-G.L*(v-E.L)-G.I*B.b21*(v-E.I)-G.E*d2*(v-E.E))./(G.NaP*mp_inf(v).*(v-E.Na));
% V-nullclines of cell 3
vfree_null3 = (-G.L*(v-E.L)-G.E*d3*(v-E.E))./(G.NaP*mp_inf(v).*(v-E.Na));
vinhibited_null3 = (-G.L*(v-E.L)-G.I*B.b21*(v-E.I)-G.E*d3*(v-E.E))./(G.NaP*mp_inf(v).*(v-E.Na));


%% d1 = 1.03;

d1 = 1.03; 
initials2 = [ -62.4338  -63.6336  -16.2800    0.7641    0.4769    0.3300];
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[~,P2] = ode15s(@relaxation,tspan,initials2,options,d1,d2,d3,theta_I);
v1_2 = P2(:,1); v2_2 = P2(:,2); v3_2 = P2(:,3); h1_2 = P2(:,4); h2_2 = P2(:,5); h3_2 = P2(:,6);

% Nullclines
vfree_null1_2 = (-G.L*(v-E.L)-G.E*d1*(v-E.E))./(G.NaP*mp_inf(v).*(v-E.Na));
vinhibited_null1_2 = (-G.L*(v-E.L)-G.I*B.b21*(v-E.I)-G.E*d1*(v-E.E))./(G.NaP*mp_inf(v).*(v-E.Na));

%% d1 = 0.97;

d1 = 0.97; 
initials3 = [ -0.9329  -62.6510  -64.1107    0.5517    0.7359    0.3147];
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[~,P3] = ode15s(@relaxation,tspan,initials3,options,d1,d2,d3,theta_I);
v1_3 = P3(:,1); v2_3 = P3(:,2); v3_3 = P3(:,3); h1_3 = P3(:,4); h2_3 = P3(:,5); h3_3 = P3(:,6);

% Nullclines
vfree_null1_3 = (-G.L*(v-E.L)-G.E*d1*(v-E.E))./(G.NaP*mp_inf(v).*(v-E.Na));
vinhibited_null1_3 = (-G.L*(v-E.L)-G.I*B.b21*(v-E.I)-G.E*d1*(v-E.E))./(G.NaP*mp_inf(v).*(v-E.Na));


%% Plot
figure  

subplot(1,3,1)  %cell 1
plot(v1_1,h1_1,'.k','MarkerSize',3); hold on
plot(v1_3,h1_3,'.b','MarkerSize',3);
plot(v1_2,h1_2,'.r','MarkerSize',3);
plot(theta_I*ones(1,length([0:.01:1])),[0:.01:1],':m','LineWidth',0.5);
plot(v,h_null,'-g','LineWidth',0.5)
plot(v,vfree_null1_1,'-k','LineWidth',0.5)
plot(v,vfree_null1_3,'-b','LineWidth',0.5)
plot(v,vfree_null1_2,'-r','LineWidth',0.5)
plot(v,vinhibited_null1_1,'--k','LineWidth',0.5); 
plot(v,vinhibited_null1_3,'--b','LineWidth',0.5);
plot(v,vinhibited_null1_2,'--r','LineWidth',0.5); hold off
axis([-68, 15, 0, 1])
xlabel('v_1 (mV)'); ylabel('h_1')
set(gca,'FontSize',11)

subplot(1,3,2) %cell 2
plot(v2_1,h2_1,'.k','MarkerSize',3); hold on
plot(v2_3,h2_3,'.b','MarkerSize',3);
plot(v2_2,h2_2,'.r','MarkerSize',3);
plot(theta_I*ones(1,length([0:.01:1])),[0:.01:1],':m','LineWidth',0.5);
plot(v,h_null,'-g','LineWidth',0.5)
plot(v,vfree_null2,'-r','LineWidth',0.5)
plot(v,vinhibited_null2,'--r','LineWidth',0.5); hold off
axis([-68, 15, 0, 1])
xlabel('v_2 (mV)'); ylabel('h_2')
set(gca,'FontSize',11)

subplot(1,3,3) %cell 3
plot(v3_1,h3_1,'.k','MarkerSize',3); hold on
plot(v3_3,h3_3,'.b','MarkerSize',3);
plot(v3_2,h3_2,'.r','MarkerSize',3);
plot(theta_I*ones(1,length([0:.01:1])),[0:.01:1],':m','LineWidth',0.5);
plot(v,h_null,'-g','LineWidth',0.5)
plot(v,vfree_null3,'-r','LineWidth',0.5)
plot(v,vinhibited_null3,'--r','LineWidth',0.5); hold off
axis([-68, 15, 0, 1])
xlabel('v_3 (mV)'); ylabel('h_3')
set(gca,'FontSize',11)