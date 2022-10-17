% Fig. 4: A solution of synaptic-escape system
% Change the following parameters in relaxation.m
% G = struct('NaP', 6.8, 'L', 3, 'I', 0.4, 'E', 0.1); % conductances, nS
% Theta = struct('h', -40, 'mp', -37); % half activations (mV)
% Sigma = struct('h', 5, 'mp', -6); % slopes

clear; clc;

d1 = 1; d2 = 1; d3 = 1;
theta_I = -62;  %synaptic escape
tF = 150;

initials = [-10.0000  -62.6063  -63.9030    0.4049    0.7455    0.3885];
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[T,P] = ode15s(@relaxation,[0 tF],initials,options,d1,d2,d3,theta_I);
v1 = P(:,1); v2 = P(:,2); v3 = P(:,3); h1 = P(:,4); h2 = P(:,5); h3 = P(:,6);

figure  %trajectories
subplot(3,1,1)
plot(T,v1,'-k','LineWidth',2); hold on
plot(T,theta_I*ones(length(T),1),'--m','LineWidth',2); hold off
ylabel('v_1 (mV)'); set(gca,'FontSize',13)
subplot(3,1,2)
plot(T,v2,'-b','LineWidth',2); hold on
plot(T,theta_I*ones(length(T),1),'--m','LineWidth',2); hold off
ylabel('v_2 (mV)'); set(gca,'FontSize',13)
subplot(3,1,3)
plot(T,v3,'-r','LineWidth',2); hold on
plot(T,theta_I*ones(length(T),1),'--m','LineWidth',2); hold off
ylabel('v_3 (mV)'); xlabel('time (ms)'); 
set(gca,'FontSize',13)

% Plot nullclines of cell 1
v = -70:0.1:16; 

E = struct('Na', 50, 'L', -65, 'I', -80, 'E', 0); % reversal potentials, mV
sigma_I = -0.01;
G = struct('NaP', 6.8, 'L', 3, 'I', 0.4, 'E', 0.1); % conductances, nS
Theta = struct('h', -40, 'mp', -37); % half activations (mV)
Sigma = struct('h', 5, 'mp', -6); % slopes
B = struct('b12', 1, 'b13', 1, 'b21', 1, 'b23', 1, 'b31', 1, 'b32', 1); % coupling constants

mp_inf = @(v) 1./(1+exp((v-Theta.mp)/Sigma.mp));  
h_null = 1./(1+exp((v-Theta.h)/Sigma.h)); 
vfree_null = (-G.L*(v-E.L)-G.E*d1*(v-E.E))./(G.NaP*mp_inf(v).*(v-E.Na));
vinhibited_null = (-G.L*(v-E.L)-G.I*B.b21*(v-E.I)-G.E*d1*(v-E.E))./(G.NaP*mp_inf(v).*(v-E.Na));

figure  %nullclines 
plot(v1,h1,'-k','LineWidth',2); hold on
plot(theta_I*ones(1,length([0:.01:1])),[0:.01:1],':m','LineWidth',2);
plot(v,h_null,'-g','LineWidth',2)
plot(v,vfree_null,'-b','LineWidth',2)
plot(v,vinhibited_null,'--b','LineWidth',2); hold off
axis([-68, 15, 0, 1])
xlabel('v_1 (mV)'); ylabel('h_1')
set(gca,'FontSize',13)