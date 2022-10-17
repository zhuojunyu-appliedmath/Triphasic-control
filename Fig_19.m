% Fig. 19: Compare the lTRC and iPRC of relaxation oscillation model
% Change the following parameters in relaxation.m
% G = struct('NaP', 6.8, 'L', 3, 'I', 0.4, 'E', 0.1); % conductances, nS
% Theta = struct('h', -40, 'mp', -37); % half activations (mV)
% Sigma = struct('h', 6, 'mp', -6); % slopes

clear; clc;

%% Parameters

d1 = 1; d2 = 1; d3 = 1;

theta_I = -43; %intrinsic release
%theta_I = -25; %synaptic release

C = 0.21; epsilon = 0.01; sigma_I = -0.01;
E = struct('Na', 50, 'L', -65, 'I', -80, 'E', 0); %reversal potentials, mV
G = struct('NaP', 6.8, 'L', 3, 'I', 0.4, 'E', 0.1); %conductances, nS
Theta = struct('h', -40, 'mp', -37); %half activations (mV)
Sigma = struct('h', 6, 'mp', -6); %slopes
B = struct('b12', 1, 'b13', 1, 'b21', 1, 'b23', 1, 'b31', 1, 'b32', 1); %coupling constants

%% Vector field and Jacobian

S_inf = @(v) 1/(1+exp((v-theta_I)/sigma_I));
h_inf = @(v) 1/(1+exp((v-Theta.h)/Sigma.h));
mp_inf = @(v) 1/(1+exp((v-Theta.mp)/Sigma.mp));
tauh = @(v) epsilon*cosh((v-Theta.h)/Sigma.h/2);

dS_inf = @(v) -1/sigma_I/(2+exp((v-theta_I)/sigma_I)+1/exp((v-theta_I)/sigma_I)); 
dh_inf = @(v) -1/Sigma.h/(2+exp((v-Theta.h)/Sigma.h)+1/exp((v-Theta.h)/Sigma.h)); 
dmp_inf = @(v) -1/Sigma.mp/(2+exp((v-Theta.mp)/Sigma.mp)+1/exp((v-Theta.mp)/Sigma.mp)); 
dtauh = @(v) epsilon*sinh((v-Theta.h)/Sigma.h/2)/Sigma.h/2;

I_NaP = @(v,h) G.NaP*mp_inf(v)*h*(v-E.Na);
I_L = @(v) G.L*(v-E.L);

dI_NaPdv = @(v,h) G.NaP*dmp_inf(v)*h*(v-E.Na)+G.NaP*mp_inf(v)*h;
dI_NaPdh = @(v) G.NaP*mp_inf(v)*(v-E.Na);

% unperturbed
F = @(v1,v2,v3,h1,h2,h3) [(-I_NaP(v1,h1)-I_L(v1)-G.I*(B.b21*S_inf(v2)+B.b31*S_inf(v3))*(v1-E.I)-G.E*d1*(v1-E.E))/C;
                          (-I_NaP(v2,h2)-I_L(v2)-G.I*(B.b12*S_inf(v1)+B.b32*S_inf(v3))*(v2-E.I)-G.E*d2*(v2-E.E))/C;
                          (-I_NaP(v3,h3)-I_L(v3)-G.I*(B.b13*S_inf(v1)+B.b23*S_inf(v2))*(v3-E.I)-G.E*d3*(v3-E.E))/C;
                          (h_inf(v1)-h1)*tauh(v1);
                          (h_inf(v2)-h2)*tauh(v2);
                          (h_inf(v3)-h3)*tauh(v3)];
A = @(v1,v2,v3,h1,h2,h3) [(-dI_NaPdv(v1,h1)-G.L-G.I*(B.b21*S_inf(v2)+B.b31*S_inf(v3))-G.E*d1)/C  -G.I*B.b21*dS_inf(v2)*(v1-E.I)/C  -G.I*B.b31*dS_inf(v3)*(v1-E.I)/C  -dI_NaPdh(v1)/C  0  0;
                          -G.I*B.b12*dS_inf(v1)*(v2-E.I)/C  (-dI_NaPdv(v2,h2)-G.L-G.I*(B.b12*S_inf(v1)+B.b32*S_inf(v3))-G.E*d2)/C  -G.I*B.b32*dS_inf(v3)*(v2-E.I)/C  0  -dI_NaPdh(v2)/C  0;
                          -G.I*B.b13*dS_inf(v1)*(v3-E.I)/C  -G.I*B.b23*dS_inf(v2)*(v3-E.I)/C  (-dI_NaPdv(v3,h3)-G.L-G.I*(B.b13*S_inf(v1)+B.b23*S_inf(v2))-G.E*d3)/C  0  0  -dI_NaPdh(v3)/C;
                          dh_inf(v1)*tauh(v1)+(h_inf(v1)-h1)*dtauh(v1)  0  0  -tauh(v1)  0  0;
                          0  dh_inf(v2)*tauh(v2)+(h_inf(v2)-h2)*dtauh(v2)  0  0  -tauh(v2)  0;
                          0  0  dh_inf(v3)*tauh(v3)+(h_inf(v3)-h3)*dtauh(v3)  0  0  -tauh(v3)];
                      
%% Limit cycle               

initials = [-43.0000  -58.1205  -61.8513    0.8814    0.7286    0.2525]; %intrinsic release
T0 = 87.97; 
%initials = [ -25.0000  -59.2005  -25.5013    0.8129    0.6439    0.2588]; %synaptic release
%T0 = 61.97;
dt = 0.01; tspan = 0:dt:8*T0; 
Num = floor(T0/3/dt);
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[~,P1] = ode15s(@relaxation,tspan,initials,options,d1,d2,d3,theta_I);
v1 = P1(:,1); v2 = P1(:,2); v3 = P1(:,3); h1 = P1(:,4); h2 = P1(:,5); h3 = P1(:,6);

FF = zeros(6,length(tspan)); DF=zeros(6,6,length(tspan));
for i = 1:length(tspan)
    FF(:,i) = F(v1(i),v2(i),v3(i),h1(i),h2(i),h3(i));
    DF(:,:,i) = A(v1(i),v2(i),v3(i),h1(i),h2(i),h3(i));
end

%% iPRC
z = zeros(6,length(tspan));
z0 = [-0.999019521011081;0.018613927685713;0.035596773660868;-0.013605971929869;0.003787541087844;0.012121060192455]; %intrinsic release
%z0 = [0.998919836456036;-0.022684105069939;-0.035793206155744;0.013612810528020;-0.005241198446635;-0.012274336368471]; %synaptic release
k0 = FF(:,end)'*z0;
z(:,length(tspan)) = z0/k0;
for i = length(tspan)-1:-1:1
    z(:,i) = z(:,i+1)+DF(:,:,i+1)'*z(:,i+1)*dt;
    k = FF(:,i)'*z(:,i);
    z(:,i) = z(:,i)/k;
end

%% lTRC
n1_out = [-1; 0; 0; 0; 0; 0]; %normal vector of surface R1_out = {v1=theta_I,dv1dt<0}
eta1 = zeros(6,Num); 
FF_sub = FF(:,1:Num); DF_sub = DF(:,:,1:Num);
eta1(:,end) = -n1_out/(n1_out'*FF_sub(:,end));
for i = Num-1:-1:1
    eta1(:,i) = eta1(:,i+1)+DF_sub(:,:,i+1)'*eta1(:,i+1)*dt;
end

%% Plot
figure
subplot(3,1,1)  %v1
plot(tspan(1:Num),v1(1:Num),'-k','LineWidth',2); hold on
plot(tspan(1:Num),theta_I*ones(Num,1),'--m','LineWidth',2); hold off
ylabel('v_1 (mV)'); grid on;
xlim([0 T0/3]);
set(gca,'FontSize',12);
subplot(3,1,2)  %lTRC
plot(tspan(1:Num),eta1(1,:),'-k','LineWidth',2);
ylabel('lTRC, \eta_{v_1}'); grid on;
xlim([0 T0/3]);
set(gca,'FontSize',12);
subplot(3,1,3)  %iPRC 
plot(tspan(1:Num),z(1,1:Num),'-k','LineWidth',2);
ylabel('iPRC, z_{v_1}'); xlabel('time');
xlim([0 T0/3]);
grid on;
set(gca,'FontSize',12);