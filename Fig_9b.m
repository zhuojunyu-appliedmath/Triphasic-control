% Fig. 9(b): lTRC for synaptic release
% Change the following parameters in relaxation.m
% G = struct('NaP', 6.8, 'L', 3, 'I', 0.4, 'E', 0.1); % conductances, nS
% Theta = struct('h', -40, 'mp', -37); % half activations (mV)
% Sigma = struct('h', 6, 'mp', -6); % slopes

clear; clc;

%% Parameters

d1 = 1; d2 = 1; d3 = 1;
delta = 0.05;  %perturbation

theta_I = -25;  %synaptic release

C = 0.21; epsilon = 0.01; sigma_I = -0.01;
E = struct('Na', 50, 'L', -65, 'I', -80, 'E', 0); %reversal potentials, mV
G = struct('NaP', 6.8, 'L', 3, 'I', 0.4, 'E', 0.1); %conductances, nS
Theta = struct('h', -40, 'mp', -37); %half activations (mV)
Sigma = struct('h', 6, 'mp', -6); %slopes
B = struct('b12', 1, 'b13', 1, 'b21', 1, 'b23', 1, 'b31', 1, 'b32', 1); %coupling constants


%% Limit cycle

tF = 200; dt = 0.01; tspan = 0:dt:tF;
options1 = odeset('RelTol',1e-8,'AbsTol',1e-8);
% unperturbed
initials = [-10.0000  -62.7983  -63.8956    0.4055    0.7024    0.3903];
[TTu,Pu] = ode15s(@relaxation,tspan,initials,options1,d1,d2,d3,theta_I);
v1 = Pu(:,1); v2 = Pu(:,2); v3 = Pu(:,3); h1 = Pu(:,4); h2 = Pu(:,5); h3 = Pu(:,6);
% perturbed
initialsp = [-10.0000  -62.8000  -63.8958    0.4054    0.7020    0.3902];


%% Switching surfaces and normal direction

n1_out = [-1; 0; 0; 0; 0; 0]; %R1_out = {v1=theta_I,dv1dt<0}
n2_out = [0; -1; 0; 0; 0; 0]; %R2_out = {v2=theta_I,dv2dt<0}
n3_out = [0; 0; -1; 0; 0; 0]; %R3_out = {v3=theta_I,dv3dt<0}


%% Find enter/exit time 

% unperturbed
options2 = odeset('Events',@event2_up,'RelTol',1e-8,'AbsTol',1e-8);
[~,~,t_in2,x_in2,~] = ode15s(@relaxation,tspan,initials,options2,d1,d2,d3,theta_I);
options2 = odeset('Events',@event2_down,'RelTol',1e-8,'AbsTol',1e-8);
[~,~,t_out2,x_out2,~] = ode15s(@relaxation,tspan,initials,options2,d1,d2,d3,theta_I);

options2 = odeset('Events',@event3_up,'RelTol',1e-8,'AbsTol',1e-8);
[~,~,t_in3,x_in3,~] = ode15s(@relaxation,tspan,initials,options2,d1,d2,d3,theta_I);
options2 = odeset('Events',@event3_down,'RelTol',1e-8,'AbsTol',1e-8);
[~,~,t_out3,x_out3,~] = ode15s(@relaxation,tspan,initials,options2,d1,d2,d3,theta_I);

options2 = odeset('Events',@event1_up,'RelTol',1e-8,'AbsTol',1e-8);
[~,~,t_in1,x_in1,~] = ode15s(@relaxation,tspan,initials,options2,d1,d2,d3,theta_I);
options2 = odeset('Events',@event1_down,'RelTol',1e-8,'AbsTol',1e-8);
[~,~,t_out1,x_out1,~] = ode15s(@relaxation,tspan,initials,options2,d1,d2,d3,theta_I);
t_out1 = t_out1(2); x_out1 = x_out1(2,:);

% perturbed
options2 = odeset('Events',@event2_up,'RelTol',1e-8,'AbsTol',1e-8);
[~,~,tp_in2,xp_in2,~] = ode15s(@relaxation,tspan,initialsp,options2,d1+delta,d2,d3,theta_I);
options2 = odeset('Events',@event2_down,'RelTol',1e-8,'AbsTol',1e-8);
[~,~,tp_out2,xp_out2,~] = ode15s(@relaxation,tspan,initialsp,options2,d1+delta,d2,d3,theta_I);

options2 = odeset('Events',@event3_up,'RelTol',1e-8,'AbsTol',1e-8);
[~,~,tp_in3,xp_in3,~] = ode15s(@relaxation,tspan,initialsp,options2,d1+delta,d2,d3,theta_I);
options2 = odeset('Events',@event3_down,'RelTol',1e-8,'AbsTol',1e-8);
[~,~,tp_out3,xp_out3,~] = ode15s(@relaxation,tspan,initialsp,options2,d1+delta,d2,d3,theta_I);

options2 = odeset('Events',@event1_up,'RelTol',1e-8,'AbsTol',1e-8);
[~,~,tp_in1,xp_in1,~] = ode15s(@relaxation,tspan,initialsp,options2,d1+delta,d2,d3,theta_I);
options2 = odeset('Events',@event1_down,'RelTol',1e-8,'AbsTol',1e-8);
[~,~,tp_out1,xp_out1,~] = ode15s(@relaxation,tspan,initialsp,options2,d1+delta,d2,d3,theta_I);
tp_out1 = tp_out1(2); xp_out1 = xp_out1(2,:);


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
% perturbed
dFp = @(v1) [-G.E*(v1-E.E)/C; 0; 0; 0; 0; 0];  %derivative of F w.r.t perturbation


%% Calculate lTRC (backward)

% During t_in1 to t_out1
tspan1 = [t_in1:dt:t_out1 t_out1];
[~,P1] = ode15s(@relaxation,t_in1:dt:t_out1,x_in1,options1,d1,d2,d3,theta_I);
P1 = [P1; x_out1];
v11 = P1(:,1); v21 = P1(:,2); v31 = P1(:,3); h11 = P1(:,4); h21 = P1(:,5); h31 = P1(:,6);
FF1 = zeros(6,length(tspan1)); DF1 = zeros(6,6,length(tspan1)); 
for i = 1:length(tspan1)
    FF1(:,i) = F(v11(i),v21(i),v31(i),h11(i),h21(i),h31(i));
    DF1(:,:,i) = A(v11(i),v21(i),v31(i),h11(i),h21(i),h31(i));
end

eta1 = zeros(6,length(tspan1)); 
eta1(:,end) = -n1_out/(n1_out'*FF1(:,end));
eta1(:,end-1) = eta1(:,end)+DF1(:,:,end)'*eta1(:,end)*(tspan1(end)-tspan1(end-1));
for i = length(tspan1)-2:-1:1
    eta1(:,i) = eta1(:,i+1)+DF1(:,:,i+1)'*eta1(:,i+1)*dt;
end

% During t_in2 to t_out2
tspan2 = [t_in2:dt:t_out2 t_out2];
[~,P2] = ode15s(@relaxation,t_in2:dt:t_out2,x_in2,options1,d1,d2,d3,theta_I);
P2 = [P2; x_out2];
v12 = P2(:,1); v22 = P2(:,2); v32 = P2(:,3); h12 = P2(:,4); h22 = P2(:,5); h32 = P2(:,6);
FF2 = zeros(6,length(tspan2)); DF2 = zeros(6,6,length(tspan2));
for i = 1:length(tspan2)
    FF2(:,i) = F(v12(i),v22(i),v32(i),h12(i),h22(i),h32(i));
    DF2(:,:,i) = A(v12(i),v22(i),v32(i),h12(i),h22(i),h32(i));
end

eta2 = zeros(6,length(tspan2)); 
eta2(:,end) = -n2_out/(n2_out'*FF2(:,end));
eta2(:,end-1) = eta2(:,end)+DF2(:,:,end)'*eta2(:,end)*(tspan2(end)-tspan2(end-1));
for i = length(tspan2)-2:-1:1
    eta2(:,i) = eta2(:,i+1)+DF2(:,:,i+1)'*eta2(:,i+1)*dt;
end

% During t_in3 to t_out3
tspan3 = [t_in3:dt:t_out3 t_out3];
[~,P3] = ode15s(@relaxation,t_in3:dt:t_out3,x_in3,options1,d1,d2,d3,theta_I);
P3 = [P3; x_out3];
v13 = P3(:,1); v23 = P3(:,2); v33 = P3(:,3); h13 = P3(:,4); h23 = P3(:,5); h33 = P3(:,6);
FF3=zeros(6,length(tspan3)); DF3=zeros(6,6,length(tspan3));
for i = 1:length(tspan3)
    FF3(:,i) = F(v13(i),v23(i),v33(i),h13(i),h23(i),h33(i));
    DF3(:,:,i) = A(v13(i),v23(i),v33(i),h13(i),h23(i),h33(i));
end

eta3 = zeros(6,length(tspan3)); 
eta3(:,end) = -n3_out/(n3_out'*FF3(:,end));
eta3(:,end-1) = eta3(:,end)+DF3(:,:,end)'*eta3(:,end)*(tspan3(end)-tspan3(end-1));
for i = length(tspan3)-2:-1:1
    eta3(:,i) = eta3(:,i+1)+DF3(:,:,i+1)'*eta3(:,i+1)*dt;
end


%% Shift in the active durations

int1 = zeros(length(tspan1),1); int2 = zeros(length(tspan2),1); int3 = zeros(length(tspan3),1);
for i = 1:length(tspan1)
    int1(i) = eta1(:,i)'*dFp(v11(i));
end
for i = 1:length(tspan2)
    int2(i) = eta2(:,i)'*dFp(v12(i));
end
for i = 1:length(tspan3)
    int3(i) = eta3(:,i)'*dFp(v13(i));
end

% linear shift
T1 = eta1(:,1)'*(xp_in1'-x_in1')/delta+dt*trapz(int1(1:end-1))+(tspan1(end)-tspan1(end-1))*trapz(int1(end-1:end));
T2 = eta2(:,1)'*(xp_in2'-x_in2')/delta+dt*trapz(int2(1:end-1))+(tspan2(end)-tspan2(end-1))*trapz(int2(end-1:end));
T3 = eta3(:,1)'*(xp_in3'-x_in3')/delta+dt*trapz(int3(1:end-1))+(tspan3(end)-tspan3(end-1))*trapz(int3(end-1:end));

% difference from unperturbed durations, using lTRC
dT1_lTRC = delta*T1;
dT2_lTRC = delta*T2;
dT3_lTRC = delta*T3;

% difference from unperturbed durations, direct calculation
dT1_direct = (tp_out1-tp_in1)-(t_out1-t_in1);
dT2_direct = (tp_out2-tp_in2)-(t_out2-t_in2);
dT3_direct = (tp_out3-tp_in3)-(t_out3-t_in3);


%% Figures and data

% integrand
T1_int = zeros(length(tspan1),1); T2_int = zeros(length(tspan2),1); T3_int = zeros(length(tspan3),1);
for i = 1:length(tspan1)
    T1_int(i) = eta1(:,1)'*(xp_in1'-x_in1')*dirac(i-1)/delta + eta1(:,i)'*dFp(v11(i));
end
for i = 1:length(tspan2)
    T2_int(i) = eta2(:,1)'*(xp_in2'-x_in2')*dirac(i-1)/delta + eta2(:,i)'*dFp(v12(i));
end
for i = 1:length(tspan3)
    T3_int(i) = eta3(:,1)'*(xp_in3'-x_in3')*dirac(i-1)/delta + eta3(:,i)'*dFp(v13(i));
end

% integral
T1_sum =  zeros(length(tspan1),1); T2_sum =  zeros(length(tspan2),1); T3_sum =  zeros(length(tspan3),1);
T1_sum(1) = eta1(:,1)'*(xp_in1'-x_in1')/delta;
T2_sum(1) = eta2(:,1)'*(xp_in2'-x_in2')/delta; 
T3_sum(1) = eta3(:,1)'*(xp_in3'-x_in3')/delta;
for i = 2:length(T1_sum)-1
    T1_sum(i) = T1_sum(i-1)+dt*trapz(int1(i-1:i));
end
for i = 2:length(T2_sum)-1
    T2_sum(i) = T2_sum(i-1)+dt*trapz(int2(i-1:i));
end
for i = 2:length(T3_sum)-1
    T3_sum(i) = T3_sum(i-1)+dt*trapz(int3(i-1:i));
end
T1_sum(end) = T1_sum(end-1)+(tspan1(end)-tspan1(end-1))*trapz(int1(end-1:end));
T2_sum(end) = T2_sum(end-1)+(tspan2(end)-tspan2(end-1))*trapz(int2(end-1:end));
T3_sum(end) = T3_sum(end-1)+(tspan3(end)-tspan3(end-1))*trapz(int3(end-1:end));

figure(5)
subplot(3,1,1)  % v trajectories
plot(TTu,v1,'*k','MarkerSize',1); hold on
plot(TTu,v2,'*b','MarkerSize',1); 
plot(TTu,v3,'*r','MarkerSize',1); 
plot(TTu,theta_I*ones(length(TTu),1),'--m','LineWidth',1.5); hold off
grid on
ylabel('v (mV)'); xlim([t_in2-1 t_out1+1]); ylim([-70 20]);
title('Synaptic release')
set(gca,'FontSize',13)
subplot(3,1,2)  % integrand
plot(tspan2,T2_int,'-b','LineWidth',2); hold on
if T2_sum(1)<0
    plot(t_in2*ones(1,length([T2_sum(1):.01:0])),[T2_sum(1):.01:0],'-b','LineWidth',2)
else
    plot(t_in2*ones(1,length([0:.01:T2_sum(1)])),[0:.01:T2_sum(1)],'-b','LineWidth',2)
end
plot(tspan3,T3_int,'-r','LineWidth',2)
if T3_sum(1)<0
    plot(t_in3*ones(1,length([T3_sum(1):.01:0])),[T3_sum(1):.01:0],'-r','LineWidth',2)
else
    plot(t_in3*ones(1,length([0:.01:T3_sum(1)])),[0:.01:T3_sum(1)],'-r','LineWidth',2)
end
plot(tspan1,T1_int,'-k','LineWidth',2); 
if T1_sum(1)<0
    plot(t_in1*ones(1,length([T1_sum(1):.01:0])),[T1_sum(1):.01:0],'-k','LineWidth',2)
else
    plot(t_in1*ones(1,length([0:.01:T1_sum(1)])),[0:.01:T1_sum(1)],'-k','LineWidth',2)
end
hold off
grid on
ylabel('integrand'); xlim([t_in2-1 t_out1+1]);
set(gca,'FontSize',13)
subplot(3,1,3)  % integral
plot(tspan2,T2_sum,'-b','LineWidth',2); hold on
plot(tspan3,T3_sum,'-r','LineWidth',2);
plot(tspan1,T1_sum,'-k','LineWidth',2); hold off
grid on
xlabel('time (ms)'); ylabel('intergal'); xlim([t_in2-1 t_out1+1]);
set(gca,'FontSize',13)

% data
fprintf('\n');
fprintf('\t\t\t\t\t t_1 \t\t\t t_2 \t\t\t t_3 \n');
fprintf('d_1=1: \t\t\t\t %0.4f \t\t %0.4f \t\t %0.4f \n', t_out1-t_in1, t_out2-t_in2, t_out3-t_in3);
fprintf('d_1=1.05: \t\t\t %0.4f \t\t %0.4f \t\t %0.4f \n', tp_out1-tp_in1, tp_out2-tp_in2, tp_out3-tp_in3);
fprintf('\n');
fprintf('\t\t\t\t\t dt_1 \t\t\t dt_2 \t\t\t dt_3 \n');
fprintf('Direct difference:\t %0.4f \t\t %0.4f \t\t %0.4f \n', dT1_direct, dT2_direct, dT3_direct);
fprintf('lTRC difference:\t %0.4f \t\t %0.4f \t\t %0.4f \n', dT1_lTRC, dT2_lTRC, dT3_lTRC);