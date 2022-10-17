% Fig. 18: lTRC for the competitive threshold-linear model

clear; clc;

%% Parameters

epsilon = 0.25; 
delta = 0.5;
theta = 1; thetap = 1.01;
theta1 = thetap; theta2 = theta; theta3 = theta;
dtheta = thetap-theta;  %perturbation

%% Limit cycle

tF = 20; dt = 0.001; tspan = 0:dt:tF;
% unperturbed
initials = [0.501193190683074,0.0139591696003762,0.425786237781174];
[~,Pu] = ode45(@threshold_linear,tspan,initials,[],epsilon,delta,theta,theta,theta);
x1 = Pu(:,1); x2 = Pu(:,2); x3 = Pu(:,3);
% perturbed
initialsp = [0.500055603116643,0.0160115420011988,0.428565412506380];
[~,Pp] = ode45(@threshold_linear,tspan,initialsp,[],epsilon,delta,thetap,theta,theta);
x1p = Pp(:,1); x2p = Pp(:,2); x3p = Pp(:,3);


%% Switching surfaces and normal direction

n1_out = [-1/sqrt(2); 1/sqrt(2); 0]; %R1_out: x1-x2=0
n2_out = [0; -1/sqrt(2); 1/sqrt(2)]; %R2_out: x2-x3=0
n3_out = [1/sqrt(2); 0; -1/sqrt(2)]; %R3_out: x1-x3=0


%% Find enter/exit time 

% unperturbed
tt = []; 

for j = 1:length(tspan)-1
    
    if ((x1(j)<x2(j)) || (x1(j)<x3(j))) && ((x1(j+1)>=x2(j+1)) && (x1(j+1)>=x3(j+1)))
        tt = [tt; (j+1)*dt 11];  % 1 up
    end
    if ((x1(j)>=x2(j)) && (x1(j)>=x3(j))) && ((x1(j+1)<x2(j+1)) || (x1(j+1)<x3(j+1)))
        tt = [tt; j*dt 12];  % 1 down
    end
    if ((x2(j)<=x1(j)) || (x2(j)<x3(j))) && ((x2(j+1)>x1(j+1)) && (x2(j+1)>=x3(j+1)))
        tt = [tt; (j+1)*dt 21];  % 2 up
    end
    if ((x2(j)>x1(j)) && (x2(j)>=x3(j))) && ((x2(j+1)<=x1(j+1)) || (x2(j+1)<x3(j+1)))
         tt = [tt; j*dt 22];  % 2 down
    end
    if ((x3(j)<=x1(j)) || (x3(j)<=x2(j))) && ((x3(j+1)>x1(j+1)) && (x3(j+1)>x2(j+1)))
        tt = [tt; (j+1)*dt 31];  % 3 up
    end
    if ((x3(j)>x1(j)) && (x3(j)>x2(j))) && ((x3(j+1)<=x1(j+1)) || (x3(j+1)<=x2(j+1)))
        tt = [tt; j*dt 32];  % 3 down
    end
    
end

i = find(tt(:,2)==21); j = find(tt(:,2)==22);
t_in2 = tt(i(1),1); t_out2 = tt(j(1),1);
x_in2 = Pu(t_in2/dt,:); x_out2 = Pu(t_out2/dt,:);
i = find(tt(:,2)==31); j = find(tt(:,2)==32);
t_in3 = tt(i(1),1); t_out3 = tt(j(1),1);
x_in3 = Pu(floor(t_in3/dt),:); x_out3 = Pu(t_out3/dt,:);
i = find(tt(:,2)==11); j = find(tt(:,2)==12);
t_in1 = tt(i(1),1); t_out1 = tt(j(2),1);
x_in1 = Pu(floor(t_in1/dt),:); x_out1 = Pu(floor(t_out1/dt),:);

% perturbed
ttp = []; 

for j = 1:length(tspan)-1
    
    if ((x1p(j)<x2p(j)) || (x1p(j)<x3p(j))) && ((x1p(j+1)>=x2p(j+1)) && (x1p(j+1)>=x3p(j+1)))
        ttp = [ttp; (j+1)*dt 11];  % 1 up
    end
    if ((x1p(j)>=x2p(j)) && (x1p(j)>=x3p(j))) && ((x1p(j+1)<x2p(j+1)) || (x1p(j+1)<x3p(j+1)))
        ttp = [ttp; j*dt 12];  % 1 down
    end
    if ((x2p(j)<=x1p(j)) || (x2p(j)<x3p(j))) && ((x2p(j+1)>x1p(j+1)) && (x2p(j+1)>=x3p(j+1)))
        ttp = [ttp; (j+1)*dt 21];  % 2 up
    end
    if ((x2p(j)>x1p(j)) && (x2p(j)>=x3p(j))) && ((x2p(j+1)<=x1p(j+1)) || (x2p(j+1)<x3p(j+1)))
        ttp = [ttp; j*dt 22];  % 2 down
    end
    if ((x3p(j)<=x1p(j)) || (x3p(j)<=x2p(j))) && ((x3p(j+1)>x1p(j+1)) && (x3p(j+1)>x2p(j+1)))
        ttp = [ttp; (j+1)*dt 31];  % 3 up
    end
    if ((x3p(j)>x1p(j)) && (x3p(j)>x2p(j))) && ((x3p(j+1)<=x1p(j+1)) || (x3p(j+1)<=x2p(j+1)))
        ttp = [ttp; j*dt 32];  % 3 down
    end
    
end

i = find(ttp(:,2)==21); j = find(ttp(:,2)==22);
tp_in2 = ttp(i(1),1); tp_out2 = ttp(j(1),1);
xp_in2 = Pp(tp_in2/dt,:); xp_out2 = Pp(tp_out2/dt,:);
i = find(ttp(:,2)==31); j = find(ttp(:,2)==32);
tp_in3 = ttp(i(1),1); tp_out3 = ttp(j(1),1);
xp_in3 = Pp(floor(tp_in3/dt),:); xp_out3 = Pp(tp_out3/dt,:);
i = find(ttp(:,2)==11); j = find(ttp(:,2)==12);
tp_in1 = ttp(i(1),1); tp_out1 = ttp(j(2),1);
xp_in1 = Pp(floor(tp_in1/dt),:); xp_out1 = Pp(floor(tp_out1/dt),:);


%% Vector field and Jacobian

W13 = -1+epsilon; W21 = -1+epsilon; W32 = -1+epsilon; 
W12 = -1-delta; W23 = -1-delta; W31 = -1-delta; 

F = @(x1,x2,x3) [-x1+max(W13*x3+W12*x2+theta1,0); -x2+max(W21*x1+W23*x3+theta2,0); -x3+max(W32*x2+W31*x1+theta3,0)];


%% Calculate lTRC

% During t_in1 to t_out1
tspan1 = [t_in1:dt:t_out1 t_out1];
[~,P1] = ode45(@threshold_linear,t_in1:dt:t_out1,x_in1,[],epsilon,delta,theta,theta,theta);
P1 = [P1; x_out1];
x11 = P1(:,1); x21 = P1(:,2); x31 = P1(:,3);
FF1 = zeros(3,length(tspan1)); DF1 = zeros(3,3,length(tspan1)); 
for i = 1:length(tspan1)
    FF1(:,i) = F(x11(i),x21(i),x31(i));
    DF1(:,:,i)  = Jacobian_TL(x11(i),x21(i),x31(i),epsilon,delta,theta,theta,theta);
end

eta1 = zeros(3,length(tspan1)); 
eta1(:,end) = -n1_out/(n1_out'*FF1(:,end));
eta1(:,end-1) = eta1(:,end)+DF1(:,:,end)'*eta1(:,end)*(tspan1(end)-tspan1(end-1));
for i = length(tspan1)-2:-1:1
    eta1(:,i) = eta1(:,i+1)+DF1(:,:,i+1)'*eta1(:,i+1)*dt;
end

% During t_in2 to t_out2
tspan2 = [t_in2:dt:t_out2 t_out2];
[~,P2] = ode45(@threshold_linear,t_in2:dt:t_out2,x_in2,[],epsilon,delta,theta,theta,theta);
P2 = [P2; x_out2];
x12 = P2(:,1); x22 = P2(:,2); x32 = P2(:,3);
FF2 = zeros(3,length(tspan2)); DF2 = zeros(3,3,length(tspan2)); 
for i = 1:length(tspan2)
    FF2(:,i) = F(x12(i),x22(i),x32(i));
    DF2(:,:,i)  = Jacobian_TL(x12(i),x22(i),x32(i),epsilon,delta,theta,theta,theta);
end

eta2 = zeros(3,length(tspan2)); 
eta2(:,end) = -n2_out/(n2_out'*FF2(:,end));
eta2(:,end-1) = eta2(:,end)+DF2(:,:,end)'*eta2(:,end)*(tspan2(end)-tspan2(end-1));
for i = length(tspan2)-2:-1:1
    eta2(:,i) = eta2(:,i+1)+DF2(:,:,i+1)'*eta2(:,i+1)*dt;
end

% During t_in3 to t_out3
tspan3 = [t_in3:dt:t_out3 t_out3];
[~,P3] = ode45(@threshold_linear,t_in3:dt:t_out3,x_in3,[],epsilon,delta,theta,theta,theta);
P3 = [P3; x_out3];
x13 = P3(:,1); x23 = P3(:,2); x33 = P3(:,3);
FF3 = zeros(3,length(tspan3)); DF3 = zeros(3,3,length(tspan3));
for i = 1:length(tspan3)
    FF3(:,i) = F(x13(i),x23(i),x33(i));
    DF3(:,:,i)  = Jacobian_TL(x13(i),x23(i),x33(i),epsilon,delta,theta,theta,theta);
end

eta3 = zeros(3,length(tspan3)); 
eta3(:,end) = -n3_out/(n3_out'*FF3(:,end));
eta3(:,end-1) = eta3(:,end)+DF3(:,:,end)'*eta3(:,end)*(tspan3(end)-tspan3(end-1));
for i = length(tspan3)-2:-1:1
    eta3(:,i) = eta3(:,i+1)+DF3(:,:,i+1)'*eta3(:,i+1)*dt;
end


%% Shift in the active durations

int1 = zeros(length(tspan1),1); int2 = zeros(length(tspan2),1); int3 = zeros(length(tspan3),1);
dFp1 = zeros(3,length(tspan1)); dFp2 = zeros(3,length(tspan2)); dFp3 = zeros(3,length(tspan3));
for i = 1:length(tspan1)
    dFp1(:,i) = dFp_TL(x21(i),x31(i),epsilon,delta,theta);
    int1(i) = eta1(:,i)'*dFp1(:,i);
end
for i = 1:length(tspan2)
    dFp2(:,i) = dFp_TL(x22(i),x32(i),epsilon,delta,theta);
    int2(i) = eta2(:,i)'*dFp2(:,i);
end
for i = 1:length(tspan3)
    dFp3(:,i) = dFp_TL(x23(i),x33(i),epsilon,delta,theta);
    int3(i) = eta3(:,i)'*dFp3(:,i);
end

% linear shift
T1 = eta1(:,1)'*(xp_in1'-x_in1')/dtheta+dt*trapz(int1(1:end));
T2 = eta2(:,1)'*(xp_in2'-x_in2')/dtheta+dt*trapz(int2(1:end));
T3 = eta3(:,1)'*(xp_in3'-x_in3')/dtheta+dt*trapz(int3(1:end));

% difference from unperturbed durations, using lTRC
dT1_lTRC = dtheta*T1;
dT2_lTRC = dtheta*T2;
dT3_lTRC = dtheta*T3;

% difference from unperturbed durations, direct calculation
dT1_direct = (tp_out1-tp_in1)-(t_out1-t_in1);
dT2_direct = (tp_out2-tp_in2)-(t_out2-t_in2);
dT3_direct = (tp_out3-tp_in3)-(t_out3-t_in3);


%% Figures and data

% integrand
T1_int = zeros(length(tspan1),1); T2_int = zeros(length(tspan2),1); T3_int = zeros(length(tspan3),1);
for i = 1:length(tspan1)
    T1_int(i) = eta1(:,1)'*(xp_in1'-x_in1')*dirac(i-1)/dtheta + eta1(:,i)'*dFp1(:,i);
end
for i = 1:length(tspan2)
    T2_int(i) = eta2(:,1)'*(xp_in2'-x_in2')*dirac(i-1)/dtheta + eta2(:,i)'*dFp2(:,i);
end
for i = 1:length(tspan3)
    T3_int(i) = eta3(:,1)'*(xp_in3'-x_in3')*dirac(i-1)/dtheta + eta3(:,i)'*dFp3(:,i);
end

% integral
T1_sum = zeros(length(tspan1),1); T2_sum = zeros(length(tspan2),1); T3_sum = zeros(length(tspan3),1);
T1_sum(1) = eta1(:,1)'*(xp_in1'-x_in1')/dtheta;
T2_sum(1) = eta2(:,1)'*(xp_in2'-x_in2')/dtheta; 
T3_sum(1) = eta3(:,1)'*(xp_in3'-x_in3')/dtheta;
for i = 2:length(T1_sum)
    T1_sum(i) = T1_sum(i-1)+dt*trapz(int1(i-1:i));
end
for i = 2:length(T2_sum)
    T2_sum(i) = T2_sum(i-1)+dt*trapz(int2(i-1:i));
end
for i = 2:length(T3_sum)
    T3_sum(i) = T3_sum(i-1)+dt*trapz(int3(i-1:i));
end

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


figure(7)
subplot(3,1,1)  % trajectories
plot(tspan,x1,'*k','MarkerSize',1); hold on
plot(tspan,x2,'*b','MarkerSize',1); 
plot(tspan,x3,'*r','MarkerSize',1); 
plot(tspan(nonsmoothpt1),x1(nonsmoothpt1),'*k','MarkerSize',8)
plot(tspan(nonsmoothpt2),x2(nonsmoothpt2),'*b','MarkerSize',8)
plot(tspan(nonsmoothpt3),x3(nonsmoothpt3),'*r','MarkerSize',8); hold off
grid on
ylabel('x_1, x_2, x_3'); xlim([t_in2-dt t_out1+dt]); 
set(gca,'FontSize',13);
subplot(3,1,2)  % integrand
plot(tspan2,T2_int,'-b','LineWidth',2); hold on
if T2_int(1)<0
    plot(t_in2*ones(1,length([T2_sum(1):.01:0])),[T2_sum(1):.01:0],'-b','LineWidth',2)
else
    plot(t_in2*ones(1,length([0:.01:T2_sum(1)])),[0:.01:T2_sum(1)],'-b','LineWidth',2)
end
plot(tspan3,T3_int,'-r','LineWidth',2)
if T3_int(1)<0
    plot(t_in3*ones(1,length([T3_sum(1):.01:0])),[T3_sum(1):.01:0],'-r','LineWidth',2)
else
    plot(t_in3*ones(1,length([0:.01:T3_sum(1)])),[0:.01:T3_sum(1)],'-r','LineWidth',2)
end
plot(tspan1,T1_int,'-k','LineWidth',2); 
if T1_int(1)<0
    plot(t_in1*ones(1,length([T1_sum(1):.01:0])),[T1_sum(1):.01:0],'-k','LineWidth',2)
else
    plot(t_in1*ones(1,length([0:.01:T1_sum(1)])),[0:.01:T1_sum(1)],'-k','LineWidth',2)
end
hold off; grid on
ylabel('integrand'); xlim([t_in2-dt t_out1+dt]);
set(gca,'FontSize',13);
subplot(3,1,3)  % integral
plot(tspan2,T2_sum,'-b','LineWidth',2); hold on
plot(tspan3,T3_sum,'-r','LineWidth',2);
plot(tspan1,T1_sum,'-k','LineWidth',2); hold off
grid on
xlabel('time'); ylabel('intergal'); xlim([t_in2-dt t_out1+dt]);
set(gca,'FontSize',13);

% data
fprintf('\n');
fprintf('\t\t\t\t\t t_1 \t\t\t t_2 \t\t\t t_3 \n');
fprintf('theta1=1: \t\t\t %0.4f \t\t %0.4f \t\t %0.4f \n', t_out1-t_in1, t_out2-t_in2, t_out3-t_in3);
fprintf('theta1=1.01: \t\t %0.4f \t\t %0.4f \t\t %0.4f \n', tp_out1-tp_in1, tp_out2-tp_in2, tp_out3-tp_in3);
fprintf('\n');
fprintf('\t\t\t\t\t dt_1 \t\t\t dt_2 \t\t\t dt_3 \n');
fprintf('Direct difference:\t %0.4f \t\t %0.4f \t\t %0.4f \n', dT1_direct, dT2_direct, dT3_direct);
fprintf('lTRC difference:\t %0.4f \t\t %0.4f \t\t %0.4f \n', dT1_lTRC, dT2_lTRC, dT3_lTRC);