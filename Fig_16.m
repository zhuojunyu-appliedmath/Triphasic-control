clear;

%% Parameters

a = 0.01; ap = 0.0105;
a1 = ap; a2 = a; a3 = a;
rho = 3;
delta = ap-a;  % perturbation

%% Limit cycle

tF = 30; dt = 0.001; tspan = 0:dt:tF;
% unperturbed
initials = [0.500342977762597 0.0172224538514646 0.237135024241301];
Pu = Aplysia(initials,0,tF,dt,a,a,a,rho);
X = Pu(:,1); Y = Pu(:,2); Z = Pu(:,3);
% perturbed
initialsp = [0.500389004636971,0.0173951947597979,0.236181333160755];
Pp = Aplysia(initialsp,0,tF,dt,a1,a2,a3,rho);
Xp = Pp(:,1); Yp = Pp(:,2); Zp = Pp(:,3);


%% Switching surfaces and normal direction

n1_out = [-1/sqrt(2); 1/sqrt(2); 0]; %R1_out: x-y-a=0
n2_out = [0; -1/sqrt(2); 1/sqrt(2)]; %R2_out: y-z-a=0
n3_out = [1/sqrt(2); 0; -1/sqrt(2)]; %R3_out: x-z+a=0


%% Find enter/exit time 

% unperturbed
tt = []; 

for j = 1:length(tspan)-1
    
    if ((X(j)<Y(j)+a) || (X(j)<Z(j)-a)) && ((X(j+1)>=Y(j+1)+a) && (X(j+1)>=Z(j+1)-a))
        tt = [tt; (j+1)*dt 11];  % 1 up
    end
    if ((X(j)>=Y(j)+a) && (X(j)>=Z(j)-a)) && ((X(j+1)<Y(j+1)+a) || (X(j+1)<Z(j+1)-a))
        tt = [tt; j*dt 12];  % 1 down
    end
    if ((Y(j)<=X(j)-a) || (Y(j)<Z(j)+a)) && ((Y(j+1)>X(j+1)-a) && (Y(j+1)>=Z(j+1)+a))
        tt = [tt; (j+1)*dt 21];  % 2 up
    end
    if ((Y(j)>X(j)-a) && (Y(j)>=Z(j)+a)) && ((Y(j+1)<=X(j+1)-a) || (Y(j+1)<Z(j+1)+a))
         tt = [tt; j*dt 22];  % 2 down
    end
    if ((Z(j)<=X(j)+a) || (Z(j)<=Y(j)-a)) && ((Z(j+1)>X(j+1)+a) && (Z(j+1)>Y(j+1)-a))
        tt = [tt; (j+1)*dt 31];  % 3 up
    end
    if ((Z(j)>X(j)+a) && (Z(j)>Y(j)-a)) && ((Z(j+1)<=X(j+1)+a) || (Z(j+1)<=Y(j+1)-a))
        tt = [tt; j*dt 32];  % 3 down
    end
    
end

i = find(tt(:,2)==21); j = find(tt(:,2)==22);
t_in2 = tt(i(1),1); t_out2 = tt(j(1),1);
x_in2 = Pu(t_in2/dt,:); x_out2 = Pu(t_out2/dt,:);
i = find(tt(:,2)==31); j = find(tt(:,2)==32);
t_in3 = tt(i(1),1); t_out3 = tt(j(1),1);
x_in3 = Pu(t_in3/dt,:); x_out3 = Pu(t_out3/dt,:);
i = find(tt(:,2)==11); j = find(tt(:,2)==12);
t_in1 = tt(i(1),1); t_out1 = tt(j(2),1);
x_in1 = Pu(t_in1/dt,:); x_out1 = Pu(t_out1/dt,:);

% perturbed
ttp = []; 

for j = 1:length(tspan)-1
    
    if ((Xp(j)<Yp(j)+(a1+a2)/2) || (Xp(j)<Zp(j)-(a1+a3)/2)) && ((Xp(j+1)>=Yp(j+1)+(a1+a2)/2) && (Xp(j+1)>=Zp(j+1)-(a1+a3)/2))
        ttp = [ttp; (j+1)*dt 11];  % 1 up
    end
    if ((Xp(j)>=Yp(j)+(a1+a2)/2) && (Xp(j)>=Zp(j)-(a1+a3)/2)) && ((Xp(j+1)<Yp(j+1)+(a1+a2)/2) || (Xp(j+1)<Zp(j+1)-(a1+a3)/2))
        ttp = [ttp; j*dt 12];  % 1 down
    end
    if ((Yp(j)<=Xp(j)-(a1+a2)/2) || (Yp(j)<Zp(j)+(a2+a3)/2)) && ((Yp(j+1)>Xp(j+1)-(a1+a2)/2) && (Yp(j+1)>=Zp(j+1)+(a2+a3)/2))
        ttp = [ttp; (j+1)*dt 21];  % 2 up
    end
    if ((Yp(j)>Xp(j)-(a1+a2)/2) && (Yp(j)>=Zp(j)+(a2+a3)/2)) && ((Yp(j+1)<=Xp(j+1)-(a1+a2)/2) || (Yp(j+1)<Zp(j+1)+(a1+a3)/2))
         ttp = [ttp; j*dt 22];  % 2 down
    end
    if ((Zp(j)<=Xp(j)+(a1+a3)/2) || (Zp(j)<=Yp(j)-(a2+a3)/2)) && ((Zp(j+1)>Xp(j+1)+(a1+a3)/2) && (Zp(j+1)>Yp(j+1)-(a2+a3)/2))
        ttp = [ttp; (j+1)*dt 31];  % 3 up
    end
    if ((Zp(j)>Xp(j)+(a1+a3)/2) && (Zp(j)>Yp(j)-(a2+a3)/2)) && ((Zp(j+1)<=Xp(j+1)+(a1+a3)/2) || (Zp(j+1)<=Yp(j+1)-(a2+a3)/2))
        ttp = [ttp; j*dt 32];  % 3 down
    end
    
end

i = find(ttp(:,2)==21); j = find(ttp(:,2)==22);
tp_in2 = ttp(i(1),1); tp_out2 = ttp(j(1),1);
xp_in2 = Pp(tp_in2/dt,:); xp_out2 = Pp(tp_out2/dt,:);
i = find(ttp(:,2)==31); j = find(ttp(:,2)==32);
tp_in3 = ttp(i(1),1); tp_out3 = ttp(j(1),1);
xp_in3 = Pp(tp_in3/dt,:); xp_out3 = Pp(tp_out3/dt,:);
i = find(ttp(:,2)==11); j = find(ttp(:,2)==12);
tp_in1 = ttp(i(1),1); tp_out1 = ttp(j(2),1);
xp_in1 = Pp(tp_in1/dt,:); xp_out1 = Pp(tp_out1/dt,:);


%% Vector field and Jacobian

A1 = [-1 -rho 0; 0 1 0; 0 0 1-rho];
A2 = [1-rho 0 0; 0 -1 -rho; 0 0 1];
A3 = [1 0 0; 0 1-rho 0; -rho 0 -1];

FF1 = @(x,y,z) [1-x-(y+a)*rho; y+a; (z-a)*(1-rho)];
FF2 = @(x,y,z) [(x-a)*(1-rho); 1-y-(z+a)*rho; z+a];
FF3 = @(x,y,z) [x+a; (y-a)*(1-rho); 1-z-(x+a)*rho];

dFp1 = [-rho; 0; 0];
dFp2 = [rho-1; 0; 0];
dFp3 = [1; 0; 0];


%% Calculate lTRC

% During t_in1 to t_out1
tspan1 = [t_in1:dt:t_out1];
eta1 = zeros(3,length(tspan1)); 
eta1(:,end) = -n1_out/(n1_out'*FF1(x_out1(1),x_out1(2),x_out1(3)));
Phi1 = @(t) [exp(t-t_out1)  0  0; 
             rho*(exp(t-t_out1)-exp(t_out1-t))/2  exp(t_out1-t)  0;
             0  0  exp((rho-1)*(t-t_out1))];   % exponential of -DF1^T
for i = 1:length(tspan1)-1
    eta1(:,i) = Phi1(tspan1(i))*eta1(:,end);
end

% During t_in2 to t_out2
tspan2 = [t_in2:dt:t_out2];
eta2 = zeros(3,length(tspan2)); 
eta2(:,end) = -n2_out/(n2_out'*FF2(x_out2(1),x_out2(2),x_out2(3)));
Phi2 = @(t) [exp((rho-1)*(t-t_out2))  0  0; 
             0  exp(t-t_out2)  0;
             0  rho*(exp(t-t_out2)-exp(t_out2-t))/2  exp(t_out2-t)];   % exponential of -DF1^T
for i = 1:length(tspan2)-1
    eta2(:,i) = Phi2(tspan2(i))*eta2(:,end);
end

% During t_in3 to t_out3
tspan3 = [t_in3:dt:t_out3];
eta3 = zeros(3,length(tspan3)); 
eta3(:,end) = -n3_out/(n3_out'*FF3(x_out3(1),x_out3(2),x_out3(3)));
Phi3 = @(t) [exp(t_out3-t)  0  rho*(exp(t-t_out3)-exp(t_out3-t))/2; 
             0  exp((rho-1)*(t-t_out3))  0;
             0  0  exp(t-t_out3)];   % exponential of -DF1^T
for i = 1:length(tspan3)-1
    eta3(:,i) = Phi3(tspan3(i))*eta3(:,end);
end

%% Shift in the active durations

int1 = zeros(length(tspan1),1); int2 = zeros(length(tspan2),1); int3 = zeros(length(tspan3),1);
for i = 1:length(tspan1)
    int1(i) = eta1(:,i)'*dFp1;
end
for i = 1:length(tspan2)
    int2(i) = eta2(:,i)'*dFp2;
end
for i = 1:length(tspan3)
    int3(i) = eta3(:,i)'*dFp3;
end

% linear shift
T1 = eta1(:,1)'*(xp_in1'-x_in1')/delta-eta1(:,end)'*(xp_out1'-x_out1')/delta+dt*trapz(int1(1:end));
T2 = eta2(:,1)'*(xp_in2'-x_in2')/delta-eta2(:,end)'*(xp_out2'-x_out2')/delta+dt*trapz(int2(1:end));
T3 = eta3(:,1)'*(xp_in3'-x_in3')/delta-eta3(:,end)'*(xp_out3'-x_out3')/delta+dt*trapz(int3(1:end));

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
    T1_int(i) = eta1(:,1)'*(xp_in1'-x_in1')*dirac(i-1)/delta - eta1(:,end)'*(xp_out1'-x_out1')*dirac(i-length(tspan1))/delta + eta1(:,i)'*dFp1;
end
for i = 1:length(tspan2)
    T2_int(i) = eta2(:,1)'*(xp_in2'-x_in2')*dirac(i-1)/delta - eta2(:,end)'*(xp_out2'-x_out2')*dirac(i-length(tspan2))/delta + eta2(:,i)'*dFp2;
end
for i = 1:length(tspan3)
    T3_int(i) = eta3(:,1)'*(xp_in3'-x_in3')*dirac(i-1)/delta - eta3(:,end)'*(xp_out3'-x_out3')*dirac(i-length(tspan3))/delta + eta3(:,i)'*dFp3;
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
T1_sum(end) = T1_sum(end-1) - eta1(:,end)'*(xp_out1'-x_out1')/delta;
T2_sum(end) = T2_sum(end-1) - eta2(:,end)'*(xp_out2'-x_out2')/delta;
T3_sum(end) = T3_sum(end-1) - eta3(:,end)'*(xp_out3'-x_out3')/delta;

%% Figure and table
figure
subplot(3,1,1)  % trajectories
plot(tspan,X,'*k','MarkerSize',1); hold on
plot(tspan,Y,'*b','MarkerSize',1); 
plot(tspan,Z,'*r','MarkerSize',1); hold off 
grid on
ylabel('x, y, z'); xlim([t_in2-dt t_out1+dt]); 
set(gca,'FontSize',13);
subplot(3,1,2)  % integrand
plot(tspan2,T2_int,'-b','LineWidth',2); hold on
if T2_int(1)<0
    plot(t_in2*ones(1,length([T2_sum(1):.01:0])),[T2_sum(1):.01:0],'-b','LineWidth',2)
else
    plot(t_in2*ones(1,length([0:.01:T2_sum(1)])),[0:.01:T2_sum(1)],'-b','LineWidth',2)
end
if T2_int(end)<0
    plot(t_out2*ones(1,length([-eta2(:,end)'*(xp_out2'-x_out2')/delta:.01:0])),[-eta2(:,end)'*(xp_out2'-x_out2')/delta:.01:0],'-b','LineWidth',2)
else
    plot(t_out2*ones(1,length([0:.01:-eta2(:,end)'*(xp_out2'-x_out2')/delta])),[0:.01:-eta2(:,end)'*(xp_out2'-x_out2')/delta],'-b','LineWidth',2)
end
plot(tspan3,T3_int,'-r','LineWidth',2)
if T3_int(1)<0
    plot(t_in3*ones(1,length([T3_sum(1):.01:0])),[T3_sum(1):.01:0],'-r','LineWidth',2)
else
    plot(t_in3*ones(1,length([0:.01:T3_sum(1)])),[0:.01:T3_sum(1)],'-r','LineWidth',2)
end
if T3_int(end)<0
    plot(t_out3*ones(1,length([-eta3(:,end)'*(xp_out3'-x_out3')/delta:.01:0])),[-eta3(:,end)'*(xp_out3'-x_out3')/delta:.01:0],'-r','LineWidth',2)
else
    plot(t_out3*ones(1,length([0:.01:-eta3(:,end)'*(xp_out3'-x_out3')/delta])),[0:.01:-eta3(:,end)'*(xp_out3'-x_out3')/delta],'-r','LineWidth',2)
end
plot(tspan1,T1_int,'-k','LineWidth',2); 
if T1_int(1)<0
    plot(t_in1*ones(1,length([T1_sum(1):.01:0])),[T1_sum(1):.01:0],'-k','LineWidth',2)
else
    plot(t_in1*ones(1,length([0:.01:T1_sum(1)])),[0:.01:T1_sum(1)],'-k','LineWidth',2)
end
if T1_int(end)<0
    plot(t_out1*ones(1,length([-eta1(:,end)'*(xp_out1'-x_out1')/delta:.01:0])),[-eta1(:,end)'*(xp_out1'-x_out1')/delta:.01:0],'-k','LineWidth',2)
else
    plot(t_out1*ones(1,length([0:.01:-eta1(:,end)'*(xp_out1'-x_out1')/delta])),[0:.01:-eta1(:,end)'*(xp_out1'-x_out1')/delta],'-k','LineWidth',2)
end
hold off; grid on
ylabel('integrand'); xlim([t_in2-dt t_out1+dt]);
set(gca,'FontSize',13);
subplot(3,1,3)  % integral
plot(tspan2,T2_sum,'-b','LineWidth',2); hold on
plot(tspan3,T3_sum,'-r','LineWidth',2);
plot(tspan1,T1_sum,'-k','LineWidth',2); hold off
grid on; 
xlabel('time'); ylabel('intergal'); xlim([t_in2-dt t_out1+dt]);
set(gca,'FontSize',13);

% data
fprintf('\n');
fprintf('\t\t\t\t\t t_1 \t\t\t t_2 \t\t\t t_3 \n');
fprintf('a1=0.01: \t\t\t %0.4f \t\t %0.4f \t\t %0.4f \n', t_out1-t_in1, t_out2-t_in2, t_out3-t_in3);
fprintf('a1=0.0105: \t\t\t %0.4f \t\t %0.4f \t\t %0.4f \n', tp_out1-tp_in1, tp_out2-tp_in2, tp_out3-tp_in3);
fprintf('\n');
fprintf('\t\t\t\t\t dt_1 \t\t\t dt_2 \t\t\t dt_3 \n');
fprintf('Direct difference:\t %0.4f \t\t %0.4f \t\t %0.4f \n', dT1_direct, dT2_direct, dT3_direct);
fprintf('lTRC difference:\t %0.4f \t\t %0.4f \t\t %0.4f \n', dT1_lTRC, dT2_lTRC, dT3_lTRC);