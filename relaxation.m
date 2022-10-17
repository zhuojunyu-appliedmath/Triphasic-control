function dPdt = relaxation(t,P,d1,d2,d3,theta_I)

C = 0.21; 
epsilon = 0.01;
E = struct('Na', 50, 'L', -65, 'I', -80, 'E', 0); % reversal potentials, mV
sigma_I = -0.01;
B = struct('b12', 1, 'b13', 1, 'b21', 1, 'b23', 1, 'b31', 1, 'b32', 1); % coupling constants
G = struct('NaP', 6.8, 'L', 3, 'I', 0.4, 'E', 0.1); % conductances, nS
Theta = struct('h', -40, 'mp', -37); % half activations (mV)
Sigma = struct('h', 6, 'mp', -6); % slopes

v1 = P(1); v2 = P(2); v3 = P(3);
h1 = P(4); h2 = P(5); h3 = P(6);

S_inf = @(v) 1/(1+exp((v-theta_I)/sigma_I));
h_inf = @(v) 1/(1+exp((v-Theta.h)/Sigma.h));
mp_inf = @(v) 1/(1+exp((v-Theta.mp)/Sigma.mp));
tauh = @(v) epsilon*cosh((v-Theta.h)/Sigma.h/2);

I_NaP = @(v,h) G.NaP*mp_inf(v)*h*(v-E.Na);
I_L = @(v) G.L*(v-E.L);

F1 = -(I_NaP(v1,h1) + I_L(v1));
F2 = -(I_NaP(v2,h2) + I_L(v2));
F3 = -(I_NaP(v3,h3) + I_L(v3));

dv1dt = (F1 - G.I*(B.b21*S_inf(v2) + B.b31*S_inf(v3))*(v1-E.I) - G.E*d1*(v1-E.E))/C;
dv2dt = (F2 - G.I*(B.b12*S_inf(v1) + B.b32*S_inf(v3))*(v2-E.I) - G.E*d2*(v2-E.E))/C;
dv3dt = (F3 - G.I*(B.b13*S_inf(v1) + B.b23*S_inf(v2))*(v3-E.I) - G.E*d3*(v3-E.E))/C;
dh1dt = (h_inf(v1)-h1)*tauh(v1);
dh2dt = (h_inf(v2)-h2)*tauh(v2);
dh3dt = (h_inf(v3)-h3)*tauh(v3);

dPdt = [dv1dt; dv2dt; dv3dt; dh1dt; dh2dt; dh3dt];

end