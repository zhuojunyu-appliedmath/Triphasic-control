function dPdt = threshold_linear(t,P,epsilon,delta,theta1,theta2,theta3)

x1 = P(1); x2 = P(2); x3 = P(3);

W13 = -1+epsilon; W21 = -1+epsilon; W32 = -1+epsilon; 
W12 = -1-delta; W23 = -1-delta; W31 = -1-delta; 

dx1dt = -x1+max(W13*x3+W12*x2+theta1,0);
dx2dt = -x2+max(W21*x1+W23*x3+theta2,0);
dx3dt = -x3+max(W32*x2+W31*x1+theta3,0);

dPdt = [dx1dt; dx2dt; dx3dt];

end