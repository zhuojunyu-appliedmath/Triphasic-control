function DF = Jacobian_TL(x1,x2,x3,epsilon,delta,theta1,theta2,theta3)

W13 = -1+epsilon; W21 = -1+epsilon; W32 = -1+epsilon; 
W12 = -1-delta; W23 = -1-delta; W31 = -1-delta; 

if W13*x3+W12*x2+theta1 > 0
    A12 = W12; A13 = W13;
else
    A12 = 0; A13 = 0;
end

if W21*x1+W23*x3+theta2 > 0
    A21 = W21; A23 = W23;
else
    A21 = 0; A23 = 0;
end

if W32*x2+W31*x1+theta3 > 0
    A31 = W31; A32 = W32;
else
    A31 = 0; A32 = 0;
end

DF = [-1 A12 A13;...
      A21 -1 A23;...
      A31 A32 -1];

end