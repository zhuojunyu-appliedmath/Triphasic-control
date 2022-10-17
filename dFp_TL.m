function dFp = dFp_TL(x2,x3,epsilon,delta,theta1)

W13 = -1+epsilon; 
W12 = -1-delta; 

if W13*x3+W12*x2+theta1 > 0
    dFp1 = 1;
else
    dFp1 = 0;
end

dFp = [dFp1; 0; 0];

end