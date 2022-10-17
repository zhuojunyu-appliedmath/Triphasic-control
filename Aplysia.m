function P = analytic_Aplysia1(initials,t0,tF,dt,a1,a2,a3,rho)

tspan = t0:dt:tF; n = length(tspan);
P = zeros(n,3); P(1,:) = initials;
X = P(:,1); Y = P(:,2); Z = P(:,3);

% region 1
expA1 = @(t,tin) [exp(tin-t) rho*(exp(tin-t)-exp(t-tin))/2 0;
                  0 exp(t-tin) 0;
                  0 0 exp((1-rho)*(t-tin))];
intA1B1 = @(t,tin) [1+(a2-a1)*rho-a2*rho*(exp(t-tin)+exp(tin-t))/2+(a1*rho-1)*exp(tin-t);
                    -a2+a2*exp(t-tin);
                    a3-a3*exp((1-rho)*(t-tin))];
P1 = @(t,tin,xin) expA1(t,tin)*xin'+intA1B1(t,tin);
% region 2
expA2 = @(t,tin) [exp((1-rho)*(t-tin)) 0 0;
                  0 exp(tin-t) rho*(exp(tin-t)-exp(t-tin))/2;
                  0 0 exp(t-tin)];
intA2B2 = @(t,tin) [a1-a1*exp((1-rho)*(t-tin));
                    1+(a3-a2)*rho-a3*rho*(exp(t-tin)+exp(tin-t))/2+(a2*rho-1)*exp(tin-t);
                    -a3+a3*exp(t-tin)];
P2 = @(t,tin,xin) expA2(t,tin)*xin'+intA2B2(t,tin);
% region 3
expA3 = @(t,tin) [exp(t-tin) 0 0;
                  0 exp((1-rho)*(t-tin)) 0;
                  rho*(exp(tin-t)-exp(t-tin))/2 0 exp(tin-t)];
intA3B3 = @(t,tin) [-a1+a1*exp(t-tin);
                    a2-a2*exp((1-rho)*(t-tin));
                    1+(a1-a3)*rho-a1*rho*(exp(t-tin)+exp(tin-t))/2+(a3*rho-1)*exp(tin-t)];
P3 = @(t,tin,xin) expA3(t,tin)*xin'+intA3B3(t,tin);

i = 1; j = 1; k = 1; flag = 0;

if (X(1)>=Y(1)+(a1+a2)/2) && (X(1)>=Z(1)-(a1+a3)/2)
    flag = 1;
elseif (Y(1)>X(1)-(a1+a2)/2) && (Y(1)>=Z(1)+(a2+a3)/2)
    flag = 2;
elseif (Z(1)>X(1)+(a1+a3)/2) && (Z(1)>Y(1)-(a2+a3)/2)
    flag = 3;
end

if flag == 1
while ( i<n && j<n && k<n )

t_in1 = tspan(k); x_in1 = P(k,:);
for i = k:n-1
    if (X(i)>=Y(i)+(a1+a2)/2) && (X(i)>=Z(i)-(a1+a3)/2)
        i = i+1;
        P(i,:) = P1(tspan(i),t_in1,x_in1);
        X(i) = P(i,1); Y(i) = P(i,2); Z(i) = P(i,3);
    else
        break
    end
end
t_out1 = tspan(i-1); x_out1 = P(i-1,:);

if i <= n-1
t_in2 = tspan(i); x_in2 = P(i,:);
for j = i:n-1
    if (Y(j)>X(j)-(a1+a2)/2) && (Y(j)>=Z(j)+(a2+a3)/2)
        j = j+1;
        P(j,:) = P2(tspan(j),t_in2,x_in2);
        X(j) = P(j,1); Y(j) = P(j,2); Z(j) = P(j,3);
    else
        break
    end
end
t_out2 = tspan(j-1); x_out2 = P(j-1,:);
end

if j <= n-1
t_in3 = tspan(j); x_in3 = P(j,:);
for k = j:n-1
    if (Z(k)>X(k)+(a1+a3)/2) && (Z(k)>Y(k)-(a2+a3)/2)
        k = k+1;
        P(k,:) = P3(tspan(k),t_in3,x_in3);
        X(k) = P(k,1); Y(k) = P(k,2); Z(k) = P(k,3);
    else
        break
    end
end
t_out3 = tspan(k-1); x_out3 = P(k-1,:);
end

end
end

if flag == 2
while ( i<n && j<n && k<n )

t_in2 = tspan(i); x_in2 = P(i,:);
for j = i:n-1
    if (Y(j)>X(j)-(a1+a2)/2) && (Y(j)>=Z(j)+(a2+a3)/2)
        j = j+1;
        P(j,:) = P2(tspan(j),t_in2,x_in2);
        X(j) = P(j,1); Y(j) = P(j,2); Z(j) = P(j,3);
    else
        break
    end
end
t_out2 = tspan(j-1); x_out2 = P(j-1,:);

if j <= n-1
t_in3 = tspan(j); x_in3 = P(j,:);
for k = j:n-1
    if (Z(k)>X(k)+(a1+a3)/2) && (Z(k)>Y(k)-(a2+a3)/2)
        k = k+1;
        P(k,:) = P3(tspan(k),t_in3,x_in3);
        X(k) = P(k,1); Y(k) = P(k,2); Z(k) = P(k,3);
    else
        break
    end
end
t_out3 = tspan(k-1); x_out3 = P(k-1,:);
end

if k <= n-1
t_in1 = tspan(k); x_in1 = P(k,:);
for i = k:n-1
    if (X(i)>=Y(i)+(a1+a2)/2) && (X(i)>=Z(i)-(a1+a3)/2)
        i = i+1;
        P(i,:) = P1(tspan(i),t_in1,x_in1);
        X(i) = P(i,1); Y(i) = P(i,2); Z(i) = P(i,3);
    else
        break
    end
end
t_out1 = tspan(i-1); x_out1 = P(i-1,:);
end

end
end

if flag == 3
while ( i<n && j<n && k<n )
    
t_in3 = tspan(j); x_in3 = P(j,:);
for k = j:n-1
    if (Z(k)>X(k)+(a1+a3)/2) && (Z(k)>Y(k)-(a2+a3)/2)
        k = k+1;
        P(k,:) = P3(tspan(k),t_in3,x_in3);
        X(k) = P(k,1); Y(k) = P(k,2); Z(k) = P(k,3);
    else
        break
    end
end
t_out3 = tspan(k-1); x_out3 = P(k-1,:);  

if k <= n-1
t_in1 = tspan(k); x_in1 = P(k,:);
for i = k:n-1
    if (X(i)>=Y(i)+(a1+a2)/2) && (X(i)>=Z(i)-(a1+a3)/2)
        i = i+1;
        P(i,:) = P1(tspan(i),t_in1,x_in1);
        X(i) = P(i,1); Y(i) = P(i,2); Z(i) = P(i,3);
    else
        break
    end
end
t_out1 = tspan(i-1); x_out1 = P(i-1,:);
end

if i <= n-1
t_in2 = tspan(i); x_in2 = P(i,:);
for j = i:n-1
    if (Y(j)>X(j)-(a1+a2)/2) && (Y(j)>=Z(j)+(a2+a3)/2)
        j = j+1;
        P(j,:) = P2(tspan(j),t_in2,x_in2);
        X(j) = P(j,1); Y(j) = P(j,2); Z(j) = P(j,3);
    else
        break
    end
end
t_out2 = tspan(j-1); x_out2 = P(j-1,:);
end

end
end

end