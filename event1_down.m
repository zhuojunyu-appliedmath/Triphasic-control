function [value,isterminal,direction] = event1_down(t,P,d1,d2,d3,theta_I)

value=P(1)-theta_I;
isterminal=0;
direction=-1;

end