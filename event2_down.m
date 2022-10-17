function [value,isterminal,direction] = event2_down(t,P,d1,d2,d3,theta_I)

value=P(2)-theta_I;
isterminal=1;
direction=-1;

end