function [X] = simulate_4state(p,x_0,u,v)

x_current = x_0;
X = zeros(4,length(u));

for k=1:length(u)

X(:,k) = x_current;
x_next(1) = x_current(1)*(1-0.001*(u(k)+p(5))) + x_current(3)*0.001*p(4) + x_current(2)*0.001*p(1);
x_next(2) = x_current(2)*(1-0.001*(p(1)+p(2))) + x_current(1)*0.001*u(k) + x_current(3)*0.001*p(3);
x_next(3) = x_current(3)*(1-0.001*(p(3)+p(4)+p(6))) + x_current(2)*0.001*p(2) + x_current(1)*0.001*p(5) + x_current(4)*0.001*v(k);
x_next(4) = x_current(4)*(1-0.001*v(k)) + x_current(3)*0.001*p(6);
x_current = x_next;

end

end

