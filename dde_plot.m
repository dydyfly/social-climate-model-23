global tv tc v0 x0
tv=15.5
tc=0.3
v0=0.1
x0=0.1

fig=figure();
set(fig, 'Position', [100 100 400 300]);

sol = dde23(@climate,10,@initial,[10,200]);
plot(sol.x, sol.y,'LineWidth',2);
ylim([0.95,1.05]);
% ylim([-0.2, 1.2]);
grid on;
legend('Vegetation','Mitigation');
xlabel('time');
title("x_0="+num2str(x0)+",v_0="+num2str(v0)+",T_v="+num2str(tv)+char(176)+" C,T_c="+num2str(tc)+char(176)+" C");
     

function dvdt=forest(t,v) % forest die-back model
global tv
dvdt=2*(1-0.01*(tv-23-5*v)^2)*v*(1-v)-0.2*v;
end

function yp=climate(t,y,Z) % coupled social-climate model with time-delay
global tv tc
yp=[2*(1-0.01*(tv-23-5*y(1))^2)*(0.2+0.4*y(2))*y(1)*(1-y(1))-0.2*y(1)
    y(2)*(1-y(2))*(2*y(2)-2+5/(1+exp(3*tc+22.5*y(1)-22.5*Z(1))))];
end

function y=initial(t) % compute the initial condition with forest die-back model
global x0 v0
if t==0
    y=[x0
        v0];
else
    [no,noo]=ode23(@forest,[0,t],v0);
    len=length(noo);
    y=[noo(len) x0];
end
end