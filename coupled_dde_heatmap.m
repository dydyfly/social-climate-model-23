global tv tc v0 x0
v0=0.1
x0=0.9

A=zeros(29,51);

for i=10:0.5:35
    tv=i;
    for j=0.2:0.1:3
        tc=j;
        r=int8(10*j-1);
        c=int8(2*i-19);
        sol = dde23(@climate,10,@initial,[10,200]);
        vegetation=sol.y(1,:);
        mitigation=sol.y(2,:);
        len=length(vegetation);
        sample=len-int32(len/4);
        maxv=max(vegetation(sample:len));
        minv=min(vegetation(sample:len));
        maxx=max(mitigation(sample:len));
        minx=min(mitigation(sample:len));
        varv=maxv-minv;
        varx=maxx-minx;
        meanx=(maxx+minx)/2;
        meanv=(maxv+minv)/2;
        if varv>0.1 | varx>0.1
            A(r,c)=3;
        else
            A(r,c)=meanx+meanv;
        end
    end
end
x=[10 35];
y=[0.2 3];
imagesc(x,y,A)
set(gcf,'Position',[100 100 240 200])
title("x_0="+num2str(x0)+",v_0="+num2str(v0));
xlabel('T_v');
% ylabel('T_c');
% colormap(hsv(512))
% colorbar('Ticks',[0,1,2,3],'TickLabels',{'Barren','Brown', ...
%     'Green','Swing'})

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