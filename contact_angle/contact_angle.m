function contact
%球心在上方为正
data=importdata('e:/brim_2.txt');
hs=size(data,1);
f=@(p,data)(data(:,1)-p(1)).^2+(data(:,2)-p(2)).^2+(data(:,3)-p(3)).^2-p(4)^2;
x0=sum(data(:,1))./hs;
y0=sum(data(:,2))./hs;
z0=1.88;
r=sqrt((max(data(:,1))-min(data(:,1))).^2);
p=nlinfit(data,zeros(hs,1),f,[x0 y0 z0 r]')
t=p(1)-abs(p(4)):0.005:p(1)+abs(p(4));
z1=sqrt(p(4)^2-(t-p(1)).^2)+p(3);
z2=-sqrt(p(4)^2-(t-p(1)).^2)+p(3);
plot(t,z1,'b','LineWidth',1.5)
hold on
plot(t,z2,'b','LineWidth',1.5)
hold on
wc=0.5;
j=0;
for i=1:hs
    if data(i,2)<(p(2)+wc) && data(i,2)>(p(2)-wc)
        j=j+1;
        q(j,1)=data(i,1);
        q(j,2)=data(i,3);
    end
end
plot(q(:,1),q(:,2),'*')
hold on
x=p(1)-2*abs(p(4)):0.01:p(1)+2*abs(p(4));
y=z0;
plot(x,y,'k','LineWidth',4)
hold on
plot(p(1),p(3),'.')
%[x,y,z]=ellipsoid(p(1),p(2),p(3),p(4),p(4),p(4),50);
%hold on
%plot3(data(:,1),data(:,2),data(:,3),'o')
%colormap(bone);
%surf(x,y,z) 
axis equal
h=p(3)-1.88;
the=acos(h./p(4))./pi*180
end

