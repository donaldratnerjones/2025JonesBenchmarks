function[P,Mz,Mx,bend_twist,Spanel]=loads(b,c,Sweep_40,D_mx,L,Izz,Z,E)

np=9;   %----number of panels per halfspan
n=90;
h=(b/2)/n;
x=[0:(b/2-h)/(n-1):b/2-h];
x1=[h:(b/2-h)/(n-1):b/2];

%----Calculate Mx, Mz, and P----%

l=[0:(b/2)/np:(b/2)-(b/2)/np];

for i=0:np-1
  clear C C1
  C=c(4)+2*((b/2-x(i*(n/np)+1:(i+1)*(n/np)))/b)*(c(1)-c(4));
  C1=c(4)+2*((b/2-x1(i*(n/np)+1:(i+1)*(n/np)))/b)*(c(1)-c(4));
  f=(3*b/10)*sqrt(1-(x(i*(n/np)+1:(i+1)*(n/np)).^2)/((b/2)^2));
  f1=(3*b/10)*sqrt(1-(x1(i*(n/np)+1:(i+1)*(n/np)).^2)/((b/2)^2));
  A_Tot=(h/4)*(C+C1).*(f+f1);
  Area(i+1)=sum(A_Tot);
  Spanel(i+1)=(h*(n/np)/2)*(C(1)+C(10));
end

p=L*Area./sum(Area);
for i=1:np
  T(i)=sum(p(i:np));
  Mb(i)=sum(p(i+1:np).*(l(i+1:np)-l(i)))./(cos(Sweep_40*pi/180));
end

P=T(1:np/3:np);
P=P(:);
Mx=P.*D_mx;
Mz=Mb(1:np/3:np);
Mz=Mz(:);

%----Calculate Wing Twist due to Bending----%

chord=c(4)+(2*(b/2-l)./b)*(c(1)-c(4));
y(1,:)=(l-.4*chord.*(tan(Sweep_40*pi/180))*((cos(Sweep_40*pi/180))^2))/(cos(Sweep_40*pi/180));
y(2,:)=(l+.6*chord.*(tan(Sweep_40*pi/180))*((cos(Sweep_40*pi/180))^2))/(cos(Sweep_40*pi/180));
y(2,1)=0;
I(1:np/3)=sqrt((Izz(1)^2+Izz(2)^2)/2);
I(np/3+1:2*np/3)=sqrt((Izz(2)^2+Izz(3)^2)/2);
I(2*np/3+1:np)=sqrt((Izz(3)^2)/2);

La=y(1,2:np)-y(1,1:np-1);
La=[0 La];
Lb=y(2,2:np)-y(2,1:np-1);
Lb=[0 Lb];
A=T.*La.^3./(3*E*I)+Mb.*La.^2./(2*E*I);
B=T.*Lb.^3./(3*E*I)+Mb.*Lb.^2./(2*E*I);
Slope_A=T.*La.^2./(2*E*I)+Mb.*La./(E*I);
Slope_B=T.*Lb.^2./(2*E*I)+Mb.*Lb./(E*I);
for i=1:np-1
  Slope_A(i+1)=Slope_A(i)+Slope_A(i+1);
  Slope_B(i+1)=Slope_B(i)+Slope_B(i+1);
  A(i+1)=A(i)+Slope_A(i)*La(i+1)+A(i+1);
  B(i+1)=B(i)+Slope_B(i)*Lb(i+1)+B(i+1);
end
bend_twist=((B-A)./chord)*180/pi;
for i=2:size(bend_twist,2)
   if bend_twist(i)<bend_twist(i-1)
      bend_twist(i)=bend_twist(i-1);
   end
end
%delta_L=bend_twist.*Spanel*q*0.1;
%p_twist=p+delta_L;
%r=1:45/30:45;
%plot(A)
%hold on
%plot(B)
