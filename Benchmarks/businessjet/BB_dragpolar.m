function[y2_output,G2]=BB_dragpolar(x_sh,x,y_input,x_sh_i0,y_input_i0,Consts)

%----Inputs----%
Z=x_sh;
X2=x;
C=Consts;
Wt=y_input(3);
twst=y_input(4);
ESF=y_input(5);
ARht=Z(8);
sweepht=X2(1);
Lw=X2(2);  
Lh=X2(3);    
S_ht=Z(7);     
Nh=C(6);
       
%-----Drag computations----%

if Z(2)<36089
     V = Z(3)*(1116.39*sqrt(1-(6.875e-06*Z(2))));
     rho = (2.377e-03)*(1-(6.875e-06*Z(2)))^4.2561;
else
     V = Z(3)*968.1;
     rho = (2.377e-03)*(.2971)*exp(-(Z(2)-36089)/20806.7);
end
q=.5*rho*(V^2);
%%%% Modified by S. Tosserams:
%%%% scale coefficients for proper conditioning of matrix A 
a=q*Z(6)/1e5;
b=Nh*q*S_ht/1e5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c=Lw;
d=(Lh)*Nh*(S_ht/Z(6));
A=[a b; c d];
%%%% Modified by S. Tosserams:
%%%% scale coefficient Wt for proper conditioning of matrix A 
B=[Wt/1e5; 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CLo=(pinv(A)*B);
delta_L=twst*q;
Lw1=CLo(1)*q*Z(6)-delta_L;
CLw1=Lw1/(q*Z(6));
CLht1=-CLw1*c/d;
%%%% Modified by S. Tosserams:
%%%% scale first coefficient of D for proper conditioning of matrix A 
D=[(Wt-CLw1*a-CLht1*b)/1e5; -CLw1*c-CLht1*d];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DCL=pinv(A)*D; % changed by DRJ
if Z(3)>= 1
     kw=Z(4)*(Z(3)^2-1)*cos(Z(5)*pi/180)/(4*Z(4)*sqrt(Z(3)^2-1)-2);
     kht=ARht*(Z(3)^2-1)*cos(sweepht*pi/180)/(4*ARht*sqrt(Z(3)^2-1)-2);
else
     kw=1/(pi*0.8*Z(4));
     kht=1/(pi*0.8*ARht);
end  

%-----Polynomial function modifying CDmin for ESF-----%

S_initial1=[y_input_i0(5)];
S1=ESF;
flag1 = [1];
bound1 = [.25];
Fo1 = PolyApprox(S_initial1,S1,flag1,bound1); 
CDmin = C(5)*Fo1 + 3.05*(Z(1)^(5/3))*((cos(Z(5)*pi/180))^(3/2));

CDw=CDmin+kw*(CLo(1)^2)+kw*(DCL(1)^2);
CDht=kht*(CLo(2)^2)+kht*(DCL(2)^2);
CD=CDw+CDht;
CL=CLo(1)+CLo(2);
y2_output(1) = Wt;
y2_output(2) = q*CDw*Z(6)+q*CDht*Z(7);
y2_output(3) = CL/CD;


%-----THIS SECTION COMPUTES CONSTRAINT POLYNOMIALS-----%

S_initial2=[x_sh_i0(1)];
S2=[Z(1)];
flag2=[1];
bound2=[.25];
G(1)=PolyApprox(S_initial2,S2,flag2,bound2);   %--adverse pressure gradient

%------Constraints------%

Pg_uA=1.1;
G2(1)=G(1)/Pg_uA-1;
if CLo(1) > 0
    G2(2)=(2*(CLo(2)))-(CLo(1));
    G2(3)=(2*(-CLo(2)))-(CLo(1));
else
    G2(2)=(2*(-CLo(2)))-(CLo(1));
    G2(3)=(2*(CLo(2)))-(CLo(1));
end
% G2(4) = y2_output(1)/1e5-1;  %output constraints on L
% G2(5) = -y2_output(1)/5e3+1;
% G2(6) = y2_output(2)/7e4-1;  %output constraints on D
% G2(7) = -y2_output(2)/100+1;
% G2(8) = y2_output(3)/10-1;  %output constraints on L/D
% G2(9) = -y2_output(3)/.1+1;






