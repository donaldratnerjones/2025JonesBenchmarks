function [y3_output,G3] =  BB_power(x_sh,x,y_input,x_sh_i0,x_loc_i0,Consts)

%-----THIS SECTION COMPUTES SFC, ESF, AND ENGINE WEIGHT-----% 
X3=x;
Z=x_sh;
C=Consts;
Thrust = y_input(6);
Dim_Throttle = X3(1)*16168;   %--non-diminsional throttle setting

%-----Surface fit to engine deck (obtained using least spuares approx)-----%

s=[1.13238425638512 1.53436586044561 -0.00003295564466 -0.00016378694115 -0.31623315541888 0.00000410691343 -0.00005248000590 -0.00000000008574 0.00000000190214 0.00000001059951];
y3_output(1)=s(1)+s(2)*Z(3)+s(3)*Z(2)+s(4)*Dim_Throttle+s(5)*Z(3)^2+2*Z(2)*Z(3)*s(6)+2*Dim_Throttle*Z(3)*s(7)+s(8)*Z(2)^2+2*Dim_Throttle*Z(2)*s(9)+s(10)*Dim_Throttle^2; %SFC

y3_output(3) = (Thrust/2)/Dim_Throttle; %ESF
y3_output(2) = C(4)*(y3_output(3)^1.05)*2;  %We

%-----THIS SECTION COMPUTES SFC, ESF, AND ENGINE WEIGHT-----% 


%-----THIS SECTION COMPUTES POLYNOMIAL CONSTRAINT FUNCTIONS-----%

G(1)=y3_output(3);   %--engine scale factor

S_initial1=[x_sh_i0(3),x_sh_i0(2),x_loc_i0(1)];
S1=[Z(3),Z(2),X3(1)];
flag1 = [2,4,2];
bound1 = [.25,.25,.25];
G(2) = PolyApprox(S_initial1,S1,flag1,bound1);   %--engine temperature

p=[11483.7822254806 10856.2163466548 -0.5080237941 3200.157926969 -0.1466251679 0.0000068572];  
Throttle_uA=p(1)+p(2)*Z(3)+p(3)*Z(2)+p(4)*Z(3)^2+2*p(5)*Z(3)*Z(2)+p(6)*Z(2)^2;
G(3)=Dim_Throttle/Throttle_uA-1;   %--throttle setting

%-----THIS SECTION COMPUTES POLYNOMIAL CONSTRAINT FUNCTIONS-----%

%-----Constrain function limits (for use in optimization)-----%

Temp_uA=1.02;


G3(1)=G(2)/Temp_uA-1;
G3(2)=G(3);
% G3(3)=y3_output(3)/1.5-1;  %output constraints on ESF 
% G3(4)=0.5-y3_output(3);
% G3(5)=y3_output(2)/30000-1;  %output constraints on We
% G3(6)=(100-y3_output(2))/10000;
        
