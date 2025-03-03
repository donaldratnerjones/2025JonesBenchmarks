function[y4_output,G4]=BB_range(x_sh,y_input)

%-----THIS SECTION COMPUTES THE A/C RANGE (Breguet)-----%
Z=x_sh;

if Z(2)<36089
     theta=1-0.000006875*Z(2);
else
     theta=.7519;
end

Wt = y_input(2)+y_input(7)+y_input(10);
y4_output(1) = Wt;

%-----THIS SECTION COMPUTES THE A/C RANGE-----%
range=((Z(3)*y_input(8))*661*sqrt(theta)/y_input(9))*log(Wt/(Wt-y_input(7)));

G4= -range/2000+1;