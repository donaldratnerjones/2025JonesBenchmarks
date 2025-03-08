
function  [Yin,Yout] = my_system_analysis(Z,X,Yin,Z0,X0,Y0,Consts)        
    y1 = BB_weight(Z,X(:,1:19),Yin,Consts);
    y2 = BB_dragpolar(Z,X(:,20:22),Yin,Z0,Y0,Consts);
    y3 = BB_power(Z,X(:,23),Yin,Z0,X0,Consts);
    y4 = BB_range(Z,Yin);
    %            1     2     3     4     5     6     7     8     9     10
    Yout    = [y2(1) y3(2) y4(1) y1(3) y3(3) y2(2) y1(2) y2(3) y3(1) y1(1)];
end
