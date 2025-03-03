
function [bjet,frelevant,grelevant,hrelevant,xl,xu,xopt] = ...
             businessjet_benchmark(UseEqualities,ScaleConstraints)

    load('businessjet.mat','Z0','Z_lb','Z_ub','X0','X_lb','X_ub','Y0','Y_lb','Y_ub',...
                           'const','cg','ch','frelevant','grelevant','hrelevant','xopt');

    if UseEqualities && ~ScaleConstraints
        bjet      = @(x) compute_IDF(x,Z0,X0,Y0,const,[],[]);
        xl        = [Z_lb,X_lb,Y_lb];
        xu        = [Z_ub,X_ub,Y_ub];
    elseif UseEqualities && ScaleConstraints
        bjet      = @(x) compute_IDF(x,Z0,X0,Y0,const,cg,ch);
        xl        = [Z_lb,X_lb,Y_lb];
        xu        = [Z_ub,X_ub,Y_ub];
    elseif ~UseEqualities && ~ScaleConstraints
        bjet      = @(x) compute_MDF(x,Z0,X0,Y0,const,[],xopt);
        frelevant = true(1,31);
        grelevant = true(78,31);
        hrelevant = [];
        xl        = [Z_lb,X_lb];
        xu        = [Z_ub,X_ub];
        xopt      = xopt(1:31);
    elseif ~UseEqualities && ScaleConstraints
        bjet      = @(x) compute_MDF(x,Z0,X0,Y0,const,cg,xopt);
        frelevant = true(1,31);
        grelevant = true(78,31);
        hrelevant = [];
        xl        = [Z_lb,X_lb];
        xu        = [Z_ub,X_ub];
        xopt      = xopt(1:31);
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [f,g,h] = compute_IDF(x,Z0,X0,Y0,const,cg,ch)    
    
    Z   = x(1:8);
    X   = x(9:31);
    Yin = x(32:41);

    f = BB_range(Z,Yin)/1000; 

    [y1,G1] = BB_weight(Z,X(:,1:19),Yin,const);
    [y2,G2] = BB_dragpolar(Z,X(:,20:22),Yin,Z0,Y0,const);
    [y3,G3] = BB_power(Z,X(:,23),Yin,Z0,X0,const);
    [y4,G4] = BB_range(Z,Yin);
    
    g = real([G1 G2 G3 G4]);

    if ~isempty(cg)
        g = 2./(1+exp(-g./cg)) - 1 ;
    end

    Yout = [y2(1) y3(2) y4(1) y1(3) y3(3) y2(2) y1(2) y2(3) y3(1) y1(1)];

    h = real( (Yout-Yin)./Y0 );

    if ~isempty(ch)
        h = 2./(1+exp(-h./ch)) - 1 ;
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [f,g,h] = compute_MDF(x,Z0,X0,Y0,const,cg,xopt)    
    
    Z   = x(1:8);
    X   = x(9:31);
    [Yin,Yerror] = converge(Z,X,Z0,X0,Y0,const,xopt);
    if Yerror > 1e-6
        f = nan;
        g = nan(1,78);
        h = [];
        return
    end
    
    f   = BB_range(Z,Yin)/1000; 

    [~,G1] = BB_weight(Z,X(:,1:19),Yin,const);
    [~,G2] = BB_dragpolar(Z,X(:,20:22),Yin,Z0,Y0,const);
    [~,G3] = BB_power(Z,X(:,23),Yin,Z0,X0,const);
    [~,G4] = BB_range(Z,Yin);
    
    g = real([G1 G2 G3 G4]);

    if ~isempty(cg)
        g = 2./(1+exp(-g./cg)) - 1 ;
    end

    h = [];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Yerror = yobj(Z,X,Yin,Z0,X0,Y0,const)

    [Yin,Yout] = my_system_analysis(Z,X,Yin,Z0,X0,Y0,const);
    Yerror = norm( real( (Yout-Yin)./Y0 ) ) ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [C,Ceq] = ynonlcon(Z,X,Yin,Z0,X0,Y0,const)
    [Yin,Yout] = my_system_analysis(Z,X,Yin,Z0,X0,Y0,const);
    C = [];
    Ceq = real( (Yout-Yin)./Y0 ) ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Yconverged,Yerror] = converge(Z,X,Z0,X0,Y0,const,xopt)
    Yout = xopt(32:41);
    n    = 100;
    Yin0 = Yout;
    Yerror0 = inf;
    for i=1:n
        Yin = Yout;
        [Yin,Yout] = my_system_analysis(Z,X,Yin,Z0,X0,Y0,const);
        Yerror = norm( real( (Yout-Yin)./Y0 ) ) ;
        % fprintf('iter = %d, norm=%e\n', i, Yerror );
        if ~any(isnan(Yin)) && ~any(isnan(Yout))
            if Yerror < Yerror0
                Yin0 = Yin;
                Yerror0 = Yerror;
            end
            if Yerror0 == 0
                break
            end                    
        else
            break
        end
    end
    Yin = Yin0;
    Yerror = Yerror0;
    if Yerror > 1e-6
        fun     = @(Yin) yobj(Z,X,Yin,Z0,X0,Y0,const);
        nonlcon = @(Yin) ynonlcon(Z,X,Yin,Z0,X0,Y0,const);
        options = optimoptions('fmincon',Display='none',Algorithm='sqp',...
                               ConstraintTolerance=2e-16,StepTolerance=1e-16,...
                               FiniteDifferenceType='central',MaxFunctionEvaluations=1000);
        Yin = xopt(32:41);
        Yin = fmincon(fun,Yin,[],[],[],[],[],[],nonlcon,options);
    end
    [Yin,Yout] = my_system_analysis(Z,X,Yin,Z0,X0,Y0,const);
    Yconverged = real(Yin+Yout)/2;
    Yerror = norm( real( (Yout-Yin)./Y0 ) ) ;
end
