
function [automotive,frel,grel,hrel,xl,xu,xopt,x0] = ...
          automotive_benchmark()

    load('automotive.mat','e2x','frel','glimit','grel','gsense',...
                          'gsurf','hrel','massmult','mat','r2x','s2x',...
                          'u2x','x0','xl','xopt','xu','lb','ub');

    automotive = @(x) mopta(x,massmult,gsurf,u2x,e2x,s2x,r2x,...
                              mat,gsense,glimit,lb,ub);    
     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f,g,h] = mopta(x,massmult,gsurf,u2x,e2x,s2x,r2x,...
                         mat,gsense,glimit,lb,ub)
    x = lb + x .* (ub-lb);
    f = mopta_objective(x,massmult);
    g = mopta_constraints(x,gsurf,u2x,e2x,s2x,r2x,mat);
    g = ( gsense .* ( g - glimit ) ./ abs(glimit) )';
    h = [];
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f,gradf] = mopta_objective(x,massmult)
    f = x * massmult;
    gradf = massmult;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function g = mopta_constraints(x,gsurf,u2x,e2x,s2x,r2x,mat)

    g = zeros(68,1);
    
    % USNCAP (g1)

    d = gsurf{1}.d;
    xus = zeros(1,d);
    for k=1:d
        i = u2x(k);
        if i > 0 
            xus(k) = x(i);
        else
            xus(k) = mat(-i);
        end
    end
    j = 1;
    s = gsurf{j}.s;
    uselog = gsurf{j}.uselog;
    offset = gsurf{j}.offset;
    beta = gsurf{j}.beta;
    theta = gsurf{j}.theta;
    pcoef = gsurf{j}.pcoef;
    g(j) = predict(xus,s,uselog,offset,beta,theta,pcoef);

    % EURONCAP (g2-g6)

    d = gsurf{2}.d;
    xeuro = zeros(1,d);
    for k=1:d
        i = e2x(k);
        if i > 0 
            xeuro(k) = x(i);
        else
            xeuro(k) = mat(-i);
        end
    end
    for j=2:6
        s = gsurf{j}.s;
        uselog = gsurf{j}.uselog;
        offset = gsurf{j}.offset;
        beta = gsurf{j}.beta;
        theta = gsurf{j}.theta;
        pcoef = gsurf{j}.pcoef;
        g(j) = predict(xeuro,s,uselog,offset,beta,theta,pcoef);
    end

    % SIDE  (g7-g9)

    d = gsurf{7}.d;
    xside = zeros(1,d);
    for k=1:d
        i = s2x(k);
        if i > 0 
            xside(k) = x(i);
        else
            xside(k) = mat(-i);
        end
    end
    for j=7:9
        s = gsurf{j}.s;
        uselog = gsurf{j}.uselog;
        offset = gsurf{j}.offset;
        beta = gsurf{j}.beta;
        theta = gsurf{j}.theta;
        pcoef = gsurf{j}.pcoef;
        g(j) = predict(xside,s,uselog,offset,beta,theta,pcoef);
    end

    % NVH (g10-g21)

    for j=10:21
        s = gsurf{j}.s;
        uselog = gsurf{j}.uselog;
        offset = gsurf{j}.offset;
        beta = gsurf{j}.beta;
        theta = gsurf{j}.theta;
        pcoef = gsurf{j}.pcoef;
        g(j) = predict(x,s,uselog,offset,beta,theta,pcoef);
    end

    % DURABILITY (g22-g55)

    for j=22:55
        s = gsurf{j}.s;
        uselog = gsurf{j}.uselog;
        offset = gsurf{j}.offset;
        beta = gsurf{j}.beta;
        theta = gsurf{j}.theta;
        pcoef = gsurf{j}.pcoef;
        g(j) = predict(x,s,uselog,offset,beta,theta,pcoef);
    end

    % REAR  (g56-g68)

    d = gsurf{56}.d;
    xrear = zeros(1,d);
    for k=1:d
        i = r2x(k);
        if i > 0 
            xrear(k) = x(i);
        else
            xrear(k) = mat(-i);
        end
    end
    for j=56:68
        s = gsurf{j}.s;
        uselog = gsurf{j}.uselog;
        offset = gsurf{j}.offset;
        beta = gsurf{j}.beta;
        theta = gsurf{j}.theta;
        pcoef = gsurf{j}.pcoef;
        g(j) = predict(xrear,s,uselog,offset,beta,theta,pcoef);
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = predict(x,s,uselog,offset,beta,theta,pcoef)
       
    n = size(s,1);
       
    one = ones(n,1);
    
    corr = exp(- abs(one*x-s).^1.9 * theta); 
    
    pcorr = pcoef .* corr'; 
    
    y = [1,x]*beta + sum(pcorr);
    
    if uselog == 1
        expy = exp(y);
        y = expy + offset;
    end            
    
end


% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % function [f,gradf] = mopta_objective(x,massmult)
% %     f = x * massmult;
% %     gradf = massmult;
% % end
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % function [g,dgdx] = mopta_constraints(x,gsurf,u2x,e2x,s2x,r2x,mat)
% % 
% %     g = zeros(68,1);
% %     dgdx = zeros(124,68);
% %     
% %     % USNCAP (g1)
% % 
% %     d = gsurf{1}.d;
% %     xus = zeros(1,d);
% %     for k=1:d
% %         i = u2x(k);
% %         if i > 0 
% %             xus(k) = x(i);
% %         else
% %             xus(k) = mat(-i);
% %         end
% %     end
% %     j = 1;
% %     s = gsurf{j}.s;
% %     uselog = gsurf{j}.uselog;
% %     offset = gsurf{j}.offset;
% %     beta = gsurf{j}.beta;
% %     theta = gsurf{j}.theta;
% %     pcoef = gsurf{j}.pcoef;
% %     [g(j),dgj] = predict(xus,s,uselog,offset,beta,theta,pcoef);
% %     for k=1:d
% %         i = u2x(k);
% %         if i > 0 
% %             dgdx(i,j) = dgj(k);
% %         end
% %     end    
% % 
% %     % EURONCAP (g2-g6)
% % 
% %     d = gsurf{2}.d;
% %     xeuro = zeros(1,d);
% %     for k=1:d
% %         i = e2x(k);
% %         if i > 0 
% %             xeuro(k) = x(i);
% %         else
% %             xeuro(k) = mat(-i);
% %         end
% %     end
% %     for j=2:6
% %         s = gsurf{j}.s;
% %         uselog = gsurf{j}.uselog;
% %         offset = gsurf{j}.offset;
% %         beta = gsurf{j}.beta;
% %         theta = gsurf{j}.theta;
% %         pcoef = gsurf{j}.pcoef;
% %         [g(j),dgj] = predict(xeuro,s,uselog,offset,beta,theta,pcoef);
% %         for k=1:d
% %             i = e2x(k);
% %             if i > 0 
% %                 dgdx(i,j) = dgj(k);
% %             end
% %         end    
% %     end
% % 
% %     % SIDE  (g7-g9)
% % 
% %     d = gsurf{7}.d;
% %     xside = zeros(1,d);
% %     for k=1:d
% %         i = s2x(k);
% %         if i > 0 
% %             xside(k) = x(i);
% %         else
% %             xside(k) = mat(-i);
% %         end
% %     end
% %     for j=7:9
% %         s = gsurf{j}.s;
% %         uselog = gsurf{j}.uselog;
% %         offset = gsurf{j}.offset;
% %         beta = gsurf{j}.beta;
% %         theta = gsurf{j}.theta;
% %         pcoef = gsurf{j}.pcoef;
% %         [g(j),dgj] = predict(xside,s,uselog,offset,beta,theta,pcoef);
% %         for k=1:d
% %             i = s2x(k);
% %             if i > 0 
% %                 dgdx(i,j) = dgj(k);
% %             end
% %         end    
% %     end
% % 
% %     % NVH (g10-g21)
% % 
% %     for j=10:21
% %         s = gsurf{j}.s;
% %         uselog = gsurf{j}.uselog;
% %         offset = gsurf{j}.offset;
% %         beta = gsurf{j}.beta;
% %         theta = gsurf{j}.theta;
% %         pcoef = gsurf{j}.pcoef;
% %         [g(j), dgdx(:,j)] = predict(x,s,uselog,offset,beta,theta,pcoef);
% %     end
% % 
% %     % DURABILITY (g22-g55)
% % 
% %     for j=22:55
% %         s = gsurf{j}.s;
% %         uselog = gsurf{j}.uselog;
% %         offset = gsurf{j}.offset;
% %         beta = gsurf{j}.beta;
% %         theta = gsurf{j}.theta;
% %         pcoef = gsurf{j}.pcoef;
% %         [g(j),dgdx(:,j)] = predict(x,s,uselog,offset,beta,theta,pcoef);
% %     end
% % 
% %     % REAR  (g56-g68)
% % 
% %     d = gsurf{56}.d;
% %     xrear = zeros(1,d);
% %     for k=1:d
% %         i = r2x(k);
% %         if i > 0 
% %             xrear(k) = x(i);
% %         else
% %             xrear(k) = mat(-i);
% %         end
% %     end
% %     for j=56:68
% %         s = gsurf{j}.s;
% %         uselog = gsurf{j}.uselog;
% %         offset = gsurf{j}.offset;
% %         beta = gsurf{j}.beta;
% %         theta = gsurf{j}.theta;
% %         pcoef = gsurf{j}.pcoef;
% %         [g(j),dgj] = predict(xrear,s,uselog,offset,beta,theta,pcoef);
% %         for k=1:d
% %             i = r2x(k);
% %             if i > 0 
% %                 dgdx(i,j) = dgj(k);
% %             end
% %         end    
% %     end
% % 
% % end
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % function [y, dydx] = predict(x,s,uselog,offset,beta,theta,pcoef)
% %        
% %     [n,d] = size(s);
% %     
% %     dydx = beta(2:1+d,1);
% %     
% %     one = ones(n,1);
% %     
% %     corr =  exp(- abs(one*x-s).^1.9 * theta); % column vector of corr
% %     
% %     pcorr = pcoef .* corr'; % row vector of pcoef(i)corr(x,si)
% %     
% %     y = [1,x]*beta + sum(pcorr);
% %     
% %     for k=1:d
% %         for i=1:n
% %             if x(k) >= s(i,k)
% %                 dydx(k) = dydx(k) - pcorr(i)*theta(k)*1.9*(x(k)-s(i,k))^0.9;
% %             else
% %                 dydx(k) = dydx(k) + pcorr(i)*theta(k)*1.9*(s(i,k)-x(k))^0.9;
% %             end
% %         end
% %     end
% %     
% %     if uselog == 1
% %         expy = exp(y);
% %         y = expy + offset;
% %         dydx = expy*dydx;
% %     end            
% %     
% % end
