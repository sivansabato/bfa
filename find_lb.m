function [final_bestobj, info] = find_lb(wg, pg_fact, pi_g_fact, beta, tolerance, divisionfactor)
%Inputs: 
%wg - column vector of population weight in each region (sums to 1)
%pg_fact - column vector of p_g^1*divisionfactor
%pi_g_fact - column vector of pi_g^1*divisionfactor
%beta - value of beta (scalar)
%tolerance - gamma parameter of the algorithm
%divisionfactor - a number that can be used to provide pg and pi_g in
%different units. Used for numerical purposes. Default value is 1.
%Outputs:
%final_bestobj: value of lower bound
%info: additional information on the solution

fprintf('running with beta = %f, tolerance = %f\n', beta, tolerance);
if (~exist('divisionfactor','var'))
    divisionfactor = 1;
end
numg = length(wg);
pg = pg_fact/divisionfactor;
pi_g = pi_g_fact/divisionfactor;
MAXVALS = 5*10^5; %number of values to handle in a single batch

rg = 1 - pg_fact./pi_g_fact;
qg = 1./pi_g - 1;

mygplus = pi_g > 0 & pi_g < 1;
pg_piisone{1+0} = pg(pi_g == 0);
pg_piisone{1+1} = pg(pi_g == 1);

aupperlimit(1+0) = min(1,max(pg(pi_g ~= 1)./(1-pi_g(pi_g ~= 1)))); %(1-rg)/qg
aupperlimit(1+1) = min(1,max((1-pg(pi_g ~= 0))./pi_g(pi_g ~= 0)));


v{1+0} = [0;1;myround01(-rg(mygplus)./qg(mygplus)); myround01((1-rg(mygplus))./qg(mygplus)); pg_piisone{1+0}];
v{1+1} = [0;1;myround01(rg(mygplus)); myround01(rg(mygplus)+qg(mygplus)); (1-pg_piisone{1+1})];
for i=[0,1]
    v{1+i} = unique(v{1+i});
    v{1+i} = v{1+i}((v{1+i} >= 0) & (v{1+i} <=aupperlimit(i+1)));
end

[v0all,v1all] = meshgrid(v{1+0}, v{1+1});
allpairs = [v0all(:), v1all(:)];


besta0 = nan;
besta1 = nan;
bestobj = inf;

obj_all = calcO2_all(allpairs);
[bestobj, bestloc] = min(obj_all);
besta0 = allpairs(bestloc,1);
besta1 = allpairs(bestloc,2);



for z = find(mygplus')
    fprintf('z = %d. ', z);
    coefs = nan(numg,4,6);
    for g = find((pi_g' == 0)  | (pi_g' == 1))
        if (pi_g(1) == 1)
            coefs(g,1,:) = [beta, 0, rg(z), qg(z), 1-pg(g),0];
        else
            coefs(g,1,:) = [beta, 0, 0, 1, pg(g),0];
        end
        coefs(g,2,:) = [0,0, inf, 0, inf, 0];
        coefs(g,3,:) = [0, 0, inf ,0, inf, 0];
        coefs(g,4,:) = [0, 0, inf, 0, inf, 0];
    end
    for g = find(mygplus') %#ok<*FXUP>
        if (rg(g) > 0)
            coefs(g,1,:) = [beta*pi_g(g),0,rg(z), qg(z),rg(g), 0];
        else
            coefs(g,1,:) = [beta*(1-pi_g(g)),0,0,1,myround01(-rg(g)/qg(g)),0];
        end
        if (rg(g) + qg(g) < 1)
            coefs(g,2,:) = [beta*(1-pi_g(g)),0,rg(z),qg(z), myround01(rg(g) + qg(g)),0];
        else
            coefs(g,2,:) = [beta*pi_g(g),0,0,1,myround01((1-rg(g))/qg(g)),0];
        end
        coefs(g,3,:) = [beta*pi_g(g), 2*(1-beta)*(1-pi_g(g)),rg(z),qg(z),rg(g), qg(g)];
        coefs(g,4,:) = [beta*(1-pi_g(g)),2*(1-beta)*pi_g(g)*qg(z),0,1,myround01((rg(z)-rg(g))/qg(g)),myround01(qg(z)/qg(g))];
    end
    
    
    Acoef = 1;
    Bcoef = 2;
    acoef = 3;
    bcoef = 4;
    ccoef = 5;
    dcoef = 6;
    
    minlim = nan(numg,4);
    maxlim = nan(numg,4);
    midpoint = inf(numg,4);
    for g= 1:numg
        for i = 1:4
            if (coefs(g,1,dcoef) == 0)
                minlim(g,i) = -coefs(g,i,acoef)/coefs(g,i,bcoef);
                maxlim(g,i) = (1-coefs(g,i,acoef))/coefs(g,i,bcoef);
            else
                minlim(g,i) = max(-coefs(g,i,acoef)/coefs(g,i,bcoef),-coefs(g,i,ccoef)/coefs(g,i,dcoef));
                maxlim(g,i) = min((1-coefs(g,i,acoef))/coefs(g,i,bcoef),(1-coefs(g,i,ccoef))/coefs(g,i,dcoef));
            end
            if (myround01(coefs(g,1,dcoef) - coefs(g,1,bcoef)) ~= 0) %d \neq b
                midpoint(g,i) = (coefs(g,1,ccoef)-coefs(g,1,acoef))/(coefs(g,1,bcoef)-coefs(g,1,dcoef));
            end
        end
    end
    mino3 = max(0,-rg(z)/qg(z));
    maxo3 = min(aupperlimit(1+0),(aupperlimit(1+1)-rg(z))/qg(z));
    Theta = [mino3,maxo3,minlim(:)', maxlim(:)', midpoint(:)'];
    Theta = unique(Theta);
    Theta = Theta(Theta >= mino3 & Theta <= maxo3);
    maxdf_g = nan(1,numg);
    for j=1:length(Theta)-1
        for g=1:numg
            if (j > 1)
                oldrefpoint = refpoint;
                refpoint = (Theta(j) + Theta(j+1))/2;
                if (all((minlim(g,:) <= Theta(j) & maxlim(g,:) >= Theta(j+1)) ==  (minlim(g,:) <= Theta(j-1) & maxlim(g,:) >= Theta(j))) && ...
                        all(sign(refpoint - midpoint(g,:)) == sign(oldrefpoint - midpoint(g,:))))
                    continue;
                end
            else
                refpoint = (Theta(j) + Theta(j+1))/2;
            end
            Ig = 1:4;
            Ig = Ig(minlim(g,:) <= Theta(j) & maxlim(g,:) >= Theta(j+1));
            if (~isempty(Ig))
                casenums = nan(size(Ig));
                case1 = sign(coefs(g,:,dcoef)-coefs(g,:,bcoef)) ~= sign(refpoint - midpoint(g,:));
                case2 = ~case1;
                case3 = myround01(coefs(g,:,dcoef) - coefs(g,:,bcoef)) == 0;
                case1 = case1 & ~case3;
                case2 = case2 & ~case3;
                case4 = case3 & myround01(coefs(g,:,acoef) - coefs(g,:,ccoef)) == 0;
                case3 = case3 & ~case4;
                casenums(case1) = 1;
                casenums(case2) = 2;
                casenums(case3) = 3;   
                casenums(case4) = 4;
                maxdf_g(g) = wg(g).*max(squeeze(coefs(g,:,Acoef))'.*maxdf(casenums', squeeze(coefs(g,:,acoef:dcoef)))+squeeze(coefs(g,:,Bcoef))');
                
            end
        end
        maxdfall = sum(maxdf_g);
        stepsize = tolerance/maxdfall;
        if (stepsize < 1000*eps)
            fprintf('warning: step size close to numeric limit\n');
        end
        gridsize = ceil((Theta(j+1)-Theta(j))/stepsize);
        num_batches = ceil(gridsize / MAXVALS);
        for batch = 1:num_batches
            starti = (batch-1)*MAXVALS;
            endi = min(batch*MAXVALS-1,gridsize-1);
            gridpart = (Theta(j)+stepsize*starti):stepsize:min(Theta(j)+stepsize*(endi-1),Theta(j+1));
            if (batch == num_batches)
                gridpart = [gridpart, Theta(j+1)];
            end
            obj_all = calcO3_all(gridpart,z);
            [bestobjpart, bestloc] = min(obj_all);
            if (bestobjpart < bestobj)
                bestobj = bestobjpart;
                besta0 = gridpart(bestloc);
                besta1 = rg(z) + qg(z)*gridpart(bestloc);
            end
        end
   
    end
end
fprintf('\n');
final_bestobj = bestobj;
info.besta0 = besta0;
info.besta1 = besta1;
info.best_a0 = besta0;%backward compatibility
info.best_a1 = besta1;%backward compatibility
info.rg = rg;
info.qg = qg;
info.v = v;
info.wg = wg;
info.pg = pg;
info.pi_g = pi_g;
info.beta = beta;
info.tolerance = tolerance;

[~, a0g, a1g] = calcO2(besta0,besta1);
info.besta0g = a0g;
info.besta1g = a1g;
a0g_nonan = a0g;
a0g_nonan(isnan(a0g)) = 0;
a1g_nonan = a1g;
a1g_nonan(isnan(a1g)) = 0;
info.totalerror = wg' * ((1-pi_g).*a0g_nonan' + pi_g.*a1g_nonan');
info.totalerror0 = (wg' * ((1-pi_g).*a0g_nonan')) / (wg' *(1- pi_g));
info.totalerror1 = (wg' * (pi_g.*a1g_nonan'))/ (wg'*pi_g);
info.unfairness = calcunfairness(wg, pi_g, a0g, a1g, besta0, besta1);


    function val = myeps(a,b) %b can be a vector
        %fix numerical issues
        a = myround01(a);
        b = myround01(b);
        
        if (a == 0)
            val = b;
            return;
        end
        if (a == 1)
            val = 1-b;
            return;
        end
        less = b<a;
        larger = ~less;
        val = nan(size(b));
        val(less) = 1-b(less)/a;
        %b >=a
        val(larger) = 1-(1-b(larger))/(1-a);
    end



    function val = myeps_all_rounded(a,b)
        %a and b are vectors or matrices of the same size
        
        val = nan(size(b));
        val(a==0) = b(a==0);
        done = (a==0);
        val(a==1) = 1-b(a==1);
        done = done | (a==1);
        less = (b<a)& (~done);
        
        val(less) = 1-b(less)./a(less);
        
        larger = (~less) & (~done);
        %b >=a
        val(larger) = 1-(1-b(larger))./(1-a(larger));
    end



    function [obj,a0g, a1g] = calcO2(a0, a1)
        %set S_g
        a0 = myround01(a0);
        a1 = myround01(a1);
        
        if ((a1 > 1) || (a1 < 0) || (a0 > 1) || (a1 < 0))
            fprintf('OUT OF BOUNDS!\n');
            return;
        end
        for g=find(mygplus')
            ming = min(1, myround01(max(0,-rg(g)/qg(g))));
            maxg = max(0,myround01(min(1,(1-rg(g))/qg(g))));
            sg{g} = [ming; maxg; a0; (a1-rg(g))/qg(g)];
            sg{g} = sg{g}((sg{g} >= ming) & (sg{g} <= maxg));
            if (isempty(sg{g}))
                fprintf('problem\n');
            end
        end
        
        val = zeros(1,numg);
        a0g = inf(1,numg);
        a1g = inf(1,numg);
        for g=find(mygplus')
            [tempval,minimizer] = min((1-pi_g(g)).*(beta*myeps(a0,sg{g})+(1-beta)*sg{g})...
                +pi_g(g).*(beta*myeps(a1,rg(g) + qg(g)*sg{g})+(1-beta)*(rg(g)+qg(g)*sg{g})));
            val(g) = wg(g)*tempval;
            a0g(g) = myround01(sg{g}(minimizer));
            a1g(g) = myround01(rg(g)+qg(g)*a0g(g));
            if (val(g) < 0)
                fprintf('negative value!!!\n');
                return;
            end
        end
        for g=find(pi_g'==1)
            val(g) = wg(g)*(beta*myeps(a1,1-pg(g))+(1-beta)*(1-pg(g)));
            a1g(g) = 1-pg(g);
            a0g(g) = nan;
        end
        for g=find(pi_g'==0)
            val(g) = wg(g)*(beta*myeps(a0,pg(g))+(1-beta)*pg(g));
            a0g(g) = pg(g);
            a1g(g) = nan;
        end
        
        obj = sum(val);
    end








    function [obj_all] = calcO2_all(pairs)
        
        a0 = pairs(:,1);
        a1 = pairs(:,2);
        a0 = myround01(a0);
        a1 = myround01(a1);
        numpairs = length(a0);
        if (any(a1 > 1) || any(a1 < 0) || any(a0 > 1) || any(a1 < 0))
            fprintf('OUT OF BOUNDS!\n');
            return;
        end
        sg = nan(numg, numpairs, 4);
        for g=find(mygplus')
            %a0 should be in the 3rd place, (a1-rg(g))/qg(g) should
            %be in the 4th place
            ming = min(1,myround01(max(0,-rg(g)/qg(g))));
            maxg = max(0,myround01(min(1,(1-rg(g))/qg(g))));
            if (ming > maxg)
                fprintf('problem\n');
            end
            assert(ming <= maxg);
            sg(g,:,1) = ming;
            sg(g,:,2) = maxg;
            sg(g,:,3) = max(min(a0,maxg),ming);
            sg(g,:,4) = max(min((a1-rg(g))/qg(g),maxg),ming);
        end
                
        obj_all = zeros(1,numpairs);
        myeps0 = nan(numg,numpairs,4);
        myeps1 = nan(numg,numpairs,4);
        myeps0small = nan(numg,numpairs);
        myeps1small = nan(numg,numpairs);
         for g=find(mygplus')
            sg_flat = nan(numpairs,4);
            sg_flat(:,:) = sg(g,:,:);
            myeps0(g,:,:) = myeps_all_rounded(repmat(a0,1,4),sg_flat);
            myeps1(g,:,:) = myeps_all_rounded(repmat(a1,1,4),myround01(rg(g) + qg(g)*sg_flat));
         end
         for g=find(pi_g'==1)
             myeps1small(g,:) = myeps_all_rounded(a1, repmat(1-pg(g),numpairs,1));
         end
         for g=find(pi_g'==0)
             myeps0small(g,:) = myeps_all_rounded(a0, repmat(pg(g),numpairs,1));
         end
        
         val = zeros(numpairs, numg);
        
            for g=find(mygplus')
                val(:,g) = wg(g)*min((1-pi_g(g)).*(beta*myeps0(g,:,:)+(1-beta)*sg(g,:,:))...
                    +pi_g(g).*(beta*myeps1(g,:,:)+(1-beta)*(rg(g)+qg(g)*sg(g,:,:))), [], 3);
                if (any(isnan(val(:,g))))
                    fprintf('NAN!\n');
                end
                   
                if (val(:,g) < 0)
                    fprintf('negative value!!!\n');
                    return;
                end
            end
            for g=find(pi_g'==1)
                val(:,g) = wg(g)*(beta*myeps1small(g,:)+(1-beta)*(1-pg(g)));
                
            end
            for g=find(pi_g'==0)
                val(:,g) = wg(g)*(beta*myeps0small(g,:)+(1-beta)*pg(g));
                
            end
            
            obj_all = sum(val');
            
        if (any(isnan(obj_all)))
            fprintf('NAN objective!!!\n');
        end
    end



    function obj = calcO3(a0, z)
        obj = calcO2(a0, rg(z) + qg(z)*a0);
        if (isnan(obj))
            fprintf('NAN objective!!!\n');
        end
    end

    function obj = calcO3_all(a0, z)
        obj = calcO2_all([a0', rg(z) + qg(z)*a0']);
    end

    function val = maxdf(casenums, allcoefs) % each row is one set of coefs
        val = nan(size(casenums));
        val(casenums == 1) = daux(allcoefs(casenums==1,:));
        val(casenums == 2) = daux([1-allcoefs(casenums == 2,1)-allcoefs(casenums == 2,2),allcoefs(casenums == 2,2),1-allcoefs(casenums == 2,3)-allcoefs(casenums == 2,4), allcoefs(casenums == 2,4)]);
        val(casenums == 3) = allcoefs(casenums == 3,4)./abs(allcoefs(casenums == 3,1)-allcoefs(casenums == 3,3));
        val(casenums == 4) = 0;
    end

    function val = daux(allcoefs)
        a = allcoefs(:,1);
        b = allcoefs(:,2);
        c = allcoefs(:,3);
        d = allcoefs(:,4);
        
        val = zeros(size(allcoefs,1),1);
        case1 = (d==0) & (c>0);
        val(case1) = b(case1)./c(case1);
        
        %%in next two cases, taking eps distance to avoid ad=bc up to
        %%numerical instability
        case2 =  (d>0) & (a.*d > c.*b+1000*eps); %(a./b > c./d);
        val(case2) = d(case2).^2./abs(b(case2).*c(case2) - a(case2).*d(case2));
        
        case3 = ((d>0) & (a.*d < c.*b-1000*eps) & (b > d));
        val(case3) = abs(b(case3).*c(case3) - a(case3).*d(case3))./...
            (a(case3)+b(case3).*(c(case3)-a(case3))./(b(case3)-d(case3))).^2;
        
    end

    function rounded = myround01(x) %for numberical stability. can work on vectors as well
        rounded = x;
        rounded(abs(x-1) <= 1000*eps) = 1;
        rounded(abs(x) <= 1000*eps) = 0;
    end

end
