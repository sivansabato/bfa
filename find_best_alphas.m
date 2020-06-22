function [alpha0,alpha1] = find_best_alphas(ag, pi_g_fact, alpha0g, alpha1g, divisionfactor)

if (~exist('divisionfactor','var'))
    divisionfactor = 1;
end
pi_g = pi_g_fact/divisionfactor;

alphayg(1,:) = alpha0g;
alphayg(2,:) = alpha1g;
endpoints{1} = [0,1, alpha0g];
endpoints{2} = [0,1,alpha1g];
for y = [0,1]
    endpoints{y+1} = unique(endpoints{y+1});
end
bestay = nan(2,1);
bestobjy = inf(2,1);
pi_g_y(1,:) = 1-pi_g;
pi_g_y(2,:) = pi_g;

for y = [0,1]
    for i = 1:length(endpoints{y+1})
        alphay = endpoints{y+1}(i);
        
        
        obj = 0;
        for g = 1:length(ag)
            obj = obj + ag(g)*calcObjective_y(alphay, alphayg(y+1,g), pi_g_y(y+1,g));
        end
        if (obj < bestobjy(y+1))
            bestobjy(y+1) = obj;
            bestay(y+1) = alphay;
        end
    end
end

alpha0 = bestay(1);
alpha1 = bestay(2);


end



function obj = calcObjective_y(alphay, myalphayg, mypi_g_y)
epsgy = calcEps(myalphayg,alphay);

if (~isnan(epsgy))
    obj = mypi_g_y * epsgy;
else
    obj = 0;
end
end



function eps = calcEps(alphag, alpha)
if (alpha == 0)
    eps = alphag;
else
    if (alpha == 1)
        eps = 1-alphag;
    else
        eps = max(1-alphag/alpha, 1- (1-alphag)/(1-alpha));
    end
end
end


