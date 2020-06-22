   fprintf('loading data file...\n');
load('USCensus1990raw.data.mat'); %%%loads data set. It can be downloaded from this link: https://archive.ics.uci.edu/ml/datasets/US+Census+Data+(1990)
    fprintf('loaded.\n');


%%%the split to train and test is random, so results could be different every run.
powstate = 93;
data_screened = data(data(:,powstate) > 0 & data(:,powstate)<90,:);
perm = randperm(size(data_screened,1));
data_screened = data_screened(perm,:);
nsamples = size(data_screened,1);
trainsize = floor(nsamples/2); 


clear catCol
clear allres
clear label
for i=1:size(data_screened,2)
    vals{i} = unique(data(:,i));
    catCol(i) = length(vals{i}) < 10;
    if (~catCol(i))
        [~,vals{i}] = hist(data(:,i),10);
    end
end




unfairness = [];
totalerror = [];


for i = 1:length(vals)
    fprintf('---------- i = %d\n', i);
    for j = 1:length(vals{i})
        
        fprintf('------------ j = %d\n', j);
        
        labelthresh = vals{i}(j);
        if catCol(i)
            Y = (data_screened(:,i) == labelthresh)*2-1;
        else
            Y = (data_screened(:,i) > labelthresh)*2-1;
        end
        if (min(mean(Y==1), 1-mean(Y==1))< 0.01)
            fprintf('skipping due to few different labels\n');
            continue;
        end
        
        X = data_screened;
        X(:,[i,powstate]) = []; %remove these columns
        catColTemp = catCol;
        catColTemp([i,powstate]) = [];
        fprintf('column %d with val %g\n', i, vals{i}(j));
alltemp = calculate_classifier(X,catColTemp, Y, trainsize, data_screened(:,powstate));
        
        allres(i,j).ag = alltemp.ag;
        allres(i,j).pg = alltemp.pg;
        allres(i,j).pi_g = alltemp.pi_g;
        allres(i,j).pg_fact = alltemp.pg_fact;
        allres(i,j).pi_g_fact = alltemp.pi_g_fact;
        allres(i,j).divisionfactor = alltemp.divisionfactor;
        allres(i,j).trainaccuracy = alltemp.trainaccuracy;
        allres(i,j).testaccuracy = alltemp.testaccuracy;
        allres(i,j).defaultaccuracy = alltemp.defaultaccuracy;
        allres(i,j).alpha0g = alltemp.alpha0g;
        allres(i,j).alpha1g = alltemp.alpha1g;
        beta = 1;
        tolerance = 1/100; %can be reduced, but will take longer to run
        [allres(i,j).min_unfairness, allres(i,j).info] = find_lb(allres(i,j).ag, allres(i,j).pg_fact, ...
            allres(i,j).pi_g_fact, beta, tolerance, allres(i,j).divisionfactor);



        if (~isnan(allres(i,j).min_unfairness))
            wg = allres(i,j).info.wg;
            pg = allres(i,j).info.pg;
            p = sum(pg .* wg);
            pi = sum(allres(i,j).info.pi_g .* wg);
            pi_g = allres(i,j).info.pi_g;
            alpha0g = allres(i,j).alpha0g; %this is the true one
            alpha1g = allres(i,j).alpha1g; %this is the true one
            alpha0g(isnan(alpha0g)) = 0; %will be multiplied by 0 anyway
            alpha1g(isnan(alpha1g)) = 0; %will be multiplied by 0 anyway
            [allres(i,j).truea0, allres(i,j).truea1] = find_best_alphas(wg, allres(i,j).pi_g_fact, alpha0g, alpha1g, allres(i,j).divisionfactor);
            
            allres(i,j).true_unfairness = max(1-alpha0g/allres(i,j).truea0, 1-(1-alpha0g)/(1-allres(i,j).truea0)) * ((1-pi_g) .* wg) +...
                max(1-alpha1g/allres(i,j).truea1, 1-(1-alpha1g)/(1-allres(i,j).truea1)) * (pi_g .* wg);
            min_unfairness = allres(i,j).min_unfairness;
            fprintf('true_unfairness: %g, min_unfairness: %g, ratio: %g\n', allres(i,j).true_unfairness, min_unfairness, min_unfairness/allres(i,j).true_unfairness);
        else
            fprintf('NaN in min_unfairness\n');
        end
    end
end
for i = 1:size(allres,1)
    for j = 1:size(allres,2)
        if (isempty(allres(i,j).trainaccuracy))
            allres(i,j).trainaccuracy = NaN;
            allres(i,j).testaccuracy = NaN;
            allres(i,j).defaultaccuracy = NaN;
        end
    end
end
testaccmat = zeros(size(allres));
defaultaccmat = zeros(size(allres));
minspmat = nan(size(allres));
alpha0upperbound = zeros(size(allres));
for i = 1:size(allres,1)
    for j = 1:size(allres,2)
        testaccmat(i,j) = allres(i,j).testaccuracy;
        defaultaccmat(i,j) = allres(i,j).defaultaccuracy;
        if (isempty(allres(i,j).min_unfairness))
            minspmat(i,j) = NaN;
        else
            minspmat(i,j) = allres(i,j).min_unfairness;
        end
    end
end

[I,J] = find(~isnan(minspmat));
inds = sortrows([I,J],1);
I = inds(:,1);
J = inds(:,2);
%%
true_unfairness = inf(1,length(I));
min_unfairness = inf(1,length(I));
true_accuracy = inf(1,length(I));

for n = 1:length(I)
    i = I(n);
    j = J(n);
    if (isempty(allres(i,j).info))
        continue;
    end
    
     
    true_unfairness(n) = allres(i,j).true_unfairness;
    min_unfairness(n) = allres(i,j).min_unfairness;
    true_accuracy(n) = allres(i,j).testaccuracy;
end

%%
filemat = zeros(length(I),6);
filemat(:,1) = true_unfairness';
filemat(:,2) = min_unfairness';
ratios= min_unfairness./true_unfairness;
ratios(true_unfairness == 0 & min_unfairness == 0) = 1;
filemat(:,3) = ratios';
filemat(:,4) = true_accuracy';
filemat(:,5) = I';
filemat(:,6) = J';
filename = sprintf('census_commands_%s.csv', datetime('now','Format','yyyy-MM-dd__HH_mm_ss'));
csvwrite(filename, filemat);

fprintf('Results for unfairness lower bound:\n');
for n = 1:length(I)
    i = I(n);
    j = J(n);
    if (isempty(allres(i,j).info))
        continue;
    end
    fprintf('i=%d,j=%d (categorical: %d), value %g, test acc %g, min_unfairness %g, true_unfairness %g, ratio %g\n', ...
	    i, j, catCol(j), vals{i}(j), allres(i,j).testaccuracy, min_unfairness(n), true_unfairness(n), min_unfairness(n)/true_unfairness(n));
end



 
