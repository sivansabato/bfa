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
true_unfairness = [];
true_error = [];
classifierindex = zeros(0,2);
betas = [0:0.1:1];
valperm = randperm(length(vals));
numtorun = 10;


for iperm = 1:length(vals)
       i = valperm(iperm);
       fprintf('---------- i = %d\n', i);
       j = randi([1,length(vals{i})]);
        
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


        ag = allres(i,j).ag;
        pg = allres(i,j).pg;
        p = sum(pg .* ag);
        pi_g = allres(i,j).pi_g;
        pi = sum(pi_g .* ag);
        alpha0g = allres(i,j).alpha0g; %this is the true one
        alpha1g = allres(i,j).alpha1g; %this is the true one
        alpha0g(isnan(alpha0g)) = 0; %will be multiplied by 0 anyway
        alpha1g(isnan(alpha1g)) = 0; %will be multiplied by 0 anyway
        [allres(i,j).truea0, allres(i,j).truea1] = find_best_alphas(ag, allres(i,j).pi_g_fact, alpha0g, alpha1g, allres(i,j).divisionfactor);
            
        tmp_true_unfairness = max(1-alpha0g/allres(i,j).truea0, 1-(1-alpha0g)/(1-allres(i,j).truea0)) * ((1-pi_g) .* ag) +...
                max(1-alpha1g/allres(i,j).truea1, 1-(1-alpha1g)/(1-allres(i,j).truea1)) * (pi_g .* ag);
        tmp_true_error = 1-allres(i,j).testaccuracy;

        if (tmp_true_error == 0 && tmp_true_unfairness ==0)
             fprintf('skipping due to triviality\n');
             continue;
       end
       tolerance = 1/100; %can be reduced, but will take longer to run

       current_index = size(unfairness,2)+1;
       true_unfairness(current_index) = tmp_true_unfairness;
       true_error(current_index) = tmp_true_error;
       for betai = 1:length(betas)
            beta = betas(betai); 
            [~, info] = find_lb(allres(i,j).ag, allres(i,j).pg_fact, ...
            allres(i,j).pi_g_fact, beta, tolerance, allres(i,j).divisionfactor);
            unfairness(betai,current_index) = info.unfairness;
            totalerror(betai,current_index) = info.totalerror;
            classifierindex(current_index,1) = i;
            classifierindex(current_index,2) = j;
       end

    
      if (current_index == numtorun)
           break;
      end
end


filemat = zeros(length(betas), numtorun+1);
filemat(:,1) = betas';
for n = 1:numtorun
    i = classifierindex(n,1);
    j = classifierindex(n,2);
    predicted = betas .* unfairness(:,n)'+ (1-betas).*totalerror(:,n)';
    real = betas * true_unfairness(n) + (1-betas)*true_error(n);
    filemat(:,n+1) = predicted./real;
end
filename = sprintf('discdata_%s.csv',  datetime('now','Format','yyyy-MM-dd__HH_mm_ss'));
csvwrite(filename, filemat);

        




