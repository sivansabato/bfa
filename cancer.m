cancer_data %load data
numprob = size(norm_data,2);
pi_g_all = ground_truth / 10^5;
pg_all = norm_data;


wg = states / sum(states);
betas = 0:0.01:1;

numbetas = length(betas);

allunfairness = nan(numbetas,numprob);
alltotalerror = nan(numbetas,numprob);
alltotalerror0 = nan(numbetas, numprob);
alltotalerror1 = nan(numbetas, numprob);
allprecision = nan(numbetas, numprob);
allrecall = nan(numbetas, numprob);
pg_total = nan(1,numprob);
min_goal = zeros(numbetas, numprob);
allmybetas = nan(numbetas, numprob);

all_pi_total = zeros(numprob, 1);


clear info
pi_total = zeros(1,numprob);


for prob = 1:numprob
    fprintf('\n\n\ncondition: %s %d\n\n', relevant_conditions{prob}, prob);
    pi_g = pi_g_all(:,prob);
    pi_total(prob) = wg' * pi_g;
    foundsomething = false;
    
    pg = pg_all(:,prob);
    pg_total(prob) = wg' * pg;
    fprintf('true positives: %f, predicted positives: %f, max ratio - %f\n', pi_total(prob), pg_total(prob), max(pi_total(prob)/pg_total(prob), pg_total(prob)/pi_total(prob)))
    if ((pg_total(prob) < 0.5*pi_total(prob)) || (pg_total(prob) > 2*pi_total(prob)))
        fprintf('prediction is too far off, skipping\n');
        continue;
    end
    foundsomething = true;

    for betai = 1:length(betas)
        beta = betas(betai);
        fprintf('setting beta to %g\n', beta);
        tolerance = pi_total(prob)/10;
        fprintf('tolerance = %g\n', tolerance);
        [min_goal(betai, prob), info(betai, prob)] = find_lb(wg, pg, pi_g, beta, tolerance);
        allunfairness(betai, prob) = info(betai, prob).unfairness;
        alltotalerror(betai, prob) = info(betai, prob).totalerror;
        alltotalerror0(betai, prob) = info(betai, prob).totalerror0;
        alltotalerror1(betai, prob) = info(betai, prob).totalerror1;
        allprecision(betai,prob) = (pi_total(prob)./pg_total(prob))*(1-alltotalerror1(betai,prob));
        allrecall(betai,prob) = 1-alltotalerror1(betai,prob);
    end
    if foundsomething
        allmybetas(:,prob) = betas';
        fprintf('true positives: %f\n', pi_total(prob))
        fprintf('beta\t   min_goal \t unfairness \t  unfairness/pos     total_error\t       error0 \t      error1\t    precision\t    recall \n');
        fprintf('%.2f\t %12f\t %12f\t %12f\t %12f\t %12f\t %12f\t %12f\t %12f\n', [allmybetas(:,prob), min_goal(:,prob), allunfairness(:, prob), allunfairness(:, prob)/pi_total(prob),...
            alltotalerror(:, prob), alltotalerror0(:, prob), alltotalerror1(:, prob), allprecision(:,prob), allrecall(:,prob)]');
    end
end


for prob = 1:numprob
    fprintf('\n\n\ncondition: %s, index %d\n', relevant_conditions{prob}, prob);
    fprintf('true positives: %f, predicted positives: %f\n', pi_total(prob), pg_total(prob))
    fprintf('beta\t   min_goal \t unfairness \t  unfairness/pos     total_error\t       error0 \t      error1\t    precision\t    recall \n');
    fprintf('%.2f\t %12f\t %12f\t %12f\t %12f\t %12f\t %12f\t %12f\t %12f\n', [allmybetas(:,prob), min_goal(:,prob), allunfairness(:, prob), allunfairness(:, prob)/pi_total(prob),...
        alltotalerror(:, prob), alltotalerror0(:, prob), alltotalerror1(:, prob), allprecision(:,prob), allrecall(:,prob)]');
end



   filename = sprintf('cancer_results_%s.csv', datetime('now','Format','yyyy-MM-dd__HH_mm_ss'));
     filemat = zeros(numbetas, 3*numprob);
    for prob = 1:numprob
        if (isnan(allunfairness(1,prob)))
            continue;
        end
        firstnan = find(isnan(allmybetas(:,prob)),1);    
        if (isempty(firstnan))
            firstnan = size(allunfairness,1)+1;
        end
        filemat(1:(firstnan-1),prob*3-2) = allmybetas(1:(firstnan-1),prob);
        filemat(1:(firstnan-1),prob*3-1) = allunfairness(1:(firstnan-1),prob)/pi_total(prob);
        filemat(1:(firstnan-1),prob*3) = alltotalerror(1:(firstnan-1),prob)/pi_total(prob);
        if (~isempty(firstnan))
            filemat(firstnan:end,[prob*3-2:prob*3]) = repmat(filemat(firstnan-1,[prob*3-2:prob*3]),size(allunfairness,1)-firstnan+1,1);
        end
        
    end
    csvwrite(filename, filemat);
