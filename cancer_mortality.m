load('cancer_mortality.mat');

cancer_names = {'Brain','Breast','Colon','Kidney','Liver','Lung','Oral','Pancreatic','Stomach','Thyroid'};

pop_total = sum(cancer_death_data(:,1));
wg_orig = cancer_death_data(:,1)/pop_total;

num_probs = (size(cancer_death_data,2)-1)/2;

clear obj
clear info
clear pi_total

betas = 0:0.01:1;
for i = 1:num_probs 
    disease = cancer_death_data(:,i*2+1);
    death = cancer_death_data(:,i*2);
    nan_in_data = isnan(disease) | isnan(death);
    disease = disease(~nan_in_data);
    death = death(~nan_in_data);
    wg = wg_orig(~nan_in_data);
    avg_disease_per_10e5 = sum(disease.*wg);
    avg_death_per_10e5 = sum(death.*wg);
    death_rate = avg_death_per_10e5/avg_disease_per_10e5;
    
    pi_g = disease * death_rate * 10e-5;
    pg = death *10e-5;
    pi_total(i) = wg'*pi_g;
    tolerance = pi_total(i)/10;
    for b = 1:length(betas)
        beta = betas(b);
        [obj(b,i), info(b,i)] = find_lb(wg, pg, pi_g, beta, tolerance, 1);
        fprintf('%s (%d) beta = %f, unfairness = %f, error = %f, unfairness/prob = %f, error/prob = %f\n', ...
            cancer_names{i}, i, info(b,i).beta, info(b, i).unfairness, info(b, i).totalerror, info(b, i).unfairness/pi_total(i), info(b,i).totalerror/pi_total(i));
    end
end

filemat = zeros(length(betas), num_probs*2);
for i = 1:num_probs
    filemat(:,i*2-1) = [info(:,i).unfairness]/pi_total(i);
    filemat(:,i*2) = [info(:,i).totalerror]/pi_total(i);
end
csvwrite(['cancer_mortality_' date '.csv'], filemat);

