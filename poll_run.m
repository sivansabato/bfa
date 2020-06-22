load('polls_data.mat');

%these are the dates of the polls in the loaded poll data, in order
dates = {'Nov-06','Oct-13','Sep-01','Nov-06','Nov-07', 'Sep-13', 'Nov-07',	'Oct-26','Sep-01', 'Oct-27'};

order_by_dates = [3,9,6,2,8,10,1,4,5,7];
polls = polls(:,[1,2,order_by_dates+2]);
dates = dates(order_by_dates);
pop_total = sum(polls(:,1));
wg_orig = polls(:,1)/pop_total;

num_probs = size(polls,2)-2;

clear obj
clear info
clear pi_total

pi_total = zeros(1,num_probs);
betas = 0:0.01:1;
obj = inf(length(betas), num_probs);
for i = 1:num_probs 
    poll = polls(:,i+2);
    actual = polls(:,2);
    nan_in_data = isnan(poll) | isnan(actual);
    poll = poll(~nan_in_data);
    actual = actual(~nan_in_data);
    wg = wg_orig(~nan_in_data);
    
    pi_g = actual/100;
    pg = poll/100;
    pi_total(i) = wg'*pi_g;
    tolerance = 1/100;
    for b = 1:length(betas)
        beta = betas(b);
        [obj(b,i), info(b,i)] = find_lb(wg, pg, pi_g, beta, tolerance, 1);
        fprintf('problem %d: beta = %f, unfairness = %f, error = %f, unfairness/prob = %f, error/prob = %f\n', ...
            i, info(b,i).beta, info(b, i).unfairness, info(b, i).totalerror, info(b, i).unfairness/pi_total(i), info(b,i).totalerror/pi_total(i));
    end
end
filemat = zeros(length(betas), num_probs*2);
for i = 1:num_probs
    filemat(:,i*2-1) = [info(:,i).unfairness];
    filemat(:,i*2) = [info(:,i).totalerror];
end
csvwrite(['polls_' date '.csv'], filemat);
