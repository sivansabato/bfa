%To obtain the data for this experiment, go to http://projects.fivethirtyeight.com/general-model/president_general_polls_2016.csv
%You will need to reorganize the data from this table into a matlab matrix where each row is a state
%and the columns are:
%1. Total votes (from: https://en.wikipedia.org/wiki/2016_United_States_presidential_election)
%2. Percentage voting for Clinton (from: https://en.wikipedia.org/wiki/2016_United_States_presidential_election)
%3. Poll results from the above-mentioned file, in the following order (dates denote end dates of the poll):
%       a. YouGov, 6/11/2016
%       b. Ipsos, 10/13/2016
%       c. Ipsos, 9/1/2016
%       d. Ipsos, 6/11/2016
%       e. Google Consumer Surveys, 11/7/2016
%       f. Google Consumer Surveys, 9/13/2016	
%       g. SurveyMonkey, 11/7/2016
%       h. SurveyMonkey, 10/26/2016
%       i. SurveyMonkey, 9/1/2016
%       j. Ipsos, 10/27/2016
%
%From each poll take the column adjpoll_clinton
% 
%Once reorganized, the data should be in a 51x12 (51 states by 12 polls) matrix called "polls"

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
