clear
clc

global ddevs

% First, we bandpass-filter the data from Russia to get the quarterly deviations
load revision2_output.mat

tlength = length(r2_russia_log_output);
ttime = transpose(0:1:tlength-1);

lin_trend = @(x) x(1) + x(2)*ttime;

const_guess = r2_russia_log_output(1);
slope_guess = (r2_russia_log_output(end)-r2_russia_log_output(1))/tlength;

squared_devs = @(x) sum( (lin_trend(x) - r2_russia_log_output).^2);

init_guess = [const_guess;slope_guess];
opt_lin = fminunc(squared_devs,init_guess);

ttrend = lin_trend(opt_lin);

lin_devs = r2_russia_log_output - ttrend;

% Now perform the band-pass filter on the log-linear deviations
ddevs = bpass(lin_devs,6,32);

save REVISION2_BP_devs_gdp.mat ddevs

% Next, we estimate an AR(1) process for the deviations

init_param_vec = [.9, .03];

ML_est = fmincon(@neg_likelihood,init_param_vec,[],[],[],[],[0.0,0.0],[1.0,0.2],[]);

% save REVISION2_BP_GDP_EST_THROW_OUT_TAILS.mat ML_est

function [negLL] = neg_likelihood(param_vec)
    global ddevs
    
    Nobservs = length(ddevs);
    
    rrhoy = param_vec(1);
    ssigmay = param_vec(2);
    
    logLL = 0.0;
    
    for i=2:Nobservs
        logLL = logLL + log( normpdf(ddevs(i),rrhoy*ddevs(i-1),ssigmay));
    end
    
    negLL = -logLL;
    
end