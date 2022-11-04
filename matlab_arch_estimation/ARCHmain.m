% This code uses quasi-MLE to derive the parameters for a GARCH(1,1)
% process 
% "The pricing of sovereign risk under costly information"
% Grace Gu and Zach Stangebye, 7/11/2018

clear
clc

% First, guess initial parameters
% AR parameters
mmu_s = .031;
rrho_s = .6;

% GARCH innovation parameters
oomega_sqrt = .006;
aalpha = .6;

x_init = [mmu_s; rrho_s; oomega_sqrt; aalpha];

xLB = [-.1; 0; 0.0; 0.0];
xUB = [.1; .99; .015; 0.999];

% Find the parameters vis MLE (given initial conditions, likelihood can be
% expressed easily)
% load data_rus1.mat
load data_rus1.mat

ddata_spread = data_rus1/10000;

obj1 = @(x) ARCHneg_LL(x,ddata_spread);

ooptions = optimset('Display','Iter','TolFun',1e-10);
[ssoln, obj2,~,~,~,~,hess1] = fmincon(obj1,x_init,[],[],[],[],xLB,xUB,[],ooptions);

SE = sqrt(diag(inv(hess1)));

conf_bands = [ssoln-1.96*SE, ssoln, ssoln+1.96*SE];

save REVISION2_one_alpha_no_gamma.mat

%% Now, we simulate the model and compute the CVR
load REVISION2_one_alpha_no_gamma.mat

rng(124);

mmu_s = .0307; % ssoln(1);
rrho_s = ssoln(2);

% GARCH innovation parameters
oomega = .00375^2; %ssoln(3)^2;
aalpha = .95 ; % ssoln(4); % Change this to run comparative statics for CVR

T =   1501000; %  10000; %  % Length of simulation

innov = zeros(T,1);
simul_sprd = zeros(T,1);
var_proc = zeros(T,1);

% Initialize variance process
simul_sprd(1) = mmu_s;
var_proc(2) = oomega/(1-aalpha);
for t=2:T-1
    innov(t) = normrnd(0,sqrt(var_proc(t)));
    simul_sprd(t) = (1-rrho_s)*mmu_s + rrho_s*simul_sprd(t-1) + innov(t);
    var_proc(t+1) = oomega + aalpha*innov(t)^2 ;
end


% Now, compute the CVR
w = 5;
simul_sprd = simul_sprd(1001:T);
T = T-1000; % burn in

crisis_ind = 2:1:T-1;
sprd_change = simul_sprd(2:T-1)-simul_sprd(1:T-2);
sorted_sprd = sort(sprd_change);
% sprd_change_thresh = quantile(sorted_sprd,.971);
sprd_change_thresh = quantile(sorted_sprd,.942);

crisis_ind = crisis_ind'.*(sprd_change > sprd_change_thresh);
crisis_ind(crisis_ind < 1) = [];

numCrises = length(crisis_ind);

roll_sum = 0.0;

if crisis_ind(numCrises)+w > T
    crisis_ind(numCrises) = [];
    numCrises = numCrises - 1;
end

if crisis_ind(1)-1-w < 1
    crisis_ind(1) = [];
    numCrises = numCrises - 1;
end

for t=1:numCrises
    roll_sum = roll_sum + std(simul_sprd(crisis_ind(t):crisis_ind(t)+w))/std(simul_sprd(crisis_ind(t)-1-w:crisis_ind(t)-1));
end
CVR = roll_sum/numCrises;

[mean(simul_sprd),std(simul_sprd),CVR]

% save REVISION2_one_alpha_no_gamma_with_CVR.mat