% This code extracts policy/pricing functions from the fortran code for the
% purpose of plotting figures and analyzing simulations

clear
clc

% orig_dir = cd;

% cd 20e-4

global nDims nInputs gridPointsX gridPointsM
nDims = 3;
nInputs = nDims*20;

% Here are the relevant model parameters we will require
r = .01;
coup = r;
llambda = .066 - coup;
ssigma_y = .0226;
rho_y = .9212;
ssigma_uncond = ssigma_y/sqrt(1.0-rho_y^2);
bbeta = .963;
sov_CRRA = 2.0;


% Load in all the necessary files to initialize the variables
fileID = fopen('bounds.txt','r');
bbounds = fscanf(fileID,'%f',[2,5]);
fclose(fileID);

bL = bbounds(1,1);
bH = bbounds(2,1);
yL = bbounds(1,2);
yH = bbounds(2,2);
xL = bbounds(1,3);
xH = bbounds(2,3);
mL = bbounds(1,4);
mH = bbounds(2,4);
VL = bbounds(1,5);
VH = bbounds(2,5);

transM = @(x) (x - [yL;bL;mL])./[yH-yL;bH-bL;mH-mL];
transX = @(x) (x - [yL;bL;xL])./  [yH-yL;bH-bL;xH-xL];
trans = @(x) (x - [yL;bL])./[yH-yL;bH-bL];

invTransM = @(x) [yL;bL;mL] + x.*[yH-yL;bH-bL;mH-mL];
invTransX = @(x) [yL;bL;xL] + x.*[yH-yL;bH-bL;xH-xL];
invTrans = @(x) [yL;bL] + x.*[yH-yL;bH-bL];

ythresh = 0.9; % break below which default is free
MP = .225; % maximum proportional default cost

ylev = -MP*ythresh/(yH-ythresh); % Default cost level parameter
ycurv = MP/(yH-ythresh); % Default cost curvature parameter
ydef = @(y) y - max(ylev*y + ycurv*y^2.0,0.0);


fileID = fopen('means.txt','r');
mmeans = fscanf(fileID,'%f',[1,8]);
fclose(fileID);

VR_mean = mmeans(1);
VD_mean = mmeans(2);
A_mean = mmeans(3);
I_mean = mmeans(4);
q_mean = mmeans(5);
qd_mean = mmeans(6);
qRN_mean = mmeans(7);
qdRN_mean = mmeans(8);

fileID = fopen('stds.txt','r');
stds = fscanf(fileID,'%f',[1,8]);
fclose(fileID);

VR_std = stds(1);
VD_std = stds(2);
A_std = stds(3);
I_std = stds(4);
q_std = stds(5);
qd_std = stds(6);
qRN_std = stds(7);
qdRN_std = stds(8);

nPlots = 101;

b_grid = bL: (bH-bL)/(nPlots-1) : bH;
y_grid = yL: (yH-yL)/(nPlots-1) : yH;
x_grid = xL: (xH-xL)/(nPlots-1) : xH;
m_grid = mL: (mH - mL)/(nPlots-1) : mH;

% Load in all the necessary files to initialize the variables
fileID = fopen('grid_points_X.txt','r');
gridPointsX = fscanf(fileID,'%f',[nDims,nInputs]);
fclose(fileID);

fileID = fopen('grid_points_M.txt','r');
gridPointsM = fscanf(fileID,'%f',[nDims,nInputs]);
fclose(fileID);

% %
fileID = fopen('VD_kern_coeffs.txt','r');
xx = fscanf(fileID,'%f',[1,2]);
fclose(fileID);
VD_kernCoeff = xx';

fileID = fopen('q_kern_coeffs.txt','r');
xx = fscanf(fileID,'%f',[1,2]);
fclose(fileID);
q_kernCoeff = xx';

fileID = fopen('qd_kern_coeffs.txt','r');
xx = fscanf(fileID,'%f',[1,2]);
fclose(fileID);
qd_kernCoeff = xx';

fileID = fopen('qRN_kern_coeffs.txt','r');
xx = fscanf(fileID,'%f',[1,2]);
fclose(fileID);
qRN_kernCoeff = xx';

fileID = fopen('qdRN_kern_coeffs.txt','r');
xx = fscanf(fileID,'%f',[1,2]);
fclose(fileID);
qdRN_kernCoeff = xx';

fileID = fopen('VR_kern_coeffs.txt','r');
xx = fscanf(fileID,'%f',[1,2]);
fclose(fileID);
VR_kernCoeff = xx';

fileID = fopen('A_kern_coeffs.txt','r');
xx = fscanf(fileID,'%f',[1,2]);
fclose(fileID);
A_kernCoeff = xx';

fileID = fopen('I_kern_coeffs.txt','r');
xx = fscanf(fileID,'%f',[1,2]);
fclose(fileID);
I_kernCoeff = xx';

fileID = fopen('VD_gpr_coeffs.txt','r');
VD_gprCoeff = fscanf(fileID,'%f',[nInputs,1]);
fclose(fileID);

fileID = fopen('q_gpr_coeffs.txt','r');
q_gprCoeff = fscanf(fileID,'%f',[nInputs,1]);
fclose(fileID);

fileID = fopen('qd_gpr_coeffs.txt','r');
qd_gprCoeff = fscanf(fileID,'%f',[nInputs,1]);
fclose(fileID);

fileID = fopen('qRN_gpr_coeffs.txt','r');
qRN_gprCoeff = fscanf(fileID,'%f',[nInputs,1]);
fclose(fileID);

fileID = fopen('qdRN_gpr_coeffs.txt','r');
qdRN_gprCoeff = fscanf(fileID,'%f',[nInputs,1]);
fclose(fileID);

fileID = fopen('VR_gpr_coeffs.txt','r');
VR_gprCoeff = fscanf(fileID,'%f',[nInputs,1]);
fclose(fileID);

fileID = fopen('A_gpr_coeffs.txt','r');
A_gprCoeff = fscanf(fileID,'%f',[nInputs,1]);
fclose(fileID);

fileID = fopen('I_gpr_coeffs.txt','r');
I_gprCoeff = fscanf(fileID,'%f',[nInputs,1]);
fclose(fileID);

% cd(orig_dir);

I_pol = @(x) 1/(1.0 + exp(-I_std*gpr_approx(trans(x),I_kernCoeff,I_gprCoeff,1)-I_mean));
A_pol = @(x) bL + (bH-bL)/(1.0 + exp(-A_std*gpr_approx(trans(x),A_kernCoeff,A_gprCoeff,1)-A_mean));

q = @(x) 1/(1.0 + exp(-q_std*gpr_approx(transX(x),q_kernCoeff,q_gprCoeff,2)-q_mean));
qd = @(x) 1/(1.0 + exp(-qd_std*gpr_approx(transM(x),qd_kernCoeff,qd_gprCoeff,3)-qd_mean));

qRN = @(x) 1/(1.0 + exp(-qRN_std*gpr_approx(transX(x),qRN_kernCoeff,qRN_gprCoeff,2)-qRN_mean));
qdRN = @(x) 1/(1.0 + exp(-qdRN_std*gpr_approx(transM(x),qdRN_kernCoeff,qdRN_gprCoeff,3)-qdRN_mean));

VR = @(x) VL + (VH-VL)/(1.0 + exp(-VR_std*gpr_approx(trans(x),VR_kernCoeff,VR_gprCoeff,1)-VR_mean));
VD = @(x) VL + (VH-VL)/(1.0 + exp(-VD_std*gpr_approx(transM(x),VD_kernCoeff,VD_gprCoeff,3)-VD_mean));

defFun = @(x) VR(x(1:2)) < VD(x);

CEC_0 = (1-sov_CRRA)*(VR([1.0;0])*(1-bbeta))^(1/(1-sov_CRRA));
CEC_erg = (1-sov_CRRA)*(VR([1.0;.61])*(1-bbeta))^(1/(1-sov_CRRA));
CEC_crisis = (1-sov_CRRA)*(VR([0.95;.71])*(1-bbeta))^(1/(1-sov_CRRA));

VR_plot = b_grid;
VD_plot = b_grid;
for ib = 1:nPlots
    VR_plot(ib) = VR([1.0;b_grid(ib)]);
    VD_plot(ib) = VD([1.0;b_grid(ib);0.651]);
end

figure
plot(b_grid,VR_plot,b_grid,VD_plot,'--')
title('Value Across b')
xlabel('B_t')
legend('VR','VD')

VR_plot2 = y_grid;
VD_plot2 = y_grid;
for iy = 1:nPlots
    VR_plot2(iy) = VR([y_grid(iy);0.6]);
    VD_plot2(iy) = VD([y_grid(iy);0.6;0.651]);
end

figure
plot(y_grid,VR_plot2,y_grid,VD_plot2,'--')
title('Value Across y')
xlabel('y_t')
legend('VR','VD')

VR_plot3 = m_grid;
VD_plot3 = m_grid;
for im = 1:nPlots
    VR_plot3(im) = VR([1.0;0.6]);
    VD_plot3(im) = VD([0.9;0.6;m_grid(im)]);
end

figure
plot(m_grid,VR_plot3,m_grid,VD_plot3,'--')
title('Value Across m')
xlabel('m_t')
legend('VR','VD')


q_plot1 = b_grid;
q_plot2 = b_grid;
q_plot3 = b_grid;

for ib = 1:nPlots
    q_plot1(ib) = q([1.0-2*ssigma_uncond;b_grid(ib);0.652]);
    q_plot2(ib) = q([1.0;b_grid(ib);0.652]);
    q_plot3(ib) = q([1.0+2*ssigma_uncond;b_grid(ib);0.652]);
end

yield_plot_q1 = (llambda + coup + (1-llambda).*q_plot1)./q_plot1;
yield_plot_q2 = (llambda + coup + (1-llambda).*q_plot2)./q_plot2;
yield_plot_q3 = (llambda + coup + (1-llambda).*q_plot3)./q_plot3;

sprd_plot1 = yield_plot_q1.^4 - (1+r)^4;
sprd_plot2 = yield_plot_q2.^4 - (1+r)^4;
sprd_plot3 = yield_plot_q3.^4 - (1+r)^4;

figure
plot(b_grid,q_plot1,'r-.',b_grid,q_plot2,'b-',b_grid,q_plot3,'g--',...
    'LineWidth',3)
xlabel('$B_{t+1}$','Interpreter','LaTeX','FontSize',30)
ylabel('$q_{t}$','Interpreter','LaTeX','FontSize',30)
title('Price Across $y_t$','Interpreter','LaTeX','FontSize',30)
legend('$y_t = 1.0-2 \sigma_{y,uncond}$','$y_t = 1.0$','$y_t = 1.0 + 2 \sigma_{y,uncond}$','Interpreter','LaTeX','FontSize',25)
axis([bL bH 0 1.1])

figure
plot(b_grid,sprd_plot1,b_grid,sprd_plot2,'-.',b_grid,sprd_plot3,'--')
xlabel('b_{t+1}')
title('Spread Across y')
legend('yL','ySS','yH')

q_plotRN = q_plot2;
for ib = 1:nPlots
    q_plotRN(ib) = qRN([1.0;b_grid(ib);0.651]);
end

yield_plot_RN = (llambda + coup + (1-llambda).*q_plotRN)./q_plotRN;
sprd_plotRN = yield_plot_RN.^4 - (1+r)^4;

figure
plot(b_grid,sprd_plot2,b_grid,sprd_plotRN,'--')
xlabel('b_{t+1}')
title('RN Share')
legend('Benchmark','Risk-Neutral ZNS')


q_plot1 = b_grid;
q_plot2 = b_grid;
q_plot3 = b_grid;

for ib = 1:nPlots
    q_plot1(ib) = q([1.0-2*ssigma_uncond;b_grid(ib);xL + (.5-2/3*.5)*(xH-xL)]);
    q_plot2(ib) = q([1.0-2*ssigma_uncond;b_grid(ib);xL + .5*(xH-xL)]);
    q_plot3(ib) = q([1.0-2*ssigma_uncond;b_grid(ib);xL + (.5+2/3*.5)*(xH-xL)]);
end

yield_plot_q1 = (llambda + coup + (1-llambda).*q_plot1)./q_plot1;
yield_plot_q2 = (llambda + coup + (1-llambda).*q_plot2)./q_plot2;
yield_plot_q3 = (llambda + coup + (1-llambda).*q_plot3)./q_plot3;

sprd_plot1 = yield_plot_q1.^4 - (1+r)^4;
sprd_plot2 = yield_plot_q2.^4 - (1+r)^4;
sprd_plot3 = yield_plot_q3.^4 - (1+r)^4;

figure
plot(b_grid,q_plot1,'r-.',b_grid,q_plot2,'b-',b_grid,q_plot3,'g--',...
    'LineWidth',3)
xlabel('$B_{t+1}$','Interpreter','LaTeX','FontSize',30)
ylabel('$q_{t}$','Interpreter','LaTeX','FontSize',30)
title('Price: $y_t = 1.0-2 \sigma_{y,uncond}$','Interpreter','LaTeX','FontSize',30)
legend('$x_t = \bar{x}-2 \sigma_{x}$','$x_t = \bar{x}$','$x_t = \bar{x} + 2 \sigma_{x}$','Interpreter','LaTeX','FontSize',25)
axis([bL bH 0 1.1])

figure
plot(b_grid,sprd_plot1,'r-.',b_grid,sprd_plot2,'b-',b_grid,sprd_plot3,'g--',...
    'LineWidth',3)
xlabel('$B_{t+1}$','Interpreter','LaTeX','FontSize',30)
ylabel('Spread$_{t}$','Interpreter','LaTeX','FontSize',30)
title('Spreads: $y_t = 1.0-2 \sigma_{y,uncond}$','Interpreter','LaTeX','FontSize',30)
legend('$x_t = \bar{x}-2 \sigma_{x}$','$x_t = \bar{x}$','$x_t = \bar{x} + 2 \sigma_{x}$','Interpreter','LaTeX','FontSize',25)

q_plot1 = x_grid;


for ix = 1:nPlots
    q_plot1(ix) = q([.85;.78;x_grid(ix)]);
end


figure
plot(x_grid,q_plot1)
xlabel('x_t')
title('Price Across x')
axis([bL bH 0 1.1])


A_plot1 = b_grid;
A_plot2 = b_grid;
A_plot3 = b_grid;
for ib = 1:nPlots
    A_plot1(ib) = A_pol([1.0-2*ssigma_uncond;b_grid(ib)]);
    A_plot2(ib) = A_pol([1.0;b_grid(ib)]);
    A_plot3(ib) = A_pol([1.0+2*ssigma_uncond;b_grid(ib)]);
end

figure
plot(b_grid,A_plot1,'r-',b_grid,A_plot3,'g--','LineWidth',3);
hold on
plot(b_grid,b_grid,'k.')
title('Policy Function Across $y_t$','Interpreter','LaTeX','FontSize',30)
xlabel('$B_t$','Interpreter','LaTeX','FontSize',30)
ylabel('$B_{t+1}$','Interpreter','LaTeX','FontSize',30)
legend('$y_t = 1.0-2 \sigma_{y,uncond}$','$y_t = 1.0 + 2 \sigma_{y,uncond}$','Interpreter','LaTeX','FontSize',25)
axis([bL bH 0 bH])

I_plot1 = b_grid;
I_plot2 = b_grid;
I_plot3 = b_grid;
for ib = 1:nPlots
    I_plot1(ib) = I_pol([1.0-2*ssigma_uncond;b_grid(ib)]);
    I_plot2(ib) = I_pol([1.0;b_grid(ib)]);
    I_plot3(ib) = I_pol([1.0+2*ssigma_uncond;b_grid(ib)]);
end

figure
plot(b_grid,I_plot1,'r-.',b_grid,I_plot2,'b-',b_grid,I_plot3,'g--',...
    'LineWidth',3)
title('Information Policy Function','Interpreter','LaTeX','FontSize',30)
xlabel('$B_{t+1}$','Interpreter','LaTeX','FontSize',30)
legend('$y_t = 1.0-2 \sigma_y$','$y_t = 1.0$','$y_t = 1.0 + 2 \sigma_y$','Interpreter','LaTeX','FontSize',25)
axis([bL bH 0 1.1])

%% Now for the simulations
% This is used to generate additional moments (such as the CVR) as well as
% explore simulated events

% cd 20e-4

Nsimul = 500000;
fileID = fopen('simulations.txt','r');
simul_mat = transpose(fscanf(fileID,'%f',[11,Nsimul]));
fclose(fileID);

% cd(orig_dir);

sim_y = simul_mat(:,1);
sim_b = simul_mat(:,2);
sim_m = simul_mat(:,3);
sim_x = simul_mat(:,4);
sim_access = simul_mat(:,5);
sim_q = simul_mat(:,6);
sim_sprd = simul_mat(:,7);
sim_rho = simul_mat(:,8);
sim_sprdRN = simul_mat(:,9);
sim_issuing = simul_mat(:,10);
sim_DP = simul_mat(:,11);

% Check the info fraction
sim_rho_issuing = sim_rho(sim_issuing > 0.5);
attn_frac = length(sim_rho_issuing(sim_rho_issuing > 0.5*max(sim_rho_issuing)))/length(sim_rho_issuing);


% First, compute a host of cyclicalities and moments
sim_sprd_issuing = sim_sprd(sim_issuing > 0.5);
sim_sprdRN_issuing = sim_sprdRN(sim_issuing > 0.5);
sim_y_issuing = sim_y(sim_issuing > 0.5);
sim_b_issuing = sim_b(sim_issuing > 0.5);
sim_q_issuing = sim_q(sim_issuing > 0.5);

RN_share_issuing = sim_sprdRN_issuing./sim_sprd_issuing;

sim_de_issuing = sim_y_issuing(1:end-1) - sim_b_issuing(1:end-1)*(llambda + coup) + ...
    + sim_q_issuing(1:end-1).*(sim_b_issuing(2:end) - sim_b_issuing(1:end-1)*(1-llambda));

sim_tb_issuing = sim_y_issuing(1:end-1) - sim_de_issuing;

avg_de_share = mean(sim_de_issuing./sim_y_issuing(1:end-1));
avg_tb_share = mean(sim_tb_issuing./sim_y_issuing(1:end-1));

tb_vol = std(sim_tb_issuing);

de_cyc = corr(sim_de_issuing,sim_y_issuing(1:end-1));

c_y_vol_ratio = std(sim_de_issuing)/std(sim_y_issuing);

sim_def = (sim_issuing(1:end-1) > 0.5).*(sim_issuing(2:end) < 0.5);
sim_def = [sim_def;0];
y_defs = sim_y( sim_def > 0.5 );

by_before_default = mean(sim_b(sim_def > 0.5)./sim_y(sim_def > 0.5));

% Share of defaults below ss-output
recession_def_share = sum(y_defs < 1.0)/length(y_defs);

% Next we plot dynamics of the `most typical' of the defaults
Ndefs = length(y_defs);
tt = transpose(1:1:Nsimul);
tdefs = tt.*(sim_def > 0.5);
tdefs = tdefs(tdefs > 0.5); % These are the periods in which the country defaults

LL_list = tdefs;

% Conditional on defaulting, we want to find the most likely default event
tdef_typical = tdefs(1);
def_window = 16; % explore def_window periods before default
Lbest = -500000;
sum_cond_def_prob = 0.0;
for t=1:Ndefs
    % In this case, the country had access def_window periods before
    % default
    sum_cond_def_prob = sum_cond_def_prob + sim_DP(tdefs(t)-1);
    if (sum( sim_issuing(tdefs(t)-def_window+1:tdefs(t)) ) == def_window) && ...
            (sum(sim_rho(tdefs(t)-def_window+1:tdefs(t)-4)) < .1)
        LL = 0.0;
        for t1=tdefs(t)-def_window+2:tdefs(t)-4

            LL = LL + log( normpdf( (log(sim_y(t1)) - rho_y*log(sim_y(t1-1)))/ssigma_y  ) );

        end
        LL_list(t) = LL;
        if LL > Lbest
            tdef_typical = tdefs(t);
            Lbest = LL;
        end
    end
end
cond_def_prob = sum_cond_def_prob/Ndefs;

tplot = -def_window+1:1:0;
def_y_path = tplot;
def_b_path = tplot;
def_info_path = tplot;
def_sprd_path = tplot;

for t=1:def_window
    def_y_path(t) = sim_y(tdef_typical-def_window+t);
    def_b_path(t) = sim_b(tdef_typical-def_window+t);
    def_info_path(t) = sim_rho(tdef_typical-def_window+t);
    def_sprd_path(t) = sim_sprd(tdef_typical-def_window+t);
end

figure
plot(tplot,def_y_path,'LineWidth',3)
xlabel('$t$ (Default at $0$)','Interpreter','LaTeX','FontSize',30)
ylabel('$y_{t}$','Interpreter','LaTeX','FontSize',30)

figure
plot(tplot,def_b_path,'LineWidth',3)
xlabel('$t$ (Default at $0$)','Interpreter','LaTeX','FontSize',30)
ylabel('$B_{t}$','Interpreter','LaTeX','FontSize',30)

figure
plot(tplot,def_info_path,'LineWidth',3)
xlabel('$t$ (Default at $0$)','Interpreter','LaTeX','FontSize',30)
ylabel('$\rho_{t}$','Interpreter','LaTeX','FontSize',30)

figure
plot(tplot,def_sprd_path,'LineWidth',3)
xlabel('$t$ (Default at $0$)','Interpreter','LaTeX','FontSize',30)
ylabel('Spread$_{t}$','Interpreter','LaTeX','FontSize',30)


% Now, compute the CVR
% First, rank the spreads
sprd_change_issuing = sim_sprd_issuing(2:end) - sim_sprd_issuing(1:end-1);

crisis_thresh = prctile(sprd_change_issuing,97.5);

wwindow = 5;
wwindow_plot = 15;

sim_sprd_issuing_with_zeros = sim_sprd.*sim_issuing;
sim_sprd_change_with_zeros = sim_sprd_issuing_with_zeros(2:end) - sim_sprd_issuing_with_zeros(1:end-1);
ttemp = (2:1:Nsimul)'.*(sim_sprd_change_with_zeros > crisis_thresh);
crisis_periods = ttemp(ttemp > 0.5);

crisis_periods = crisis_periods(crisis_periods > wwindow_plot);
crisis_periods = crisis_periods(crisis_periods < Nsimul-wwindow_plot);

typical_crisis_risk_neutral_share = zeros(2*wwindow_plot,length(crisis_periods));
typical_crisis_risk_premium = zeros(2*wwindow_plot,length(crisis_periods));
typical_crisis_info = zeros(2*wwindow_plot,length(crisis_periods));

% Compute the CVR AND search for the most typical crisis
best_LL = -1.0e8;

LL_list = crisis_periods;
most_likely_ind = 1;

sum_ratios = 0;
num_crises = 0;
for t=1:length(crisis_periods)
    if(sum( sim_issuing( crisis_periods(t)-wwindow:crisis_periods(t)+wwindow-1) ) == 2*wwindow) % Make sure country is always issuing
        num_crises = num_crises + 1;

        sum_ratios = sum_ratios + std(sim_sprd(crisis_periods(t):crisis_periods(t)+wwindow-1))/...
            std(sim_sprd(crisis_periods(t)-wwindow:crisis_periods(t)-1));

        typical_crisis_risk_neutral_share(:,num_crises) = sim_sprdRN(crisis_periods(t)-wwindow_plot:crisis_periods(t)+wwindow_plot-1)./...
            sim_sprd(crisis_periods(t)-wwindow_plot:crisis_periods(t)+wwindow_plot-1);
        typical_crisis_risk_premium(:,num_crises) = (sim_sprd(crisis_periods(t)-wwindow_plot:crisis_periods(t)+wwindow_plot-1) - sim_sprdRN(crisis_periods(t)-wwindow_plot:crisis_periods(t)+wwindow_plot-1) ); 
        typical_crisis_info(:,num_crises) = sim_rho(crisis_periods(t)-wwindow_plot:crisis_periods(t)+wwindow_plot-1);

        LL = 0.0;
        for t1 = crisis_periods(t)-wwindow+1:crisis_periods(t)+wwindow-1

            LL = LL + log( normpdf( (log(sim_y(t1)) - rho_y*log(sim_y(t1-1)))/ssigma_y  ) );
        end

        LL_list(t) = LL;
        if LL > Lbest
            most_likely_ind = num_crises;
            tcrisis_typical = crisis_periods(t);
            Lbest = LL;
        end
    end
end
CVR = sum_ratios/num_crises;

RN_share_uncond = median(RN_share_issuing);
RN_share_crisis = median(typical_crisis_risk_neutral_share(wwindow_plot+1,:));
RN_share_crisis_year_later = median(typical_crisis_risk_neutral_share(wwindow_plot+5,:));

average_crisis_risk_neutral_share = typical_crisis_risk_neutral_share(:,most_likely_ind);
average_crisis_risk_premium = typical_crisis_risk_premium(:,most_likely_ind);
average_crisis_info = typical_crisis_info(:,most_likely_ind);

for t=1:length(average_crisis_risk_neutral_share)
    average_crisis_risk_neutral_share(t) = median(typical_crisis_risk_neutral_share(t,:));
    average_crisis_risk_premium(t) = median(typical_crisis_risk_premium(t,:));
    average_crisis_info(t) = median(typical_crisis_info(t,:));
end


event_axis = -wwindow_plot:1:wwindow_plot-1;

figure
plot(event_axis,average_crisis_risk_neutral_share)
xlabel('Time (Crisis at $0$)','Interpreter','LaTeX');
ylabel('Average Risk-Neutral Share','Interpreter','LaTeX')
title('Typical Crisis','Interpreter','LaTeX')

figure
plot(event_axis,average_crisis_info)
xlabel('Time (Crisis at $0$)','Interpreter','LaTeX');
ylabel('Average Information','Interpreter','LaTeX')
title('Typical Crisis','Interpreter','LaTeX')

figure
plot(event_axis,average_crisis_risk_premium)
xlabel('Time (Crisis at $0$)','Interpreter','LaTeX');
ylabel('Average Risk Premium','Interpreter','LaTeX')
title('Typical Crisis','Interpreter','LaTeX')

