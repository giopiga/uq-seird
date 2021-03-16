%%
% In this main you can find a complete workflow for the analysis 
% of SEIRD model for epidemics:
%   - FwdUQ by MC analysis, based on uniform prior
%   - SA by VBSA (Saltelli) and PAWN (KS norm) methods for some QoIs
%                                       (states, peak time and peak value)
%   - InvUQ by MCMC algorithm (DRAM)
%   - FwdUQ by MC analysis, based on the posterior
%
%   Each subfolder could also be run independently, paying attention to
%   define all the needed variables
%
% Our project was based on the paper: 
%   C. Piazzola, L. Tamellini, R. Tempone, 
%   "A note on tools for predictionunder uncertainty and identifiability
%   of SIR-like dynamical systems forepidemiology"
%
%   suggesting us some aspects of the workflow and the data of the
%   parameters
%
%   M.Cesaratto, G.Pigani - march, 2021
%

run SEIRDcode_start.m

clear
close all
clc

%% General settings
% model setting
data.y0=[0.95;0.05;0;0;0];  % useful only if y0 are not random params
data.tRange=(0:0.5:150);
data.nStates=length(data.y0);
data.z=0;                % perc. intensity of a possible lockdown
data.tLoc=0;             % day of the lockdown

% parameters bounds (from the paper)
beta_min=0.25; beta_max=0.35;
r_min=0.06; r_max=0.18;
d_min=0.01; d_max=0.02;
i_min=0.14; i_max=0.33;
pDomain=[beta_min r_min d_min i_min; beta_max r_max d_max i_max];

% bounds of IC as params (centered in 0.95, avoiding =1)
% in all the folders, the functions terminating by _ic.m are supposed to
% receive random initial conditions, here defined between 0.9 and 0.99
S0_min=0.9; S0_max=0.99;
icDomain=[S0_min; S0_max];

% QoIs
QoIs={'allStates','peakTime','peakValue'};        % Possibility: stateX (or allStates), peakTime, peakValue 
                                                    % stateX, where X=S,E,I,R or D
                                                    % default: allStates

%% Forward analysis (MC)
% first visualization of SEIRD trajectories
disp(" --- System estimation (MC) ---");
disp(" ");

nMC=150;
data.z=0;                % perc. intensity
data.tLoc=0;             % time
[~,SEIRDmean,Max,Min]=SEIRDmc_ic(nMC,5,pDomain,icDomain,data,QoIs);
% [~,SEIRDmean]=SEIRDmc(nMC,5,pDomain,data,QoIs);

%%
% possibility to add a lockdown

data.z=0.7;                % perc. intensity
data.tLoc=15;             % time

toDisp=strcat(" --- System estimation (MC) : lockdown at day ",num2str(data.tLoc)," ---");
disp(toDisp);
disp(" ");
% [~,SEIRDmean]=SEIRDmc(nMC,5,pDomain,data,QoIs);
[~,SEIRDmean]=SEIRDmc_ic(nMC,5,pDomain,icDomain,data,QoIs);

%% Sensitivity analysis        
tcdf=[15,50,data.tRange(end)];                    % set of time at which
                                                  % PAWN idx are evaluated
SEIRD_PAWN_ic(data,pDomain,icDomain,tcdf,QoIs);
% SEIRD_PAWN(data,pDomain,tcdf,QoIs)

%% Inverse analysis (MCMC)
my_dir = fileparts(which("SEIRDcode_start.m")); 
addpath(strcat(my_dir,"\Code\InvUQ\MCMCstat"));

disp(" --- MCMC inverse analysis ---");
disp(" ");

data.sigma2=0.01;           % variance of syntethic data
data.y0=[0.95;0.05;0;0;0];  % better to set again their values after PAWN_ic

beta0=0.28;
r0=0.11;
d0=0.018;
i0=0.18;
param0=[beta0,r0,d0,i0];   % from paper

statesSynt={'I','D'};      % possibilities: 'S','E','I','R','D'

% Generation of synthetic data
disp("   > STEP 0: generation of synthetic data");
[ySol_synt,idxSynt,peakV_synt,peakT_synt]=SEIRDsynt(param0,data,statesSynt);
    
% Optimal estimators - LSQ
disp("   > STEP 1: looking for optimal estimator (LSQ)");
disp(" ");

param0=pDomain(1,:);
data.ySol_synt=ySol_synt;
data.idxSynt=idxSynt;
paramLSQ=SEIRDlsq(param0,data);    

% MCMC algorithm
nMCMC=1e4;
C=data.sigma2*eye(4);
disp(" ");
disp("   > STEP 2: MCMC algorithm");
data.Max=Max(:,idxSynt);
data.Min=Min(:,idxSynt);

% here MCMCstat package is required
nMC=100;
[yPred,chain]=SEIRDmcmc_ic(paramLSQ,data,nMCMC,C,pDomain,nMC);

rmpath(strcat(my_dir,"\Code\InvUQ\MCMCstat"));

disp("   > STEP 4.a: Posterior pdfs of the peaks");

[peakV_post,idx_peakT_post]=max(mySqueeze(yPred(:,1,:),2),[],1);
peakT_post=data.tRange(idx_peakT_post);

[xT_post,pdfT_post,muT_post,sigmaT_post]=PDFemp(peakT_post);
[xV_post,pdfV_post,muV_post,sigmaV_post]=PDFemp(peakV_post);

SEIRDy=SEIRDmc(nMC,5,pDomain,data,QoIs,0,0);
[peakV_seird,idx_peakT_seird]=max(mySqueeze(SEIRDy(:,3,:),2),[],1);
peakT_seird=data.tRange(idx_peakT_seird);

[xT_seird,pdfT_seird,muT_seird,sigmaT_seird]=PDFemp(peakT_seird);
[xV_seird,pdfV_seird,muV_seird,sigmaV_seird]=PDFemp(peakV_seird);

figure()
subplot(1,2,1)
plot(xT_post,pdfT_post,'-','LineWidth',2);
hold on
plot(xT_seird,pdfT_seird,'-','LineWidth',2);
xline(data.tRange(peakT_synt),'k--','LineWidth',2);
title('peakTime');
axis on; grid on;
legend('posterior','prior','real value')
subplot(1,2,2)
plot(xV_post,pdfV_post,'-','LineWidth',2);
hold on
plot(xV_seird,pdfV_seird,'-','LineWidth',2);
xline(peakV_synt,'k--','LineWidth',2);
title('peakValue');
axis on; grid on;
legend('posterior','prior','real value')
sgtitle('Peak PDF - posterior');

disp("   > STEP 4.b: Posterior pdfs of the params");
idx_post=ceil(linspace(3000,size(chain,1),7000));
param_post_value=chain(idx_post,:);

beta_post=param_post_value(:,1);
r_post=param_post_value(:,2);
d_post=param_post_value(:,3);
i_post=param_post_value(:,4);

[x_post_beta,pdf_post_beta,mu_beta,s_beta]=PDFemp(beta_post);
[x_post_r,pdf_post_r,mu_r,s_r]=PDFemp(r_post);
[x_post_d,pdf_post_d,mu_d,s_d]=PDFemp(d_post);
[x_post_i,pdf_post_i,mu_i,s_i]=PDFemp(i_post);

figure()
subplot(2,2,1)
plot(x_post_beta,pdf_post_beta,'-','LineWidth',2);
hold on
plot([beta_min,beta_max],[1/(beta_max-beta_min),1/(beta_max-beta_min)],'-','LineWidth',2);
xline(beta0,'k--','LineWidth',2);
xline(mu_beta,'k-','LineWidth',2);
grid on;
legend('posterior','prior','real value','post mean');
title('beta');

subplot(2,2,2)
plot(x_post_r,pdf_post_r,'-','LineWidth',2);
hold on
plot([r_min,r_max],[1/(r_max-r_min),1/(r_max-r_min)],'-','LineWidth',2);
xline(r0,'k--','LineWidth',2);
xline(mu_r,'k-','LineWidth',2);
grid on;
legend('posterior','prior','real value','post mean');
title('r');

subplot(2,2,3)
plot(x_post_d,pdf_post_d,'-','LineWidth',2);
hold on
plot([d_min,d_max],[1/(d_max-d_min),1/(d_max-d_min)],'-','LineWidth',2);
xline(d0,'k--','LineWidth',2);
xline(mu_d,'k-','LineWidth',2);
grid on;
legend('posterior','prior','real value','post mean');
title('d');

subplot(2,2,4)
plot(x_post_i,pdf_post_i,'-','LineWidth',2);
hold on
plot([i_min,i_max],[1/(i_max-i_min),1/(i_max-i_min)],'-','LineWidth',2);
xline(i0,'k--','LineWidth',2);
xline(mu_i,'k-','LineWidth',2);
grid on;
legend('posterior','prior','real value','post mean');
title('i');
sgtitle('Params PDF - posterior');
