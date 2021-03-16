function [yPred,chain]=SEIRDmcmc_ic(paramLSQ,data,nMCMC,C,pDomain,nPred,plots)
    arguments
        paramLSQ,data,nMCMC,C,pDomain,nPred double {mustBeInteger} = 100, 
        plots double {mustBeInteger} = 1, 
    end
    
% MCMC model
    model.ssfun=@SEIRDss;
    model.sigma2=data.sigma2;

% MCMC params
%     params = {
%     {'\beta',paramLSQ(1),pDomain(1,1),pDomain(2,1)}
%     {'r',paramLSQ(2),pDomain(1,2),pDomain(2,2)}
%     {'d',paramLSQ(3),pDomain(1,3),pDomain(2,3)}
%     {'i',paramLSQ(4),pDomain(1,4),pDomain(2,4)} };
    params = {
        {'\beta',paramLSQ(1),0,Inf}
        {'r',paramLSQ(2),0,Inf}
        {'d',paramLSQ(3),0,Inf}
        {'i',paramLSQ(4),0,Inf} };


% MCMC options
    options.nsimu=nMCMC; % number of samples
    options.qcov=C; % (initial) proposal covariance
    options.method='dram'; % method: DRAM
    options.adaptint=100; % adaptation interval
    options.verbosity=0;
    options.burnintime=3000;

    [result,chain]=mcmcrun(model,data,params,options);
    
% prediction
%     chainBurnt=chain(ceil(length(chain)/2):end,:);
    chainBurnt=chain(:,:);
    if nPred>0
        disp("   > STEP 3: Prediction");
        
        idx=ceil(linspace(1,length(chainBurnt),nPred));
        for j=1:nPred
            y=SEIRDsolver(chainBurnt(idx(j),:),data);
            yPred(:,:,j)=y(:,data.idxSynt);            
        end
        
    else
        disp("   .Predictions missed by the user (pred=0)");
    end
    
% plots
    if plots>0
        figure(201);
        mcmcplot(chain,[],result.names);
        figure(202);
        mcmcplot(chain,[],result.names,'pairs');
        
        if nPred>0
            st=["S","E","I","R","D"];
            figure(203)
            SEIRDplot(data.tRange,yPred,'',0);
            hold on
            %plot(data.tRange,data.SEIRDmean(:,data.idxSynt),'k-','LineWidth',2);
            for j=1:length(data.idxSynt)
                plot(data.tRange,data.Max(:,j),'k--','LineWidth',2);
                plot(data.tRange,data.Min(:,j),'k--','LineWidth',2);
            end
            legend(st(data.idxSynt));
            title('MC evaluatioin of SEIRD model  -  posterior');      
        end
        
    else
        disp("   .Plots missed by the user (plots=0)");
    end
    %chainstats(chain,result)
    
end