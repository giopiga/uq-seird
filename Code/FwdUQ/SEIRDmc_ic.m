function [SEIRDsol_ic,SEIRDmean,Max,Min]=SEIRDmc_ic(nMC,nP,pDomain,icDomain,data,QoIs,plots,pdfs,Legend)
    arguments
        nMC,nP,pDomain,icDomain,data,QoIs {mustBeText} = ''
        plots double {mustBeInteger} = 1,
        pdfs double {mustBeInteger} = 1,
        Legend {mustBeText} = '';
    end
    
    % params and IC sampling
    unit_set=rand(nMC,nP);
    
    MC(:,1)=pDomain(1,1)+(pDomain(2,1)-pDomain(1,1))*unit_set(:,1);
    MC(:,2)=pDomain(1,2)+(pDomain(2,2)-pDomain(1,2))*unit_set(:,2);
    MC(:,3)=pDomain(1,3)+(pDomain(2,3)-pDomain(1,3))*unit_set(:,3);
    MC(:,4)=pDomain(1,4)+(pDomain(2,4)-pDomain(1,4))*unit_set(:,4);
    
    IC(:,1)=icDomain(1,1)+(icDomain(2,1)-icDomain(1,1))*unit_set(:,5);
    IC(:,2)=1-IC(:,1);
    IC(:,3)=0*IC(:,1); IC(:,4)=0*IC(:,1); IC(:,5)=0*IC(:,1);

    % SEIRD evaluations
    for j=1:nMC
        data.y0=IC(j,:);
        [SEIRDsol_ic(:,:,j)]=SEIRDsolver(MC(j,:),data);
    end
    
    SEIRDmean=mean(SEIRDsol_ic,3);
    
    % pdfs of QoIs : peakTime and pealValue
    s_labels_QoIpeak={};
    peakT=0;
    peakV=0;
    if sum(strcmp(QoIs,'peakTime'))>0
        peakT=1;
        s_labels_QoIpeak=[s_labels_QoIpeak,'peakTime'];
    end
    if sum(strcmp(QoIs,'peakValue'))>0
        peakV=1;
        s_labels_QoIpeak=[s_labels_QoIpeak,'peakValue'];
    end
    
    if length(s_labels_QoIpeak)>0
        if peakV==1 && peakT==1
            [peakV,idx_peakT]=max(mySqueeze(SEIRDsol_ic(:,3,:),2));  
            peakT=data.tRange(idx_peakT);
            [xT,fT,muT,sT]=PDFemp(peakT);
            [xV,fV,muV,sV]=PDFemp(peakV);
            
            x=nan(2,max(length(xT),length(xV)));
            x(1,1:length(xT))=xT;
            x(2,1:length(xV))=xV;
            
            f=nan(2,max(length(fT),length(fV)));
            f(1,1:length(fT))=fT;
            f(2,1:length(fV))=fV;
            
        elseif peakT==0
            peakV=max(mySqueeze(SEIRDsol_ic(:,3,:),2));
            [x,f,mu,s]=PDFemp(peakV);
            
        else
            [~,idx_peakT]=max(mySqueeze(SEIRDsol_ic(:,3,:),2));
            peakT=data.tRange(idx_peakT);
            [x,f,mu,s]=PDFemp(peakT);
        end
    end    
    
    % plots
    if plots>0
        figure()
        SEIRDplot(data.tRange,SEIRDsol_ic,["S","E","I","R","D"],'MC evaluatioin of SEIRD model  -  prior: unif',1);
    end

    if length(s_labels_QoIpeak)>0 && pdfs>0
        figure(2)
        hold on
        PDFplot(x,f,s_labels_QoIpeak,Legend);
    end
    
    Max=max(SEIRDsol_ic,[],3);
    Min=min(SEIRDsol_ic,[],3);
    
end
