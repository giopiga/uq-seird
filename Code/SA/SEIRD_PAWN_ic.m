function [NySol_peak]=SEIRD_PAWN_ic(data,pDomain,icDomain,tcdf,QoIs,plots)
    arguments
        data,pDomain,icDomain,tcdf,QoIs {mustBeText} = {'allStates'},
        plots double {mustBeInteger} = 1,        
    end
    
    tRange=data.tRange;
    nStates=data.nStates;
    tLoc=data.tLoc;
    z=data.z;
    
    if tLoc>max(tRange) || tLoc<0
        error("*** Error: tLoc must be in the considered time interval.");
    elseif z<0 || z>1
        error("*** Error: z must be in [0,1]");
    else

        disp(" --- PAWN sensitivity analysis (IC as random params) ---");
        disp('  QOIs:');
        for i=1:length(QoIs)
            disp(['     ',QoIs{i}]);
        end
        disp(" ");

        p_labels={'beta','r','d','i','S0'};
        s_labels={'S','E','I','R','D'};
        nP=length(p_labels);
        nMC=100;

        % sample uniformly
        SEIRDmc_ic(nMC,nP,pDomain,icDomain,data,QoIs,0);

%% Step 1: params distributions
        disp("   > STEP 1: params definition");
        Nu=1e4;

        xmin=[pDomain(1,:),icDomain(1,:)];
        xmax=[pDomain(2,:),icDomain(2,:)];
        
        % collection of actual QoIs
        s_labels_QoIstates={};
        all='allStates';
        k=1;
        for j=1:length(s_labels)
           S=s_labels{j};
           str='state';
           if sum(strcmp(QoIs,join([str,S])))==1 || sum(strcmp(QoIs,all))==1
               s_labels_QoIstates=[s_labels_QoIstates,S];
               QoI_idx(k)=j;
               k=k+1;
           end
        end
        
        peakT=0;
        peakV=0;
        
        s_labels_QoIpeak={};
        if sum(strcmp(QoIs,'peakTime'))>0
            s_labels_QoIpeak=[s_labels_QoIpeak,'peakTime'];
            peakT=1;
        end
        if sum(strcmp(QoIs,'peakValue'))>0
            s_labels_QoIpeak=[s_labels_QoIpeak,'peakValue'];
            peakV=1;
        end

%% Step 2-3-4: CDFs - KSnorm - plots
        disp("   > STEP 2: evaluation of CDFs");
        disp("   > STEP 3: evaluation of KS distances");

        Nu=1000; % number of samples to estimate unconditional CDF
        Nc=100;  % number of samples to estimate conditional CDF
        n=10;    % number of values to which we fix each params
        
        Nc=100 ; % number of samples to estimate conditional CDFs
        n=10 ; % number of conditioning points
        
        data.tRange=data.tRange(data.tRange<=max(tcdf));
        [NySol]=SEIRDmc_ic(Nu,nP,pDomain,icDomain,data,QoIs,0,0);
        
        if length(s_labels_QoIstates)>0
            for t=1:length(tcdf)
                T=num2str(tcdf(t));

                NySol_t=mySqueeze(NySol(tcdf(t),QoI_idx,:),1);

                if plots==0
                    disp(' ');
                    s=['           - t = ',T,' -         '];
                    for j=1:length(s_labels_QoIstates)
                        s=[s,'<strong>',s_labels_QoIstates{j},'</strong>           '];
                    end

                    disp(s)
                end

                for p=1:nP-1
                    P=p_labels{p};
                    
                    for j=1:n

                        pFix=(pDomain(2,p)-pDomain(1,p))*rand(1)+pDomain(1,p);

                        for i=1:Nc
                            p0(1)=pDomain(1,1)+(pDomain(2,1)-pDomain(1,1))*rand(1);
                            p0(2)=pDomain(1,2)+(pDomain(2,2)-pDomain(1,2))*rand(1);
                            p0(3)=pDomain(1,3)+(pDomain(2,3)-pDomain(1,3))*rand(1);
                            p0(4)=pDomain(1,4)+(pDomain(2,4)-pDomain(1,4))*rand(1);
                            ic0(1)=icDomain(1,1)+(icDomain(2,1)-icDomain(1,1))*rand(1);
                            ic0(2)=1-ic0(:,1);
                            ic0(3)=0*ic0(:,1);
                            ic0(4)=0*ic0(:,1);
                            ic0(5)=0*ic0(:,1);

                            p0(p)=pFix;

                            data.y0=ic0;
                            NySol_c(:,:,i,j)=SEIRDsolver(p0,data);	
                        end

                        NySol_c_t(:,:,j)=mySqueeze(NySol_c(tcdf(t),QoI_idx,:,j),1);

                    end      

                    [X,Fu,Fc]=CDFemp(NySol_t,NySol_c_t);

                    if plots>0
                        KS=KSnorm(Fc,Fu,s_labels_QoIstates,P,0);

                        figure(20+nP*(t-1)+(p-1))
                        CDFplot(X,Fc,'k',1,'',s_labels_QoIstates,KS);

                        hold on                            
                        eval(['Title="cdf - t=',T,' : ',P,'";']);
                        CDFplot(X,Fu,'r',2,Title,s_labels_QoIstates,KS);
                    else
                        KSnorm(Fc,Fu,s_labels_QoIstates,P);
                    end
                end

                P='S0';
                pFix=(icDomain(2)-icDomain(1))*rand(n,1)+icDomain(1);
                
                for j=1:n
                    
                    for i=1:Nc
                        p0(1)=pDomain(1,1)+(pDomain(2,1)-pDomain(1,1))*rand(1);
                        p0(2)=pDomain(1,2)+(pDomain(2,2)-pDomain(1,2))*rand(1);
                        p0(3)=pDomain(1,3)+(pDomain(2,3)-pDomain(1,3))*rand(1);
                        p0(4)=pDomain(1,4)+(pDomain(2,4)-pDomain(1,4))*rand(1);
                        ic0(1)=icDomain(1,1)+(icDomain(2,1)-icDomain(1,1))*rand(1);
                        ic0(2)=1-ic0(:,1);
                        ic0(3)=0*ic0(:,1);
                        ic0(4)=0*ic0(:,1);
                        ic0(5)=0*ic0(:,1);

                        ic0(1)=pFix(j);

                        data.y0=ic0;
                        NySol_c(:,:,i,j)=SEIRDsolver(p0,data);	
                    end

                    NySol_c_t(:,:,j)=mySqueeze(NySol_c(tcdf(t),QoI_idx,:,j),1);
                end

                [X,Fu,Fc]=CDFemp(NySol_t,NySol_c_t);

                if plots>0

                    KS=KSnorm(Fc,Fu,s_labels_QoIstates,P,0);

                    figure(19+nP*(t))
                    CDFplot(X,Fc,'k',1,'',s_labels_QoIstates,KS);

                    hold on                            
                    eval(['Title="cdf - t=',T,' : ',P,'";']);
                    CDFplot(X,Fu,'r',2,Title,s_labels_QoIstates,KS);
                else
                    KSnorm(Fc,Fu,s_labels_QoIstates,P);
                    disp(" ");
                end            
            end
        end
        
        if length(s_labels_QoIpeak)>0
            
            if peakT==0
                [NySol_pV]=max(mySqueeze(NySol(:,3,:),2));
                NySol_peak=NySol_pV;
            elseif peakV==0
                [~,idxNySol_pT]=max(mySqueeze(NySol(:,3,:),2));
                NySol_pT=data.tRange(idxNySol_pT);
                NySol_peak=NySol_pT;
            else
                [NySol_pV,idxNySol_pT]=max(mySqueeze(NySol(:,3,:),2));
                NySol_pT=data.tRange(idxNySol_pT);
                NySol_peak=[NySol_pT;NySol_pV];
            end
            
            if plots==0
                disp(' ');
                
                for j=1:length(s_labels_QoIpeak)
                    s=['<strong>',s_labels_QoIpeak{j},'</strong>',repelem(' ',1,13)];
                end

                disp([repelem(' ',1,20),s]);
            end
            
            for p=1:nP-1
                P=p_labels{p};
                    
                for j=1:n
                    
                    pFix=(pDomain(2,p)-pDomain(1,p))*rand(1)+pDomain(1,p);

                    for i=1:Nc
                        p0(1)=pDomain(1,1)+(pDomain(2,1)-pDomain(1,1))*rand(1);
                        p0(2)=pDomain(1,2)+(pDomain(2,2)-pDomain(1,2))*rand(1);
                        p0(3)=pDomain(1,3)+(pDomain(2,3)-pDomain(1,3))*rand(1);
                        p0(4)=pDomain(1,4)+(pDomain(2,4)-pDomain(1,4))*rand(1);
                        ic0(1)=icDomain(1,1)+(icDomain(2,1)-icDomain(1,1))*rand(1);
                        ic0(2)=1-ic0(:,1);
                        ic0(3)=0*ic0(:,1);
                        ic0(4)=0*ic0(:,1);
                        ic0(5)=0*ic0(:,1);

                        p0(p)=pFix;

                        data.y0=ic0;
                        NySol_c(:,:,i,j)=SEIRDsolver(p0,data);	
                    end
                    
                    if peakT==0
                        [NySol_c_pV]=max(mySqueeze(NySol_c(:,3,:,j),2));
                        NySol_c_peak(:,:,j)=NySol_c_pV;
                    elseif peakV==0
                        [~,idxNySol_c_pT]=max(mySqueeze(NySol_c(:,3,:,j),2));
                        NySol_c_pT=data.tRange(idxNySol_c_pT);
                        NySol_c_peak(:,:,j)=NySol_c_pT;
                    else
                        [NySol_c_pV,idxNySol_c_pT]=max(mySqueeze(NySol_c(:,3,:,j),2));
                        NySol_c_pT=data.tRange(idxNySol_c_pT);
                        NySol_c_peak(:,:,j)=[NySol_c_pT;NySol_c_pV];
                    end

                end      
                
                [X,Fu,Fc]=CDFemp(NySol_peak,NySol_c_peak);

                if plots>0
                    KS=KSnorm(Fc,Fu,s_labels_QoIpeak,P,0);

                    figure(100+p)
                    CDFplot(X,Fc,'k',1,'',s_labels_QoIpeak,KS);

                    hold on                            
                    eval(['Title="cdf - peak time and value : ',P,'";']);
                    CDFplot(X,Fu,'r',2,Title,s_labels_QoIpeak,KS);
                else
                    KSnorm(Fc,Fu,s_labels_QoIpeak,P,'');
                end
            end
            
            P='S0';
            pFix=(icDomain(2)-icDomain(1))*rand(n,1)+icDomain(1);
            
            for j=1:n
                
                for i=1:Nc
                    p0(1)=pDomain(1,1)+(pDomain(2,1)-pDomain(1,1))*rand(1);
                    p0(2)=pDomain(1,2)+(pDomain(2,2)-pDomain(1,2))*rand(1);
                    p0(3)=pDomain(1,3)+(pDomain(2,3)-pDomain(1,3))*rand(1);
                    p0(4)=pDomain(1,4)+(pDomain(2,4)-pDomain(1,4))*rand(1);
                    ic0(1)=icDomain(1,1)+(icDomain(2,1)-icDomain(1,1))*rand(1);
                    ic0(2)=1-ic0(:,1);
                    ic0(3)=0*ic0(:,1);
                    ic0(4)=0*ic0(:,1);
                    ic0(5)=0*ic0(:,1);

                    ic0(1)=pFix(j);

                    data.y0=ic0;
                    NySol_c(:,:,i,j)=SEIRDsolver(p0,data);	
                end

                if peakT==0
                    [NySol_c_pV]=max(mySqueeze(NySol_c(:,3,:,j),2));
                    NySol_c_peak(:,:,j)=NySol_c_pV;
                elseif peakV==0
                    [~,idxNySol_c_pT]=max(mySqueeze(NySol_c(:,3,:,j),2));
                    NySol_c_pT=data.tRange(idxNySol_c_pT);
                    NySol_c_peak(:,:,j)=NySol_c_pT;
                else
                    [NySol_c_pV,idxNySol_c_pT]=max(mySqueeze(NySol_c(:,3,:,j),2));
                    NySol_c_pT=data.tRange(idxNySol_c_pT);
                    NySol_c_peak(:,:,j)=[NySol_c_pT;NySol_c_pV];
                end

            end
                
            [X,Fu,Fc]=CDFemp(NySol_peak,NySol_c_peak);

            if plots>0
                
                KS=KSnorm(Fc,Fu,s_labels_QoIpeak,P,0);
                
                figure(100+nP+1)
                CDFplot(X,Fc,'k',1,'',s_labels_QoIpeak,KS);

                hold on                            
                eval(['Title="cdf - peak time and value : ',P,'";']);
                CDFplot(X,Fu,'r',2,Title,s_labels_QoIpeak,KS);
            else
                KSnorm(Fc,Fu,s_labels_QoIpeak,P);
                disp(" ");
            end           
        end
                 
        if plots>0
            disp("   > STEP 4: plots");
        else
            disp("   > STEP 4: plots - missed by user options (plots=0)");
        end

%% End
        disp(" ");
        disp(" --- End PAWN sensitivity analysis --- ");
        disp(" ");
    end
end