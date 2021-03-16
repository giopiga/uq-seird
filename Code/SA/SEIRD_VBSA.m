function [Si,STi]=SEIRD_VBSA(data,pDomain,QoIs,plots)
    arguments
        data,pDomain,QoIs {mustBeText} = {'allStates'},
        plots double {mustBeInteger} = 1,
        
    end
    
    tRange=data.tRange;
    nStates=data.nStates;
    y0=data.y0;
    tLoc=data.tLoc;
    z=data.z;
    
    if tLoc>max(tRange) || tLoc<0
        error("*** Error: tLoc must be in the considered time interval.");
    elseif z<0 || z>1
        error("*** Error: z must be in [0,1]");
    else
    
        addpath("C:\Users\meces\Desktop\COMPUTATIONAL STATISTICS\Project\safe_R1.1\safe_R1.1\VBSA");
        addpath("C:\Users\meces\Desktop\COMPUTATIONAL STATISTICS\Project\safe_R1.1\safe_R1.1\sampling");

        disp(" --- VBSA sensitivity analysis (Sobol idx) --- ");
        disp('  QOIs:');
        for i=1:length(QoIs)
            disp(['     ',QoIs{i}]);
        end
        disp(" ");

        nP=size(pDomain,2);
        nMC=100;

        % sample uniformly
        SEIRDmc(nMC,nP,pDomain,data,QoIs,0,0);

        %% Step 1: params distributions
        disp(" > STEP 1: params definition");

        DistrFun='unif';
        xmin=pDomain(1,:);
        xmax=pDomain(2,:);
        p_labels={'beta','r','d','i'};
        s_labels={'S','E','I','R','D'};

        warmup=1;% time from which indices are computed
        GSA_met='VBSA';
        %% Step 2: params sampling (build matrices X)
        disp(" > STEP 2: params sampling (build matrices X)");

        N=10000; % sample size
        SampStrategy='lhs'; % Sampling strategy

        DistrPar=cell(nP,1);
        for i=1:nP
            DistrPar{i}=[xmin(i) xmax(i)];
        end

        X=AAT_sampling(SampStrategy,nP,DistrFun,DistrPar,2*N); % Sampling
        [XA,XB,XC]=vbsa_resampling(X) ;

        %% Step 3: model evaluation (on X, to build matrices Y)
        disp(" > STEP 3: model evaluation (on X, to build matrices Y)");

        for j=1:N
            ySol=SEIRDsolver(XA(j,:),data);

            for i=1:length(s_labels)
                S=s_labels{i};
                eval(['YA_',S,'(j,:)=transpose(ySol(:,i));']);
            end
        end

        disp("       .Matrices YA: evaluated");

        for j=1:N
            ySol=SEIRDsolver(XB(j,:),data);

            for i=1:length(s_labels)
                S=s_labels{i};
                eval(['YB_',S,'(j,:)=transpose(ySol(:,i));']);
            end
        end

        disp("       .Matrices YB: evaluated");

        for j=1:nP*N
            ySol=SEIRDsolver(XC(j,:),data);

            for i=1:length(s_labels)
                S=s_labels{i};
                eval(['YC_',S,'(j,:)=transpose(ySol(:,i));']);
            end
        end

        disp("       .Matrices YC: evaluated");

        %% Step 4: compute time - varying Sensitivity Indices
        disp(" > Step 4: compute time - varying Sensitivity Indices (first order)");

        T=size(tRange,2);
        for i=1:length(s_labels)
            S=s_labels{i};
            eval(['Si_',S,'=nan(T,nP);']);
            eval(['STi_',S,'=nan(T,nP);']);
        end
        
        for t= warmup:T
            for i=1:length(s_labels)
                S=s_labels{i};
                eval(['[Si_',S,'(t,:),STi_',S,'(t,:)]=vbsa_indices(YA_',S,'(:,t),YB_',S,'(:,t),YC_',S,'(:,t));']);
            end
        end

        for j=1:length(p_labels)
            P=p_labels{j};
            for i=1:length(s_labels)
                S=s_labels{i};
                eval(['Si.',P,'(:,i)=Si_',S,'(:,j);']);
                eval(['STi.',P,'(:,i)=STi_',S,'(:,j);']);
            end
        end

        %% Step 5: plots
        disp(" > Step 5: plots");

        if plots>0
            lw='LineWidth';

            for i=1:length(s_labels)
                S=s_labels{i};
                Title1='Sobol indices ';
                Title2='Sobol Total indices ';

                for j=1:nP
                    figure (2)
                    hold on
                    eval(['plot(tRange,Si_',S,'(:,j),lw,2);']);
                    
                    figure (3)
                    hold on
                    eval(['plot(tRange,STi_',S,'(:,j),lw,2);']);
                end
                
                figure(2)
                legend ('beta','r','d','i');
                eval(['title(Title1,S);']);
                figure(3)
                legend ('beta','r','d','i');
                eval(['title(Title2,S);']);
            end

            
        end

%% End

        disp(" ");
        disp(" --- End VBSA sensitivity analysis (Sobol idx) --- ");
        disp(" ");
        disp(" ");
    end
end