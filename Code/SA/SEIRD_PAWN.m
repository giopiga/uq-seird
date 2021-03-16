function SEIRD_PAWN(data,pDomain,tcdf,QoIs,plots)
    arguments
        data,pDomain,tcdf,QoIs {mustBeText} = {'allStates'},
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
        
        disp(" --- PAWN sensitivity analysis ---");
        disp('  QOIs:');
        for i=1:length(QoIs)
            disp(['     ',QoIs{i}]);
        end
        disp(" ");

        nStates=length(y0);

        data.y0=y0;
        data.tRange=tRange;
        data.nStates=nStates;

        nP=size(pDomain,2);
        nMC=100;

        % sample uniformly
        SEIRDmc(nMC,nP,pDomain,data,QoIs,0,0);

%% Step 1: params distributions
        disp("   > STEP 1: params definition");
        N=1e4;

        xmin=pDomain(1,:);
        xmax=pDomain(2,:);
        nP=length(xmin); % number of inputs
        p_labels={'beta','r','d','i'};
        s_labels={'S','E','I','R','D'};
        tcdf=tcdf(tcdf<=max(tRange));        % avoiding tcdf out of range

        z=0;    

%% Step 2: unconditined CDF
        disp("   > STEP 2: evaluation of unconditioned CDFs");

        N=1000; % number of samples to estimate unconditional CDF

        [NySol]=SEIRDmc(N,nP,pDomain,data,QoIs,0);

        % restriction to the istant we want to consider (tcdf)
        for j=1:length(tcdf)
            T=num2str(tcdf(j));
            eval(['NySol_',T,'=squeeze(NySol(',T,',:,:));']);
        end

        % evaluation of empirical unconditioned cdf
        for i=1:length(tcdf)
            T=num2str(tcdf(i));
            for j=1:length(s_labels)
                S=s_labels{j};

                eval(['[cdf_',S,'_',T,',X_',S,'_',T,']=ecdf(NySol_',T,'(',num2str(j),',:));']);
                eval(['Fu(:,j,i)=cdf_',S,'_',T,';']);
                eval(['Xu(:,j,i)=X_',S,'_',T,';']);
            end
        end

%% Step 3: conditined CDF
        disp("   > STEP 3: evaluation of conditioned CDFs");

        Nc=100 ; % number of samples to estimate conditional CDFs
        n=10 ; % number of conditioning points

        % sampling n values to fix each parameter (one at time)
        for j=1:nP
           eval(['p',num2str(j),'=(pDomain(2,',num2str(j),')-pDomain(1,',num2str(j),'))*rand(n,1)+pDomain(1,',num2str(j),');']);
        end

        unit_set=rand(Nc,nP);
        p0(:,1)=pDomain(1,1)+(pDomain(2,1)-pDomain(1,1))*unit_set(:,1);
        p0(:,2)=pDomain(1,2)+(pDomain(2,2)-pDomain(1,2))*unit_set(:,2);
        p0(:,3)=pDomain(1,3)+(pDomain(2,3)-pDomain(1,3))*unit_set(:,3);
        p0(:,4)=pDomain(1,4)+(pDomain(2,4)-pDomain(1,4))*unit_set(:,4);

        % evaluations for each fixed param
        for i=1:n
            for j=1:Nc
                p_condBeta=[p1(i),p0(j,2),p0(j,3),p0(j,4)];
                p_condR=[p0(j,1),p2(i),p0(j,3),p0(j,4)];
                p_condD=[p0(j,1),p0(j,2),p3(i),p0(j,4)];
                p_condI=[p0(j,1),p0(j,2),p0(j,3),p4(i)];

                ySolc_beta(:,:,j,i)=SEIRDsolver(p_condBeta,data);
                ySolc_r(:,:,j,i)=SEIRDsolver(p_condR,data);
                ySolc_d(:,:,j,i)=SEIRDsolver(p_condD,data);
                ySolc_i(:,:,j,i)=SEIRDsolver(p_condI,data);
            end
        end

%% Step 4: plots
        if plots>0
            disp("   > STEP 4: plots");

            for i=1:nP
                I=p_labels{i};
                for j=1:n   

                    % evaluation of empirical conditioned cdf
                    % cdf_S_T_I stands for cdf of S at time T conditioned to param I
                    for k=1:length(tcdf)
                        T=num2str(tcdf(k));
                        eval(['ySolc_',T,'_',I,'=squeeze(ySolc_',I,'(',T,',:,:,j));']);

                        for w=1:length(s_labels)
                            S=s_labels{w};

                            eval(['[cdf_',S,'_',T,'_',I,',X_',S,'_',T,'_',I,']=ecdf(ySolc_',T,'_',I,'(',num2str(w),',:));']);
                            eval(['Fc(:,w)=cdf_',S,'_',T,'_',I,';']);
                            eval(['Xc(:,w)=X_',S,'_',T,'_',I,';']);
                        end

                        % plot of conditioned cdfs
                        Title='';
                        col='k';
                        LW=1;

                        figure(length(tcdf)*(i-1)+k+3)
                        CDFplot(Xc,Fc,col,LW,Title,s_labels);
                    end
                end
            end

            % plot of unconditioned cdfs
            col='r';
            LW=2;
            for i=1:nP
                I=p_labels{i};

                for k=1:length(tcdf)
                    T=num2str(tcdf(k));
                    figure(length(tcdf)*(i-1)+k+3)
                    eval(['Title="cdf - t=',T,' : ',I,'";']);
                    
                    sXu=squeeze(Xu(:,:,k));
                    sFu=squeeze(Fu(:,:,k));
                    CDFplot(sXu,sFu,col,LW,Title,s_labels);
                end
            end
        
        else
            disp("   > STEP 4: plots - missed by user options (plots=0)");
        end

%% End
        disp(" ");
        disp(" --- End PAWN sensitivity analysis --- ");
        disp(" ");
        disp(" ");
    end
end
