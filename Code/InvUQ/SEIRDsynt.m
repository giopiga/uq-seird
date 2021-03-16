function [ySol_synt,idxSynt,pV,pT]=SEIRDsynt(param,data,statesSynt,plots)
    arguments
        param,data,statesSynt,plots double {mustBeInteger} = 1, 
    end
    
    ySol_synt0=SEIRDsolver(param,data);
    
    sigma2=data.sigma2;
    
    states={'S','E','I','R','D'};
    for j=1:length(statesSynt)
        [v,idxSynt(j)]=max(strcmp(states,statesSynt{j}));
    end
    
    [pV,pT]=max(ySol_synt0(:,idxSynt(1)));
    
    r=size(ySol_synt0,1);
    ySol_synt=ySol_synt0(:,idxSynt)+normrnd(0,sigma2,r,1);
    
%     i=1;
%     for j=1:length(ySol_synt0)
%         sigma2=0.1*ySol_synt0(j,idxSynt(i));
%         ySol_synt(j,i)=ySol_synt0(j,idxSynt(i))+normrnd(0,sigma2,1);
%     end
%     
%     i=2;
%     for j=2:length(ySol_synt0)
%         sigma2=8*(ySol_synt0(j,idxSynt(i))-ySol_synt0(j-1,idxSynt(i)));
%         ySol_synt(j,i)=ySol_synt0(j,idxSynt(i))+normrnd(0,sigma2,1);
%     end
    
    if plots>0
        figure(200)
        
        hold on
        for i=1:size(ySol_synt,2)
            plot(data.tRange,ySol_synt(:,i),'.','LineWidth',3);
        end 
        
        plot(data.tRange,ySol_synt0(:,idxSynt),'k-','LineWidth',2);
           
        legend(statesSynt);
        title('Synthetic data');
        grid on
        set(gca,'FontSize',10);
    end    
    
end