function SEIRDplot(x,y,Legend,Title,meanPlot)
        arguments
             x,y,Legend = ''
             Title = ''
             meanPlot double {mustBeInteger} = 0;
        end
        
        %plot  MC realizations
        nMC=size(y,3);
        
        for j=1:nMC
            plot(x,y(:,:,j),'LineWidth',2);
            set (gca,'ColorOrderIndex',1)
            
            hold on
        end

        if meanPlot>0
            yMean=mean(y,3);
            % plot MC expected value
            hold on
            plot(x,yMean,'k-','LineWidth',4);
        end
        
        legend(Legend);
        title(Title);
        set(legend,'Location','East');
        grid on
        set(gca,'FontSize',10);
end