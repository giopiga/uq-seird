function PDFplot(X,F,names,Legend)
    n=size(X,1);

    for j=1:n
        
        x=X(j,:);
        for k=1:size(F,3)
            f=F(j,:,k);
            subplot(1,n,j)
            plot(x(~isnan(x)),f(~isnan(f)),'--','LineWidth',2);
            hold on 
        end
        title(names(j));
        if strcmp(Legend,"")==0
            legend(Legend);
        end
        grid on
        
    end
    sgtitle('Peak PDF');
    set(gca ,'FontSize',10);
%     set (gca ,'ColorOrderIndex' ,1);

end