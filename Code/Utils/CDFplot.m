function CDFplot(X,F,Col,LW,Title,states,KS)
    n=size(X,1);

    for j=1:n
        S=states{j};
        
        x=X(j,:);
        for k=1:size(F,3)
            f=F(j,:,k);
            subplot(1,n,j)
            plot(x(~isnan(x)),f(~isnan(f)),'col',Col,'LineWidth',LW);
            hold on
        end
        title(S);
        
        a=['KS=',num2str(KS(j))];
%         text(2*max(x)/3,0.05,a,'FontSize',20/n);
        xlabel(a);
    end
    sgtitle(Title);
end