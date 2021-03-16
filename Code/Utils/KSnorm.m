function KS=KSnorm(Fc,Fu,s_labels_QoIstates,P,disp)
    arguments
        Fc,Fu,s_labels_QoIstates,P,disp double {mustBeInteger} = 1,
    end
    
    format bank
    
    for s=1:length(s_labels_QoIstates)
        for k=1:size(Fc,3)
            maxErr(k)=max(abs(Fu(s,:)-Fc(s,:,k)));
        end
        
        x=mean(maxErr);
        KS(s)=round(x - sign(x)*.5/10^3,3);
    end
    
    if disp>0
        l=length(P);
        space=repelem(' ',1,21-l);

        disp([space,'<strong>',P,'</strong>       ',num2str(KS)]);
    end
end