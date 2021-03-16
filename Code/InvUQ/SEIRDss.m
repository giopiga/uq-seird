function [ss]=SEIRDss(param,data)
    idxSynt=data.idxSynt;
    
    y_mod=SEIRDsolver(param,data);
    
    y=y_mod(:,idxSynt);
    
    ss=sum(sum((data.ySol_synt-y).^2));
end