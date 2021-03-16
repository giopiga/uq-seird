function paramLSQ=SEIRDlsq(param0,data,d)
    arguments
        param0,data,d double {mustBeInteger} = 1, 
    end
    
    idxSynt=data.idxSynt;

    [paramLSQ,ssLSQ]=fminsearch(@SEIRDss,param0,[],data);
   
    if d>0
        disp("        LSQ estimates for params:");
        fprintf('           .beta = %.4g\n',paramLSQ(1));
        fprintf('           .  r  = %.4g\n',paramLSQ(2));
        fprintf('           .  d  = %.4g\n',paramLSQ(3));
        fprintf('           .  i  = %.4g\n',paramLSQ(4));
    end

end