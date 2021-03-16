function [x,f,mu,sigma]=PDFemp(Y)
    [f,x]=ksdensity(Y,'Support','positive');
    mu=mean(Y);
    sigma=std(Y);
end