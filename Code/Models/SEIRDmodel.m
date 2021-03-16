function F=SEIRDmodel(t,y,p)
    N_pop=1;

    beta=p(1)*(1-p(5));         %  Eventually, in p(5) is stored the lockdown param
    r=p(2);
    d=p(3);
    i=p(4);
    
    S=y(1);
    E=y(2);
    I=y(3);
    R=y(4);
    D=y(5);
    

    F=[ -beta/N_pop*S*I;
        beta/N_pop*S*I-i*E;
        i*E-d*I-r*I;
        r*I;
        d*I];
return