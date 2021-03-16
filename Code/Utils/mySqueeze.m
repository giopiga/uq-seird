function Sq=mySqueeze(M,d)
    l=size(M);
    
    k=1;
    L=zeros(1,length(l)-length(d));
    for j=1:length(l)
        if sum(d==j)==0
            L(k)=l(j);
            k=k+1;
        end
    end
   
    Sq=reshape(M,L);
end