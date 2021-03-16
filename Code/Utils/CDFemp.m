function [X,Fu,Fc]=CDFemp(Yu,Yc)

    r=size(Yc,1);  % number of states of QoI
    cu=size(Yu,2);   % number of uncond. simulations (Nu)
    cc=size(Yc,2);   % number of cond. simulations (Nc)
    h=size(Yc,3);   % number of repetitions (n)
    
    C=sort(Yc,2);
    U=sort(Yu,2);
%     
    for i=1:r
        xMax(i)=max(max(C(i,end,:)),U(i,end));
        xMin(i)=min(min(C(i,1,:)),U(i,1));
        
        X(i,:)=linspace(xMin(i),xMax(i),floor(100/r));
    end
    
    
    for i=1:r
        for j=1:length(X(i,:))
            Fu(i,j)=sum(U(i,:)<=X(i,j));
            
            for k=1:h
                Fc(i,j,k)=sum(C(i,:,k)<=X(i,j));
            end
        end
    end
% 
    Fu=Fu/cu;
    Fc=Fc/cc;
end