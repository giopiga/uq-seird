function [ySol]=SEIRDsolver(param,data)
    
    tRange=data.tRange;
    y0=data.y0;
    tLoc=data.tLoc;
    z=data.z;
    
    param=[param,0];
    
    if tLoc>0 && tLoc<max(tRange)           % if lockdown is considered
        [v,idx]=max(tRange==tLoc);
        
        tRange1=tRange(1:idx);
        [~,y1]=ode45(@(t,y)SEIRDmodel(t,y,param),tRange1,y0);
        
        tRange2=tRange(idx:end);
        param(5)=z;
        y0=y1(end,:);
        [~,y2]=ode45(@(t,y)SEIRDmodel(t,y,param),tRange2,y0);
        
        y=[y1;y2(2:end,:)];                  % if lockdown is not considered
    else
        [~,y]=ode45(@(t,y)SEIRDmodel(t,y,param),tRange,y0 );
    end
    
    ySol=y;    
end