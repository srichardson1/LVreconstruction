time=zeros(65536,1);
aortic=zeros(65536,1);
external=zeros(65536,1);


for j=2:65536
   time(j,1)=time(j-1,1)+(0.8/65537); 
end


for j=1:65536
    
    if time(j,1)<0.2
    external(j,1)=13*(time(j,1)/0.2);
    aortic(j,1)=((1730)*((time(j,1)+0.6)^3))-((3625)*((time(j,1)+0.6)^2))+((2454)*(time(j,1)+0.6))-(440);
    end
    
    if time(j,1)>=0.2 && time(j,1)<0.5
    external(j,1)=13;
    end

    if time(j,1)>=0.5 && time(j,1)<0.7
    external(j,1)=13+((140)*(1-(exp(-(((time(j,1)-0.5)*(time(j,1)-0.5))/(0.007))))));
    end
    
    if time(j,1)>=0.7
    external(j,1)=152.53*(1-(exp(-(((0.8-time(j,1))*(0.8-time(j,1)))/(0.0015)))));
    end    
    
    
    if time(j,1)>=0.2 && time(j,1)<0.575
    aortic(j,1)=-((48.68)*(time(j,1)-0.2))+87.75;
    end
    

    if time(j,1)>=0.575 && time(j,1)<0.75
    aortic(j,1)=((5800)*((time(j,1)-0.2)*(time(j,1)-0.2)*(time(j,1)-0.2)))-((13288)*((time(j,1)-0.2)*(time(j,1)-0.2)))+((8726)*((time(j,1)-0.2)))-(1634);
    end
    
    
    if time(j,1)>=0.75
    aortic(j,1)=((1730)*((time(j,1)-0.2)*(time(j,1)-0.2)*(time(j,1)-0.2)))-((3625)*((time(j,1)-0.2)*(time(j,1)-0.2)))+((2454)*((time(j,1)-0.2)))-(440);
    end
    
end



