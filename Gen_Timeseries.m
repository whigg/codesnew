function [T_series] = Gen_Timeseries(S_num, E_num)
%%  Generate timeseries of DOY from S_num(start doy) and E_num (end doy)


S_num=1979243;
E_num=2017365; 
  
V_int(:,1)=[S_num:1:E_num]'; 

year=floor(V_int(:,1)*0.001); % very bad!
doy=V_int-(year*1000); 
  % if (isleap(year)==1); dayInYear=366;  else; dayInYear=365; end
   
 for i=1:length(V_int); 
      
    if isleap(year(i,1))==1
        dayInYear=366;  
    else
        dayInYear=365;
    end
      
    if doy(i,1)<=dayInYear
        n(i,1)=1;
    else
        n(i,1)=0;
    end
     % if((i/100)==floor(i/100)); n(i,1)=0; end
    if(doy(i,1)==0)
        n(i,1)=0; 
    end
 end
% for k=1:floor(length(n)/1000); n(k*1000,1)=0; end; ??? 
 
 I=find(n==1);  
 T_series=year(I,1)*1000+ doy(I,1); 

end