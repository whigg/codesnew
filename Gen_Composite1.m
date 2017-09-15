function [Daily_V] = Gen_Composite(MODIS_t, ts, bad_ts, D_t); 
%%  Generate timeseries of DOY from S_num(start doy) and E_num (end doy)

% for i=1:length(MODIS_t)-1;  
%     I=find((D_t>=MODIS_t(i,1))&(D_t<MODIS_t(i+1,1))); 
%     Daily_V(I,1)=ts(i,1); 
%     Daily_V(I,2)=bad_ts(i,1); 
% end
% Daily_V(length(D_t),1)=ts(end,1); 
% Daily_V(length(D_t),2)=bad_ts(end,1); 
% 
% end



%%%%%%%  modified with new 8 days composite definition 
Y =find( diff(MODIS_t)>600);  % end date of year % this is when end of the year in the modis-t changes so [2011361;2012001] is in Y

for i=1:length(MODIS_t)-1  
    if any(i==Y)   
           I=find((D_t>=MODIS_t(i,1))&(D_t<MODIS_t(i+1,1)+3));      %  DOY 361 includes 361- 003 
               Daily_V(I,1)=ts(i,1); 
          %      Daily_V(I,2)=i; %bad_ts(i,1); 
    else
    I=find((D_t>=MODIS_t(i,1))&(D_t<MODIS_t(i+1,1))); 
    Daily_V(I,1)=ts(i,1); 
   %  Daily_V(I,2)=i;  %bad_ts(i,1); 
    end
    
end
Daily_V(length(D_t),1) = ts(end,1); 
%Daily_V(length(D_t),2) = 713; %bad_ts(end,1); 

end

