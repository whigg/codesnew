

% Load the data from Main_correlations.m 
% The input data - C and Area_and_altimetry, or Altimetry and Area,
% including C (time)

   VolData=[C Area(:,2) Altimetry(:,2)];   % C: time, Area: Surface Area Timeseries, Depth: Altimetry height timeseries
   Iv=find((VolData(:,2)>0)&(VolData(:,3)>-100));  
   VolData=VolData(Iv,:);  
   VolData(:,2)=VolData(:,2)*463.3*463.3;   % Convert to meter
   Sort_VolData=  sortrows(VolData,2);
     
     diffVol(1,1)=(0.5*(Sort_VolData(1,3)).*Sort_VolData(1,2))  ; 
     for k=2:length(Sort_VolData); 
         diffArea(k,1)=Sort_VolData(k,2)-Sort_VolData(k-1,2); 
         diffDepth(k,1)=Sort_VolData(k,3)-Sort_VolData(k-1,3); 
         diffVol(k,1)=diffVol(k-1,1)+ ( (diffDepth(k,1)) *  (Sort_VolData(k-1,2)+0.5*diffArea(k,1)) ); 
     end
    % Volume=Sort_Data(:,3).*Sort_Data(:,2)./2000; 
     DDv=[Sort_VolData(:,1) diffVol(:,1)]; 
     DDv2=sortrows(DDv,1); 
        
 Vol_est=[VolData(:,1) DDv2(:,2)-min(DDv2(:,2))];    % 1st colume: Time, 2nd colume:  relative volume change (m^3) 

 scatter (Data(:,1),DD2(:,2)-min(DD2(:,2)));
 
 % write output in text file 
%  ipath = 'D:\aqua\results'; 
 xlswrite('D:\aqua\results\vol.xlsx',Vol_est,'Sheet1','A1');
 xlswrite('D:\aqua\results\vol.xlsx',diffArea,'Sheet1','D1');
 xlswrite('D:\aqua\results\vol.xlsx',diffDepth,'Sheet1','E1');
 xlswrite('D:\aqua\results\vol.xlsx',diffVol,'Sheet1','F1');
%  text_file=([ipath 'Storage_',num2str(Target_ID(kk,1)),'.txt']);
%  dlmwrite(text_file,Vol_est,'delimiter',' ','precision','%.3f')
 
       