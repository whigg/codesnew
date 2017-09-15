% this is the main code for caculation of correlaion between altimetry
% height and surface area, and lakes storage change.
% clear all

%%% 1) Correlations %%%

% In Command Window, type Target_ID = 0; to copy in lakes IDs.
clearvars -except Target_ID

%%%%%%%%%%   update on 9/14/2017  %%%%%%%%%%%%%

% load MODIS dates, 'modis_dates.mat' 
load('D:\aqua\codesnew\modis_dates.mat');

% make a matrix DV
DV = datevec(dates);

% Format dates to DOY
for i = 1:length(DV)
    MODIS_t(i,1) = yyyymmdd2doy(DV(i,1)*10000+DV(i,2)*100+DV(i,3)); 
end
 

% Read Area data to Target_ID
for jj=1:length(Target_ID)
    jj
    
    clearvars -except DataExist Targets Target_ID  MODIS_t jj Corr    Max_Data Min_Data Num_Data Corr_1 Corr_2 Corr_org Obs RR   
    % Read area data from folder in 'D:\...' 
    ipath = 'D:\aqua\area\';
    
    % look for Surface Area values in all 'ts' files in columns 2; if there is an
    % empty value, look for in if column 3.
    ifiles=dir([ipath 'ts_','*','_', num2str(Target_ID(jj,3)),'_','*']); 
    if length(ifiles)==0 
         ifiles=dir([ipath 'ts_','*','_', num2str(Target_ID(jj,2)),'_','*']);
    end
    number_of_area_files = length(ifiles);
    
    if number_of_area_files>0
    
        %%%%%%   Extract Surface Area data  (ts: timeseries of surface area) 
        for k = 1:length(ifiles) 
            filename = [ipath char(ifiles(k,1).name)]; 
            load(filename);    
            I = find(ts_water<0); 
            ts_water(I) = NaN; 
            t(:,k) = ts_water'; % t_bad(:,k)=ts_bad';  max_v(1,k)=max_val; 
        end
        
        if k>=2
            ts = t(:,1)+t(:,2);  
        else
            ts = t; 
        end 
    end         
    
    %%%%%% Read and Extract the Altimetry Data %%%%%%%                
    
    ipath = 'D:\aqua\alt\';
    
    % look for files in Target_ID in column 1
    ifiles = dir([ipath num2str(Target_ID(jj,1)),'_*']);
      
    % iterate  
    TYPE = {'GRLM10', 'GRLM35','LEGOS', 'DAHITI','ESA'};    
    for type_iterator = 1:5 
        s = strfind({ifiles.name}, TYPE{1,type_iterator});
        if sum(~cellfun(@isempty,s))==0 
            ss(type_iterator,1) = 0; 
        else  
            ss(type_iterator,1) = 1; 
        end
        DataExist(jj, type_iterator) = ss(type_iterator,1);
    end
    % Create a table with indicators for existing data source
    DataExist(jj, length(TYPE)+1) = sum(DataExist(jj,1:5)); % Altimetry Availability
    DataExist(jj, length(TYPE)+2) = number_of_area_files;
    I = find(ss>0);

    [GRLM10Raw, GRLM10Smooth, GRLM35Raw, GRLM35Smooth, LEGOS, DAHITI] = Allaltimetry_read (ss, ipath, ifiles); 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    if sum(ss)>0
        for kk = 1:6     
            clear Data

            if (kk==1)&&(ss(1,1)==1)
                Data(:,1) = GRLM10Raw(:,1); 
                Data(:,2) = GRLM10Raw(:,2)-nanmean(GRLM10Raw(:,2)); 
            end 

%             if (kk==2)&&(ss(1,1)==1)
%                 Data(:,1) = GRLM10Smooth(:,1);
%                 Data(:,2) = GRLM10Smooth(:,2)-nanmean(GRLM10Smooth(:,2));
%             end 

%             if (kk==3)&&(ss(2,1)==1)
%                 Data(:,1)=GRLM35Raw(:,1); 
%                 Data(:,2)=GRLM35Raw(:,2)-nanmean(GRLM35Raw(:,2)); 
%             end     

%             if (kk==4)&&(ss(2,1)==1)
%                 Data(:,1) = GRLM35Smooth(:,1); 
%                 Data(:,2) = GRLM35Smooth(:,2)-nanmean(GRLM35Smooth(:,2)); 
%             end  

%             if (kk==5)&&(ss(3,1)==1) 
%                 Data(:,1) = LEGOS(:,1); 
%                 Data(:,2) = LEGOS(:,2)-nanmean(LEGOS(:,2)); 
%             end

%             if (kk==6)&&(ss(4,1)==1) 
%                 Data(:,1) = DAHITI(:,1); 
%                 Data(:,2) = DAHITI(:,2)-nanmean(DAHITI(:,2)); 
%             end 

  
            if exist('Data','var')&& exist('ts','var')
                D_t = Gen_Timeseries(MODIS_t(1,1), MODIS_t(end,1));  % overwrites with S_num=1979243; E_num=2017365; in function
                D_t(:,2) = Gen_Composite1(MODIS_t, ts, ts, D_t(:,1));  
               
                % remove duplicates:
                index_duplicates = find(hist(Data(:,1),unique(Data(:,1)))>1);
                DD = Data;
                DD(index_duplicates,:) = [];
                C = intersect(DD(:,1),D_t(:,1));  
                I = find(ismember(DD(:,1),C));  
                II= find(ismember(D_t(:,1), C));
                
                %%%  Calculate the correlation 
                Area = D_t(II,1:2);
                Altimetry = DD(I,1:2);                
                Area_and_altimetry = [Area(:,2), Altimetry(:,2)];
                
                I = find(Area_and_altimetry(:,1)>0);
                % Corr(jj,kk) = nancorr(D_t(II,2),Data(I,2));
                Corr(jj,kk) = nancorr(Area_and_altimetry(I,1),Area_and_altimetry(I,2));
            end
        end
    end
end
 
%%% 2) Calculate Storage change %%%

 VolData=[C Area(:,2) Altimetry(:,2)];   % C: time, Area: Surface Area Timeseries, Depth: Altimetry height timeseries
   Iv=find((VolData(:,2)>0)&(VolData(:,3)>-100));  
   VolData=VolData(Iv,:);  
   VolData(:,2)=VolData(:,2)*463.3*463.3;   % Convert to meter
  
 
   Sort_VolData=  sortrows(VolData,2); %   Time, Area, Altimetry - sorted by AREA:
     
     diffVol(1,1)=(0.5*(Sort_VolData(1,3)).*Sort_VolData(1,2))  ; 
     for k=2:length(Sort_VolData); 
         diffArea(k,1)=Sort_VolData(k,2)-Sort_VolData(k-1,2); 
         diffAlt(k,1)=Sort_VolData(k,3)-Sort_VolData(k-1,3); 
         diffVol(k,1)=diffVol(k-1,1)+ ( (diffAlt(k,1)) *  (Sort_VolData(k-1,2)+0.5*diffArea(k,1)) ); 
     end
     
     % Table: Time   Area(sorted)   Altimetry   diffArea   diffAltimetry
     % diff Volume
 DiffData = [Sort_VolData diffArea diffAlt diffVol];
     % Volume=Sort_Data(:,3).*Sort_Data(:,2)./2000;
     DDv=[Sort_VolData(:,1) diffVol(:,1)]; 
     DDv2=sortrows(DDv,1); 
         
 Vol_est=[VolData(:,1) DDv2(:,2)-min(DDv2(:,2))];    % 1st colume: Time, 2nd colume:  relative volume change (m^3) 

 % write output in text file 
 

% % % dlmwrite('D:\aqua\results\Hypsometry_509_GRLM10.txt',VolData,'delimiter','\t','precision',3);
% % % dlmwrite('D:\aqua\results\DiffData_509_GRLM10.txt',VolData,'delimiter','\t','precision',3);
% % % dlmwrite('D:\aqua\results\VolumeEst_509_GRLM10.txt',Vol_est,'delimiter','\t','precision',3);
% % % fig_Corr=scatter((Area_and_altimetry(:,1)),(Area_and_altimetry(:,2))); xlabel ('Area'); ylabel ('Altimetry');
% % % saveas(gca,'D:\aqua\results\fig_Corr_345_GRLM10','fig');
% % % saveas(gca,'D:\aqua\results\fig_Corr_345_GRLM10','jpeg');
% % % fig_Vol=plot((VolData(:,1)),DDv2(:,2)-min(DDv2(:,2)));xlabel ('Time'); ylabel ('Est Volume');
% % % saveas(gca,'D:\aqua\results\fig_Vol_345_GRLM10','fig');
% % % saveas(gca,'D:\aqua\results\fig_Vol_345_GRLM10','jpeg');

% xlswrite('D:\aqua\results\res1.xlsx',Vol_est,'Sheet3','A1');
% xlswrite('D:\aqua\results\res1.xlsx',Corr,'Sheet1','A2');
% xlswrite('D:\aqua\results\res1.xlsx',C,'Sheet2','A2');
% xlswrite('D:\aqua\results\res1.xlsx',Area_and_altimetry,'Sheet2','B2');
