function [GRLM10Raw, GRLM10Smooth, GRLM35Raw, GRLM35Smooth, LEGOS, DAHITI] = Allaltimetry_read (ss, ipath, ifiles); 

for k=1:length(ss)     
    if (k==1)&&(ss(k,1)==1) 
        for kk=1:2 % kk=1 for GRLM10 using "GRLM10." for search; kk=2 GRLM10Smooth using "GRLM10_" for search
             if kk==1 
                s=strfind({ifiles.name}, 'GRLM10.');  % always assumes that if GRLM exist in name, then we have both files GRLM10 and GRLM10_Smooth
                filename=[ipath char(ifiles(find(~cellfun(@isempty,s))).name)]; 
                [TPJO, GRLM10_Altimetry] = GRLM10_importfile(filename); % removes bad dates 
             else 
                 s=strfind({ifiles.name}, 'GRLM10_');  
                 filename=[ipath char(ifiles(find(~cellfun(@isempty,s))).name)]; 
                 GRLM10_Altimetry = GRLM10Smooth_importfile(filename); %   
             end
             GRLM10_Altimetry(:,1) =  yyyymmdd2doy(GRLM10_Altimetry(:,1)); 
             I=find(GRLM10_Altimetry(:,2)>900); % index nan altimetry data
             GRLM10_Altimetry(I,2)=NaN;   
                  
             clear T1 T2 T3 G1 G2 G3
             missions={'TOPEX', 'POSDN','Jason','OSTM'};
             t1= find(cellfun('length',regexp(TPJO,missions{1})) == 1); % index of all TOPEX mission data in GRLM10 
             t2= find(cellfun('length',regexp(TPJO,missions{2})) == 1); % index of all POSDN mission data in GRLM10 
             T1=union(t1,t2); %    index of all TOPEX and POSDN mission data in GRLM10 
             clear t1 t2 
             T2=find(cellfun('length',regexp(TPJO,missions{3})) == 1);  % index for Jason
             T3=find(cellfun('length',regexp(TPJO,missions{4})) == 1);  % index for OSTM
                 
             G1=GRLM10_Altimetry(T1,1:2); 
             G2=GRLM10_Altimetry(T2,1:2); 
             G3=GRLM10_Altimetry(T3,1:2); 
  
             if ~isempty(G1) && ~isempty(G2) && ~isempty(G3)  % each instruments were used during observed period at least once
                if (nansum(G1(:,2))>-9999)&&(nansum(G2(:,2))>-9999)&&(nansum(G3(:,2))>-9999)  % if we have more good data vs bad data?
                    %%%%  Processing GRLM10 to combine G1, G2 and G3 
                    [C1,ia,ib]=intersect(G1(:,1),G2(:,1)); % data for two instruments for the same date for (TOPEX or POSDN) and Jason
                    % all three are empty if no intersection - is it
                    % possible?
                    % ia is index of such dates in TOPEX/POSDN
                    % ib is index of such date in Jason
                    C1=find(abs(G1(ia,2)- G2(ib,2))==nanmin(abs(G1(ia,2)- G2(ib,2)))); % now C1 is overwritten with index of raws within G1, G2 and previous C1 where absolute difference between two instruments is the minimum
                    if isempty(C1)  
                        C1=ib(end,1); % use latest Jason and (Topex or posdn) data as C1
                    end 
                    if(ib(end,1)~=length(ib))  
                        C1=length(ib); % use latest Jason and (Topex or posdn) data as C1
                    end 
                    [C2,ic,id]=intersect(G2(:,1),G3(:,1)); 
                    C2=find(abs(G2(ic,2)- G3(id,2))==nanmin(abs(G2(ic,2)- G3(id,2)))); 
                    if isempty(C2)  
                        C2=1; % use first Jason and OSTM data as C1 on;y if minimum absolute difference not found not found
                    end 
                    GRLM10_n=[G1(1:ia(C1,1),1:2); G2(ib(C1,1)+1: ic(C2,1), 1:2); G3( id(C2,1)+1:end,1:2)]; % this is where duplicate dates because of different instruments are removed
                    % G1(1:ia(C1,1),1:2) - first we take TOPEX/POSDN data from start to the
                    % first date when TOPEX/POSDN altimentry data is
                    % closest to jason
                    % G2(ib(C1,1)+1: ic(C2,1), 1:2) - jason data from one day after
                    % topex/posdn alt data matches jason the best till when
                    % jason matches ostm
                    % G3( id(C2,1)+1:end,1:2) - ostm data from the best match between jason and ostm untill the end of the file 
                else  
                    GRLM10_n(:,1)=G3(:,1); % use only ostm data if there is no good data for topex/posdn and jason and ostm
                    GRLM10_n(:,2)=G3(:,2);  
                end
             else 
                 GRLM10_n(:,1)=GRLM10_Altimetry(:,1); % otherwise do not filter for duplicates
                 GRLM10_n(:,2)=GRLM10_Altimetry(:,2); 
             end

             if kk==1 
                 GRLM10Raw=GRLM10_n;  
             else 
                 GRLM10Smooth=GRLM10_n;  
             end
             clear GRLM10_Altimetry I startRow endRow filename D GRLM10_n GRLM10_int G1 G2 G3 GRLM10 
        end
     elseif (k==1)&&(ss(k,1)==0) % files with "GRLM10" in the name do not exist for this current lake
         GRLM10Raw=[];
         GRLM10Smooth=[]; 
     end

          
%%%%%       2)  Import GRLM35      %%%%%%%%
     if((k==2)&(ss(k,1)==1)); 
               for kk=1:2;   
                      if(kk==1); 
                            s=strfind({ifiles.name}, 'GRLM35.');  
                            filename=[ipath char(ifiles(find(~cellfun(@isempty,s))).name)]; 
                            startRow = 17;  endRow = inf;
                            GRLM35_Altimetry = GRLM35_importfile(filename, startRow, endRow); % 
                      else; 
                             s=strfind({ifiles.name}, 'GRLM35_');  
                             filename=[ipath char(ifiles(find(~cellfun(@isempty,s))).name)]; 
                             startRow = 12;  endRow = inf;
                             GRLM35_Altimetry = GRLM35Smooth_importfile(filename, startRow, endRow); % 
                      end
                 
                       I=find(GRLM35_Altimetry(:,1)<99999990);  
                       GRLM35=GRLM35_Altimetry(I,[1,4]);   
                       I=find((GRLM35(:,2)>900)); GRLM35(I,2)=NaN; 
                       GRLM35(:,1)= yyyymmdd2doy(GRLM35(:,1)); 

%                       GRLM35_int(:,1)=Gen_Timeseries(GRLM35(1,1), GRLM35(end,1));
%                      GRLM35_int(:,2)=interp1(GRLM35(:,1),GRLM35(:,2)-nanmean(GRLM35(:,2)), GRLM35_int(:,1)); 
 
                 if(kk==1); GRLM35Raw=GRLM35;   else; GRLM35Smooth=GRLM35;  end
                  clear GRLM35_Altimetry I startRow endRow filename  GRLM35 GRLM35_int  
               end
     elseif((k==2)&(ss(k,1)==0));  GRLM35Raw=[]; GRLM35Smooth=[]; 
     end

%%%%%      3)  Import LEGOS       %%%%%%%%
    if (k==3)&(ss(k,1)==1)
        s=strfind({ifiles.name}, 'LEGOS');  
        filename=[ipath char(ifiles(find(~cellfun(@isempty,s))).name)]; 
        %%%% new legos
        startRow = 1;
        endRow = inf;
        LEGOS_Altimetry = LEGOS_importfile(filename, startRow, endRow); % 
        I=find((LEGOS_Altimetry(:,2)>1000)&(LEGOS_Altimetry(:,3)>1));  
        LEGOS_Altimetry=LEGOS_Altimetry(I(1,1):end,:); % jj=85 error dimentions compare with jj=81
        LEGOS(:,1)=ConvertFOYtoDOY(LEGOS_Altimetry(:,1)); 
        LEGOS(:,2:3)=LEGOS_Altimetry(:,4:5);  
        %%%  old legos 
        % startRow = 7;endRow = inf;
        %  LEGOS_Altimetry = oldLEGOS_importfile(filename, startRow, endRow); % 
        % LEGOS(:,1)=ConvertFOYtoDOY(LEGOS_Altimetry(:,1)); 
        % LEGOS(:,2:3)=LEGOS_Altimetry(:,2:3);  

        %%%%%%%%%% remove the duplicate  
        [b,order]=unique(LEGOS(:,1),'first'); 
        v=LEGOS(order,:);  clear LEGOS; LEGOS=v; 
        % LEGOS_int(:,1)=Gen_Timeseries(LEGOS(1,1), LEGOS(end,1));
        % LEGOS_int(:,2)=interp1(LEGOS(:,1),LEGOS(:,2)-nanmean(LEGOS(:,2)), LEGOS_int(:,1)); 
        clear LEGOS_Altimetry startRow endRow  filename vector percntiles outlierIndex
    elseif((k==3)&(ss(k,1)==0));  
        LEGOS=[]; 
    end
 
%%%%%      4)  Import DAHITI       %%%%%%%%
 if((k==4)&(ss(k,1)==1)); 
s=strfind({ifiles.name}, 'DAHITI');  
filename=[ipath char(ifiles(find(~cellfun(@isempty,s))).name)]; 
 startRow = 15;  endRow = inf;
DAHITI_Altimetry = DAHITI_RTimportfile(filename, startRow, endRow);  % date: cell,  
%%%% convert cell to number 
D=  yyyymmdd(datetime(datenum(DAHITI_Altimetry(:,1)),'ConvertFrom','datenum'));  
for i=1:length(D); DAHITI(i,1) =  yyyymmdd2doy(D(i,1)); end
%DAHITI(:,2)=DAHITI_RTimportfile(filename, startRow, endRow)
DAHITI(:,2) = str2double(DAHITI_Altimetry(:,2));
clear DAHITI_Altimetry endRow startRow filename D
% DAHITI_int(:,1)=Gen_Timeseries(DAHITI(1,1), DAHITI(end,1));
% DAHITI_int(:,2)=interp1(DAHITI(:,1),DAHITI(:,2)-nanmean(DAHITI(:,2)), DAHITI_int(:,1)); 
 elseif((k==4)&(ss(k,1)==0));  DAHITI=[]; 
 end
 
end



