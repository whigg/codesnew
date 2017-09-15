% Converts date YYYYMMDD format 20170817 to DOY 2017229
function [doy] = yyyymmdd2doy(yyyymmdd)

if (size(yyyymmdd)==1) 
%--year
year=floor(yyyymmdd/10000);
%--month
month=floor((yyyymmdd-10000*year)/100);
%--day
day=yyyymmdd-year*10000-month*100;

%-- doy of the first day for each month 
if(isleap(year)==1)
    idoy=[1 32 61 92 122 153 183 214 245 275 306 336];
else
    idoy=[1 32 60 91 121 152 182 213 244 274 305 335];
end

doy=(year*1000)+idoy(month)-1+day;
else 
    
    for i=1:length(yyyymmdd) 
        year=floor(yyyymmdd(i,1)/10000);
%--month
month=floor((yyyymmdd(i,1)-10000*year)/100);
%--day
day=yyyymmdd(i,1)-year*10000-month*100;

%-- doy of the first day for each month 
if(isleap(year)==1)
    idoy=[1 32 61 92 122 153 183 214 245 275 306 336];
else
    idoy=[1 32 60 91 121 152 182 213 244 274 305 335];
end

doy(i,1)=(year*1000)+idoy(month)-1+day;
    end
    



end