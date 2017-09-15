function [num] = ConvertDOYtoFOY( yy )
%%%%%%  Convert from DOY (2009001) to Fraction of Year (2009.0027) 
y=yy*0.001; 
%  year = floor(y);
%   partialYear = mod(y,1);
%   date0 = datenum(num2str(year),'yyyy');
%   date1 = datenum(num2str(year+1),'yyyy');
%   daysInYear = date1 - date0;
%   num =year+  ((partialYear * daysInYear)*0.01);  

 year=floor(y); 
 
 if (isleap(year)==1); dayInYear=366; 
else; dayInYear=365; 
 end
 
 partialYear=(yy-(year*1000));
 num=year+((partialYear/dayInYear)); 
end


