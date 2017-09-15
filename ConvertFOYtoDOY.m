function [num] = ConvertFOYtoDOY( yy )
%%%%%%  Convert from  Fraction of Year (2009.0027) to DOY (2009001) 
%  year = floor(y);
%   partialYear = mod(y,1);
%   date0 = datenum(num2str(year),'yyyy');
%   date1 = datenum(num2str(year+1),'yyyy');
%   daysInYear = date1 - date0;
%   num =year+  ((partialYear * daysInYear)*0.01);  

 year=floor(yy);
 
%  if mod(year,400)==0; dayInYear=366; 
% elseif mod(year, 4) == 0 && mod(year, 100) ~= 0; dayInYear=366;
% else; dayInYear=365; 
% end

if (isleap(year)==1); dayInYear=366;  else; dayInYear=365; end

partialYear=(yy-(year));
 num=(year*1000)+((partialYear*dayInYear))+1; 
 num=round(num); 
end