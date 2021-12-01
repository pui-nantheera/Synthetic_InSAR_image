function[year] = yyyymmdd2year(yyyymmdd)

year = floor(yyyymmdd/10000);
day = mod(yyyymmdd,100);
month = floor((yyyymmdd-year*10000)/100);

if mod(year,4)==0
  daylist=[31 60 91 121 152 181 213 244 274 305 335 366];
  totday = 366;
else
  daylist=[31 59 90 120 151 180 212 243 273 304 334 365];
  totday = 365;
end

doy=day;
if month>1
  doy=day+daylist(month-1);
end
year = year+doy/totday;
