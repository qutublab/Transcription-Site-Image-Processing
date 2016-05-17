function minutes = extractTimeInMinutes(s)
minutes = 0;

timeTags(1).tag = 'hrs_';
timeTags(1).minutes = 60;

timeTags(2).tag = 'hr_';
timeTags(2).minutes = 60;

timeTags(3).tag = 'mins_';
timeTags(3).minutes = 1;

timeTags(4).tag = 'min_';
timeTags(4).minutes = 1;

for i = 1:numel(timeTags)
    k = strfind(s, timeTags(i).tag);
    if numel(k) > 1
        error('[extractTime] Tag %s occurs %d times in%s', timeTags(i).tag, numel(k), s);
    else
        if numel(k) == 1
            tm = getNumberBefore(s, k);
            if tm < 0
                error('Unable to find time value in %s', s);
            end
            minutes = tm * timeTags(i).minutes;
            return;
        end
    end
end

end

% Returns the non-negative integer immediately preceeding position k in
% string s.  If no such number exists, -1 is returned.
function n = getNumberBefore(s, k)
n = 0;
mult = 1;
i = k - 1;
found = false;
while i > 0
   c = s(i);
   if c >= '0' && c <= '9'
       found = true;
       digit = c - '0';
       n = n + (digit * mult);
       mult = mult * 10;
       i = i - 1;
   else
       break;
   end
end
if ~found
    n = -1;
end
end