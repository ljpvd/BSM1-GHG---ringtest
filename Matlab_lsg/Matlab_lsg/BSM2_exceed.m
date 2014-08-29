% calculation of length and frequency of exceedances in BSM2
k = 1;
exc = [];
tag = 0;
for j = 2:lt
    if ((vec(j) > th) && (vec(j-1) <= th))
        k = j;
    else
        if  ((vec(j) <= th) && (vec(j-1) > th))
            exc = [exc (t(j) - t(k))];
            tag = 1;
        end
    end
end
if ((k == 1) && (vec(j) > th)) || ((k > 1) && (tag == 0))
    exc = [exc (t(j) - t(k))];
end
freq = length(exc);
sexc = sum(exc);
