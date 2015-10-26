v = -4;
b = 15;
a = [v, v+b, v-b];
if ( a(1)>=a(2) && a(1)>=a(3))
    amax = a(1)
else
    if a(2)>=a(3)
        amax = a(2)
    else
        amax = a(3)
    end
end

if v >= 0
    absv = v;
else
    absv = -v;
end

amax2 = absv + b
