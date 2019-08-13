function coded = caesar(txt,n)
y = double(txt);
if n > 94
while n > 94
    n = (n-95);
end
elseif n < -94
while n < -95
    n = -(abs(n) - 95);
break
end
end
r = y + n
for a = 1:length(r)
if r(a) < 32
    r(a) = 126-(31-r(a))
elseif r(a) > 126
    r(a) = 31 + (r(a)-126)
else
end
end
coded = char(r)
end