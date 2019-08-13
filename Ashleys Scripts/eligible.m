function admit = eligible(v,q)
M = [v q];
if v >= 88 & q >= 88 & mean(M)>92
admit = true
else
    admit = false
end