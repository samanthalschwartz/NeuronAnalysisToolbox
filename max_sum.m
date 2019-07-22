function [summa, index] = max_sum(v, n)
[row col] = size(v);
for c = 1:col
    summa(c) = sum(v(c:row*n));
end
index = n
end
