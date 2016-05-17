% Round x to n decimal places
function x2 = roundN(x, n)
f = 10^n;
x2 = round(x * f) / f;
end