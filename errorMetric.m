function s = errorMetric()
syms rA1 rA2 rA3 rA4 thetaA1 thetaA2 thetaA3 thetaA4
syms rB1 rB2 rB3 rB4 thetaB1 thetaB2 thetaB3 thetaB4
syms a


totalDistSqrd = distSqrd(xCoord(rA1, thetaA1), yCoord(rA1, thetaA1), xCoord(rB1, thetaB1+a), yCoord(rB1, thetaB1+a));
totalDistPrime = diff(totalDistSqrd, a)
s = solve(totalDistSqrd, a, 'ReturnConditions', true)

end

function d = distSqrd(x1, y1, x2, y2)
d = (x1 - x2)^2 + (y1 - y2);
end

function x = xCoord(r, theta)
x = r * cos(theta);
end

function y = yCoord(r, theta)
y = r * sin(theta);
end