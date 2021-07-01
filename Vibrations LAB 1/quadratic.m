function [roots] = quadratic(A,B,C)
    % Finding the roots using quadratic formula
    % roots = (-B +/- sqrt(-B^2 - 4ac)) / (2a)
    root1 = (-B + sqrt(B^2 - 4 * A * C)) / (2 * A);
    root2 = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A);
    roots = [root1 root2];
end

