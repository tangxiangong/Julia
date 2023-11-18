clc; clear

n = 100000;
x = rand(n, 1);
tic; sin(x); toc
f = @(x) sin(x);
tic; map(f, x); toc

function y = map(f, x)
    n = length(x);
    y = zeros(n, 1);
    for k = 1:length(x)
        y(k) = f(x(k));
    end
end
