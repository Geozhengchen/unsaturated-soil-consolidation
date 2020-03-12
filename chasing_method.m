function x = chasing_method(V,W,d)
global E J
aa = cell(J+1,1);
bb = cell(J+1,1);
cc = cell(J+1,1);
alpha = cell(J,1);
beta = cell(J+1,1);
y = cell(J+1,1);
% define aa, bb, and cc
bb{1} = E;
cc{1} = zeros(2,2);
for i = 2:J
    aa{i} = -V;
    bb{i} = W;
    cc{i} = -V;
end
aa{J+1} = -2*V;
bb{J+1} = W;
% solve alpha and beta
beta{J+1} = bb{J+1};
for i = J:-1:1
    alpha{i} = cc{i}*(beta{i+1})^(-1);
    beta{i} = bb{i}-alpha{i}*aa{i+1};
end
% solve y
y{J+1} = d{J+1};
for i = J:-1:1
    y{i} = d{i}-alpha{i}*y{i+1};
end

x{i} = ((beta{1})^(-1))*y{1};
for i = 2:J+1
    x{i} = ((beta{i})^(-1))*(y{i}-aa{i}*x{i-1});
end
