function[a] = reassigncells(a,b,c)
    a(ismember(a,b)) = c;
end