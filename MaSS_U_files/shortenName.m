function out = shortenName(n)
out = n;
if length(out) > 63
    out = extractBefore(out,64);
end
end