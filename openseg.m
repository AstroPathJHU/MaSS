function sim = openseg(simf)
for i1 = 1:3
    sim1 = imread(simf, i1 + 1);
    sim(:,i1) = reshape(sim1, [], 1);
    sim(:,i1) = double(sim(:,i1));
end

end