function md = calcemd(data, mu, SIGMA)

nobservations = size(data, 1);
md = zeros(nobservations, 1);
for m=1 : nobservations
    md(m) = (data(m,:) - mu)*inv(SIGMA)*(data(m,:) - mu)';
end

end