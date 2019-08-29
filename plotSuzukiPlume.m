
H = 15;

k = 4:9;



h       = linspace(0,H,50);
leg     = cell(numel(k),1);

figure, hold on, box on

for ik = 1:numel(k)
    tmp = zeros(numel(h),1);
    
    for ih = 1:numel(h)
        tmp(ih) = (1-h(ih)/H) * exp(-k(ik)*(1-h(ih)/H));
    end
    plot(tmp, h)
    leg{ik} = ['k = ', num2str(k(ik))];
end

legend(leg)

% ds = ( (k.^2 .* (1-z./H) .* exp(k(z./H-1))) ) ./ ( H(1-(1+k) .* exp(-k)) );
