% WEIGHT: NUMBER OF GRAINS
data = GRAINS.aspectRatio;             % generate normal random 
h = histogram(data,'Normalization','probability','FaceColor','w')   % Probability normalization
title('RD-TD')
xlabel('Grain shape aspect ratio,-') 
ylabel('Number fraction,-')

Nbins = h.NumBins;
edges = h.BinEdges; 
x = zeros(1,Nbins);
for counter=1:Nbins
    midPointShift = abs(edges(counter)-edges(counter+1))/2;
    x(counter) = edges(counter)+midPointShift;
end

mu = mean(data);
sigma = std(data);

f = exp(-(x-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));

%Normalize from pdf to probability
f=f/sum(f);

hold on;
plot(x,f,'LineWidth',1.5,'color','b')