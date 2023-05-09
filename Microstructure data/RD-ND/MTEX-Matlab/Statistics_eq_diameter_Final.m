% WEIGHT: NUMBER OF GRAINS
figure(1)
data = GRAINS.equivalentDiameter;             % generate normal random 
h = histogram(data,'Normalization','probability')   % Probability normalization (no PDF)
title('RD-TD')
xlabel('Grain equivalent diameter, µm') 
ylabel('Number fraction,-')


% WEIGHT: GRAIN AREA
vv = GRAINS.equivalentDiameter; % values
ww = GRAINS.area; % weights
nbins = 100;
figure(2)
[histw, intervals] = histwc(vv, ww, nbins);
%Example Visualise (Normalized to probability):
bar(intervals, histw/sum(histw),'FaceColor','w')
title('RD-TD')
xlabel('Grain equivalent diameter, µm') 
ylabel('Area fraction,-')


Nbins = h.NumBins;
edges = h.BinEdges; 
x = zeros(1,Nbins);
for counter=1:Nbins
    midPointShift = abs(edges(counter)-edges(counter+1))/2;
    x(counter) = edges(counter)+midPointShift;
end

pd = fitdist(GRAINS.equivalentDiameter,'lognormal','frequency',GRAINS.area_int)

mu = pd.mu;
sigma = pd.sigma;
average=exp(mu+sigma^2/2)

f = exp(-(log(x)-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi)*x);

%Normalize from pdf to probability
f=f/sum(f);

hold on;
plot(x,f,'LineWidth',1.5,'color','b')