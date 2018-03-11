T = 128;
n = 1000;
name = 'UncorNormal';
t = linspace(-3,3,T);
[pdfxy,x] = MakeDensity2d(name,T,n);
figure
subplot(1,2,1)
contour3(t,t,pdfxy,20);hold on
plot(x(:,1),x(:,2),'o','MarkerSize',5)
caxis([min(pdfxy(:))-.5*range(pdfxy(:)),max(pdfxy(:))]);
title(name)
view(2)
subplot(1,2,2)
contour3(t,t,pdfxy,20);hold on
plot(x(:,1),x(:,2),'o','MarkerSize',5)
caxis([min(pdfxy(:))-.5*range(pdfxy(:)),max(pdfxy(:))]);
title(name)