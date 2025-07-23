% Color scheme from my web page, adapted from arXiv
x = linspace(-9,9,181);
g = exp(-x.^2/10);
c = cos(3*x);
figure(1)
plot(x,g,'LineWidth',7,'Color','#000000');
hold on
plot(x,g.*c*.986,'LineWidth',5,'Color','#B31B1B')
hold off
axis([-9,9,-1.1, 1.1])
set(gca,'visible','off')

