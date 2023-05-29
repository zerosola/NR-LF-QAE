function [alpha leftstd rightstd] = estimateaggdparam(vec)


gam   = 0.2:0.001:10;
r_gam = ((gamma(2./gam)).^2)./(gamma(1./gam).*gamma(3./gam));
throwAwayThresh = 0.0;
leftstd            = sqrt(mean((vec(vec<-throwAwayThresh)).^2));
rightstd           = sqrt(mean((vec(vec>throwAwayThresh)).^2));
gammahat           = leftstd/rightstd;

vec1=vec;
rhat               = (mean(abs(vec1)))^2/mean((vec1).^2);
rhatnorm           = (rhat*(gammahat^3 +1)*(gammahat+1))/((gammahat^2 +1)^2);
[min_difference, array_position] = min((r_gam - rhatnorm).^2);
alpha              = gam(array_position);

%% -------------------------------------------------------------------------------------
% structdis1 = vec;
% i = linspace(-5,5, 50);
% for j=1:49              %统计MSCN系数出现的频数
% n(j)=length(find(structdis1>=i(j) & structdis1<i(j+1)));
% x(j) = 0.5*(i(j)+i(j+1));
% end
% %     n=mapminmax(n,0,1);
% n = n./sum(n);
% 
% hold on
% plot(x,n,'-*b');
% set(gca,'XTick',[-5:1:5]);%x轴范围，间隔0.5
% legend('org');   %右上角标注
% xlabel('3D Log-Gabor Coefficients')  %x轴坐标描述
% ylabel('Probability') %y轴坐标描述
% grid on;  %加网格线
% set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1); %加网格虚线
% 
%             