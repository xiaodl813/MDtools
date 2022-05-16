function plot_kappa(t,rtc,filename)
% plot t-kappa
%   Input: t vector (double), rtc matrix (double) & filename (str)
%   Output: t-kappa plot (file)

font_size=15;
rtc_avg = mean(rtc, 2);

figure
plot(t, rtc, 'color', 0.5 * [1 1 1]);
hold on
plot(t, rtc_avg, 'linewidth', 3);
xlabel('Correlation Time (ps)', 'fontsize', font_size);
ylabel('$\kappa\ ({\rm Wm^{-1}K^{-1}})$', 'fontsize', font_size,...
                                          'interpreter','latex');
set(gca, 'fontsize', font_size);
saveas(gcf, filename);
close
end

