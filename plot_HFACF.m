function plot_HFACF(t,hac,filename)
% plot t-HFACF
%   Input: t vector (double), hac matrix (double) & filename (str)
%   Output: t-HFACF plot (file)

font_size=15;
hac_avg = mean(hac, 2);

figure
plot(t, hac ./ hac(1, :), 'color', 0.5 * [1 1 1]);
hold on
plot(t, hac_avg ./ hac_avg(1), 'linewidth',3);
xlabel('Correlation Time (ps)', 'fontsize', font_size);
ylabel('HCACF (Normalized)', 'fontsize', font_size);
set(gca, 'fontsize', font_size);
saveas(gcf, filename);
close
end



