function [] = R1_plot_fig5(STFT, IF1, R1, IF2, R2, fname)

[Nfft, L] = size(STFT);

t = (0:L-1)/L;
M_indices = 1:round(L/30):L;

fname_noise = strcat(fname, '_noise');

figure;
imagesc(t, (0:Nfft-1)*L/Nfft, abs(STFT));
set(gca,'ydir','normal');
colormap(flipud(gray));
axis square;
xlabel('time', 'interpreter', 'latex');
ylabel('frequency', 'interpreter', 'latex');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 26);
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 26);
pbaspect([1 1 1]);
set(gcf, 'Position',  [0, 0, 1000, 1000])
savefig(fname_noise);
saveas(gcf, fname_noise, 'epsc');

set(groot, 'defaultLegendInterpreter', 'latex');

figure;
imagesc(t, (0:Nfft-1)*L/Nfft, abs(STFT));
set(gca,'ydir','normal');
colormap(flipud(gray));
axis square;
xlabel('time', 'interpreter', 'latex');
ylabel('frequency', 'interpreter', 'latex');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 26);
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 26);
hold on;
Y = [0.9290 0.6940 0.1250];
% plot(t, IF1 - R1, 'r-^',...
%     'Linewidth', 1, 'MarkerIndices', M_indices,...
%     'DisplayName', '$I_1^-$');
% plot(t, IF1 + R1, 'r-v',...
%     'Linewidth', 1, 'MarkerIndices', M_indices,...
%     'DisplayName', '$I_1^+$');
plot(t, IF1, 'k',...
    'Linewidth', 3,...
    'DisplayName', '$s_1^{fin}$');
% plot(t, IF2 - R2, 'r--^',...
%     'Linewidth', 1, 'MarkerIndices', M_indices,...
%     'DisplayName', '$I_2^-$');
% plot(t, IF2 + R2, 'r--v',...
%     'Linewidth', 1, 'MarkerIndices', M_indices,...
%     'DisplayName', '$I_2^+$');
plot(t, IF2, 'k--',...
    'Linewidth', 3,...
    'DisplayName', '$s_2^{fin}$');
hold off;
lgd = legend('Location', 'southeast');
lgd.FontSize = 24;
pbaspect([1 1 1]);
set(gcf, 'Position',  [0, 0, 1000, 1000])
savefig(fname);
saveas(gcf, fname, 'epsc');
% pause;
close all

end

