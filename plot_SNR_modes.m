function plot_SNR_modes(SNR_IN, SNR_NEW_LCR, SNR_NEW_MR, SNR_S_MR, SNR_MB_MR)
    c_red = [0.6350 0.0780 0.1840];
    c_blue = [0 0.4470 0.7410];
    c_green = [60/255 179/255 113/255];

    set(groot, 'defaultLegendInterpreter', 'latex');
    figure;
    hold on;
    plot(SNR_IN, SNR_NEW_MR(1, :), 'k-^',...
        'DisplayName', 'RRP-MR $f_1$',...
        'MarkerSize', 10, 'LineWidth', 2);
    
    plot(SNR_IN, SNR_MB_MR(1, :), '-o', 'Color', c_blue,...
        'DisplayName', 'MB-MR $f_1$',...
        'MarkerSize', 10, 'LineWidth', 2);
    
    plot(SNR_IN, SNR_S_MR(1, :), '-s', 'Color', c_red,...
        'DisplayName', 'S-MR $f_1$',...
        'MarkerSize', 10, 'LineWidth', 2);
    
    plot(SNR_IN, SNR_NEW_MR(2, :), 'k--^',...
        'DisplayName', 'RRP-MR $f_2$',...
        'MarkerSize', 10, 'LineWidth', 2);
    
    plot(SNR_IN, SNR_MB_MR(2, :), '--o', 'Color', c_blue,...
        'DisplayName', 'MB-MR $f_2$',...
        'MarkerSize', 10, 'LineWidth', 2);
    
    plot(SNR_IN, SNR_S_MR(2, :), '--s', 'Color', c_red,...
        'DisplayName', 'S-MR $f_2$',...
        'MarkerSize', 10, 'LineWidth', 2);
    
    plot(SNR_IN, SNR_NEW_LCR(1, :), '-', 'Color', c_green,...
        'DisplayName', 'RRP-MR-LCR $f_1$',...
        'MarkerSize', 10, 'LineWidth', 2);
    
    plot(SNR_IN, SNR_NEW_LCR(2, :), '--', 'Color', c_green,...
        'DisplayName', 'RRP-MR-LCR $f_2$',...
        'MarkerSize', 10, 'LineWidth', 2);
    
    hold off;
    xlabel('input SNR', 'interpreter', 'latex');
    ylabel('output SNR', 'interpreter', 'latex');
    xlim([SNR_IN(1), SNR_IN(length(SNR_IN))]);

    lgd = legend('Location', 'northwest');
    lgd.FontSize = 24;
    xAX = get(gca,'XAxis');
    set(xAX,'FontSize', 26);
    yAX = get(gca,'YAxis');
    set(yAX,'FontSize', 26);
    pbaspect([1 1 1]);
    set(gcf, 'Position',  [0, 0, 1000, 1000])
end