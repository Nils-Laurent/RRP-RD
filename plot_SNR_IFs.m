function plot_SNR_IFs(SNR_IN, SNR_IF_NEW, SNR_IF_C_RD, SNR_IF_MB_RD)
    c_red = [0.6350 0.0780 0.1840];
    c_blue = [0 0.4470 0.7410];
    
    set(groot, 'defaultLegendInterpreter', 'latex');
    figure;
    hold on;
    plot(SNR_IN, SNR_IF_NEW(1, :), 'k-^',...
        'DisplayName', 'RRP-RD $f_1$',...
        'MarkerSize', 10, 'LineWidth', 2);
    
    plot(SNR_IN, SNR_IF_MB_RD(1, :), '-o', 'Color', c_blue,...
        'DisplayName', 'MB-RD $f_1$',...
        'MarkerSize', 10, 'LineWidth', 2);
    
    plot(SNR_IN, SNR_IF_C_RD(1, :), '-s', 'Color', c_red,...
        'DisplayName', 'S-RD $f_1$',...
        'MarkerSize', 10, 'LineWidth', 2);
    
    plot(SNR_IN, SNR_IF_NEW(2, :), 'k--^',...
        'DisplayName', 'RRP-RD $f_2$',...
        'MarkerSize', 10, 'LineWidth', 2);
    
    plot(SNR_IN, SNR_IF_MB_RD(2, :), '--o', 'Color', c_blue,...
        'DisplayName', 'MB-RD $f_2$',...
        'MarkerSize', 10, 'LineWidth', 2);
    
    plot(SNR_IN, SNR_IF_C_RD(2, :), '--s', 'Color', c_red,...
        'DisplayName', 'S-RD $f_2$',...
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

