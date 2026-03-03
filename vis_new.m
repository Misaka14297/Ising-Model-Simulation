% 读取数据
data = readtable('ising_data_combined.csv');
L_values = unique(data.L);
algorithms = unique(data.Algorithm);
colors = [0.1, 0.4, 0.8; 0.9, 0.6, 0.2; 0.5, 0.2, 0.7]; % 蓝、橙、紫
markers = {'-', '--'}; % Metropolis: 实线, Wolff: 虚线
line_width = 1.5;

% 子图标题和标签
titles = {'平均能量 per site vs 温度', '平均绝对磁化强度 vs 温度', ...
          '比热 vs 温度', '磁化率 vs 温度'};
y_labels = {'E / L^2', '|m|', 'C_v', 'χ'};
quantities = {'Avg_energy', 'Avg_abs_magnetization', ...
              'Energy_fluctuation', 'Magnetization_fluctuation'};
std_fields = {'Std_energy', 'Std_abs_magnetization', ...
              'Std_energy_fluctuation', 'Std_magnetization_fluctuation'};

% 循环绘制四个物理量
for p_idx = 1:4
    % 新建图形并清除状态
    figure('Position', [100, 100, 800, 600]);
    clf;
    hold on;
    % 存储 plot 句柄用于图例
    plot_handles = [];
    plot_labels = {};
    for alg_idx = 1:length(algorithms)
        algorithm = algorithms{alg_idx};
        for l_idx = 1:length(L_values)
            L = L_values(l_idx);
            % 筛选数据
            idx = strcmp(data.Algorithm, algorithm) & (data.L == L);
            if ~any(idx)
                warning('No data for Algorithm=%s, L=%d, Quantity=%s', ...
                        algorithm, L, quantities{p_idx});
                continue;
            end
            T = data.Temperature(idx);
            quantity = data.(quantities{p_idx})(idx);
            std_quantity = data.(std_fields{p_idx})(idx);
            % 归一化处理
            if strcmp(quantities{p_idx}, 'Avg_energy')
                quantity = quantity / (L * L);
                std_quantity = std_quantity / (L * L);
            end
            % 绘制曲线
            h = plot(T, quantity, 'Color', colors(l_idx,:), ...
                     'LineStyle', markers{alg_idx}, 'LineWidth', line_width, ...
                     'DisplayName', sprintf('%s (L=%d)', algorithm, L));
            plot_handles = [plot_handles, h];
            plot_labels = [plot_labels, sprintf('%s (L=%d)', algorithm, L)];
            % 绘制误差阴影
            fill([T; flip(T)], [quantity + std_quantity; flip(quantity - std_quantity)], ...
                 colors(l_idx,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
                 'HandleVisibility', 'off');
        end
    end
    hold off;
    % 设置图形属性
    title(titles{p_idx}, 'FontSize', 14);
    xlabel('温度 T', 'FontSize', 12);
    ylabel(y_labels{p_idx}, 'FontSize', 12);
    legend(plot_handles, plot_labels, 'Location', 'best', 'FontSize', 10);
    grid on;
    set(gca, 'FontSize', 10);
    % 保存为 PDF
    print(gcf, sprintf('ising_%s.pdf', quantities{p_idx}), '-dpdf', '-r600');
    close(gcf);
end