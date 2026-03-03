% 读取数据
data = readtable('ising_data_combined.csv');
L_values = unique(data.L);
T_values = unique(data.Temperature);
algorithms = unique(data.Algorithm);

% 初始化绘图
figure('Position', [100, 100, 1200, 900]);
colors = [0.1, 0.4, 0.8; 0.9, 0.6, 0.2; 0.5, 0.2, 0.7]; % 蓝、橙、紫
markers = {'-', '--'}; % Metropolis: 实线, Wolff: 虚线
line_width = 1.5;

% 子图标题和标签
titles = {'平均能量 per site vs 温度', '平均绝对磁化强度 vs 温度', ...
          '能量波动 vs 温度', '磁化强度波动 vs 温度'};
y_labels = {'E / L^2', '|m|', 'C_v', 'χ'};
quantities = {'Avg_energy', 'Avg_abs_magnetization', ...
              'Energy_fluctuation', 'Magnetization_fluctuation'};
std_fields = {'Std_energy', 'Std_abs_magnetization', ...
              'Std_energy_fluctuation', 'Std_magnetization_fluctuation'};

% 循环绘制四个物理量
for p_idx = 1:4
    subplot(2, 2, p_idx);
    hold on;
    for alg_idx = 1:length(algorithms)
        algorithm = algorithms{alg_idx};
        for l_idx = 1:length(L_values)
            L = L_values(l_idx);
            % 筛选数据
            idx = strcmp(data.Algorithm, algorithm) & (data.L == L);
            T = data.Temperature(idx);
            quantity = data.(quantities{p_idx})(idx);
            std_quantity = data.(std_fields{p_idx})(idx);
            % 归一化处理
            if strcmp(quantities{p_idx}, 'Avg_energy')
                quantity = quantity / (L * L);
                std_quantity = std_quantity / (L * L);
            end
            % 绘制曲线和误差阴影
            plot(T, quantity, 'Color', colors(l_idx,:), ...
                 'LineStyle', markers{alg_idx}, 'LineWidth', line_width, ...
                 'DisplayName', sprintf('%s (L=%d)', algorithm, L));
            fill([T; flip(T)], [quantity + std_quantity; flip(quantity - std_quantity)], ...
                 colors(l_idx,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
                 'HandleVisibility', 'off');
        end
    end
    hold off;
    title(titles{p_idx}, 'FontSize', 12);
    xlabel('温度 T', 'FontSize', 10);
    ylabel(y_labels{p_idx}, 'FontSize', 10);
    legend('Location', 'best', 'FontSize', 8);
    grid on;
end
sgtitle('2D Ising Model - Metropolis vs Wolff', 'FontSize', 14);
saveas(gcf, 'ising_results.png');


% 生成表格
fprintf('Table: Peak Values of Specific Heat and Susceptibility\n');
fprintf('-------------------------------------------------\n');
fprintf('Algorithm | L  | T_peak(C_v) | C_v_peak | Std_C_v | T_peak(χ) | χ_peak | Std_χ\n');
for alg_idx = 1:length(algorithms)
    algorithm = algorithms{alg_idx};
    for l_idx = 1:length(L_values)
        L = L_values(l_idx);
        idx = strcmp(data.Algorithm, algorithm) & (data.L == L);
        T = data.Temperature(idx);
        Cv = data.Energy_fluctuation(idx);
        Chi = data.Magnetization_fluctuation(idx);
        Std_Cv = data.Std_energy_fluctuation(idx);
        Std_Chi = data.Std_magnetization_fluctuation(idx);
        [Cv_max, Cv_idx] = max(Cv);
        [Chi_max, Chi_idx] = max(Chi);
        fprintf('%s | %d | %.2f | %.4f | %.4f | %.2f | %.4f | %.4f\n', ...
                algorithm, L, T(Cv_idx), Cv_max, Std_Cv(Cv_idx), ...
                T(Chi_idx), Chi_max, Std_Chi(Chi_idx));
    end
end