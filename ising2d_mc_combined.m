function ising2d_mc_combined(L_values, Tmin, Tmax, Tstep, steps_per_T_base, averaging_runs, equilibration_steps, algorithms)
    % 二维Ising模型的蒙特卡洛模拟（Metropolis和Wolff算法），支持多个系统尺寸和算法选择
    % 参数：
    % L_values - 系统尺寸数组，例如 [16, 32, 64]
    % Tmin - 最低温度
    % Tmax - 最高温度
    % Tstep - 温度步长
    % steps_per_T_base - 基础蒙特卡洛步数（低温区会自适应减少）
    % averaging_runs - 平均模拟次数
    % equilibration_steps - 热平衡步骤数
    % algorithms - 字符串或字符串数组，指定要运行的算法，例如 {'metropolis', 'wolff'}

    T_values = linspace(Tmax, Tmin, round((Tmax - Tmin) / Tstep) + 1);
    n_T = length(T_values);
    n_L = length(L_values);
    n_algorithms = length(algorithms);

    % 初始化存储所有结果的结构体
    results = struct();

    tic
    for alg_idx = 1:n_algorithms
        algorithm = lower(algorithms{alg_idx}); % 转换为小写
        results.(algorithm) = struct();
        fprintf('开始模拟算法: %s\n', algorithm);

        for l_idx = 1:n_L
            L = L_values(l_idx);
            results.(algorithm).(['L', num2str(L)]) = struct(...
                'Avg_energy', zeros(1, n_T), ...
                'Avg_abs_magnetization', zeros(1, n_T), ...
                'Energy_fluctuation', zeros(1, n_T), ...
                'Magnetization_fluctuation', zeros(1, n_T), ...
                'Std_energy', zeros(1, n_T), ...
                'Std_abs_magnetization', zeros(1, n_T), ...
                'Std_energy_fluctuation', zeros(1, n_T), ...
                'Std_magnetization_fluctuation', zeros(1, n_T));

            fprintf('  系统尺寸 L = %d\n', L);

            for run = 1:averaging_runs
                fprintf('   运行第 %d / %d 次\n', run, averaging_runs);
                spins = randi([0, 1], L, L) * 2 - 1;
                total_energy = calculate_total_energy_combined(spins, L); % 使用统一的能量计算函数

                for t = 1:n_T
                    T = T_values(t);
                    fprintf('    模拟温度 T = %.2f\n', T);

                    steps_per_T = steps_per_T_base;
                    if T < 2.0
                        steps_per_T = round(steps_per_T_base * 2);
                    end

                    % 热平衡
                    for step = 1:equilibration_steps
                        if strcmp(algorithm, 'metropolis')
                            [spins, delta_energy] = metropolis_update_combined(spins, L, T); % 使用统一的更新函数
                            total_energy = total_energy + delta_energy;
                        elseif strcmp(algorithm, 'wolff')
                            [spins, delta_energy] = wolff_cluster_update_combined(spins, L, T); % 使用统一的更新函数
                            total_energy = total_energy + delta_energy;
                        end
                    end

                    energy_sum = 0;
                    energy_sq_sum = 0;
                    abs_mag_sum = 0;
                    mag_sq_sum = 0;

                    for step = 1:steps_per_T
                        if strcmp(algorithm, 'metropolis')
                            [spins, delta_energy] = metropolis_update_combined(spins, L, T);
                            total_energy = total_energy + delta_energy;
                        elseif strcmp(algorithm, 'wolff')
                            [spins, delta_energy] = wolff_cluster_update_combined(spins, L, T);
                            total_energy = total_energy + delta_energy;
                        end
                        total_magnetization = sum(spins(:)) / (L * L);
                        abs_magnetization = abs(total_magnetization);

                        energy_sum = energy_sum + total_energy;
                        energy_sq_sum = energy_sq_sum + total_energy^2;
                        abs_mag_sum = abs_mag_sum + abs_magnetization;
                        mag_sq_sum = mag_sq_sum + total_magnetization^2;
                    end

                    avg_energy_run = energy_sum / steps_per_T;
                    avg_abs_mag_run = abs_mag_sum / steps_per_T;
                    avg_mag_sq_run = mag_sq_sum / steps_per_T;

                    energy_fluctuation_run = (energy_sq_sum / steps_per_T - avg_energy_run^2) / (T^2 * L^2);
                    mag_fluctuation_run = (avg_mag_sq_run - avg_abs_mag_run^2) / T;

                    results.(algorithm).(['L', num2str(L)]).Avg_energy(t) = results.(algorithm).(['L', num2str(L)]).Avg_energy(t) + avg_energy_run;
                    results.(algorithm).(['L', num2str(L)]).Avg_abs_magnetization(t) = results.(algorithm).(['L', num2str(L)]).Avg_abs_magnetization(t) + avg_abs_mag_run;
                    results.(algorithm).(['L', num2str(L)]).Energy_fluctuation(t) = results.(algorithm).(['L', num2str(L)]).Energy_fluctuation(t) + energy_fluctuation_run;
                    results.(algorithm).(['L', num2str(L)]).Magnetization_fluctuation(t) = results.(algorithm).(['L', num2str(L)]).Magnetization_fluctuation(t) + mag_fluctuation_run;

                    results.(algorithm).(['L', num2str(L)]).Std_energy(t) = results.(algorithm).(['L', num2str(L)]).Std_energy(t) + avg_energy_run^2;
                    results.(algorithm).(['L', num2str(L)]).Std_abs_magnetization(t) = results.(algorithm).(['L', num2str(L)]).Std_abs_magnetization(t) + avg_abs_mag_run^2;
                    results.(algorithm).(['L', num2str(L)]).Std_energy_fluctuation(t) = results.(algorithm).(['L', num2str(L)]).Std_energy_fluctuation(t) + energy_fluctuation_run^2;
                    results.(algorithm).(['L', num2str(L)]).Std_magnetization_fluctuation(t) = results.(algorithm).(['L', num2str(L)]).Std_magnetization_fluctuation(t) + mag_fluctuation_run^2;
                end
            end

            % 计算平均值和标准差
            results.(algorithm).(['L', num2str(L)]).Avg_energy = results.(algorithm).(['L', num2str(L)]).Avg_energy / averaging_runs;
            results.(algorithm).(['L', num2str(L)]).Avg_abs_magnetization = results.(algorithm).(['L', num2str(L)]).Avg_abs_magnetization / averaging_runs;
            results.(algorithm).(['L', num2str(L)]).Energy_fluctuation = results.(algorithm).(['L', num2str(L)]).Energy_fluctuation / averaging_runs;
            results.(algorithm).(['L', num2str(L)]).Magnetization_fluctuation = results.(algorithm).(['L', num2str(L)]).Magnetization_fluctuation / averaging_runs;

            results.(algorithm).(['L', num2str(L)]).Std_energy = sqrt(results.(algorithm).(['L', num2str(L)]).Std_energy / averaging_runs - (results.(algorithm).(['L', num2str(L)]).Avg_energy).^2) / sqrt(averaging_runs);
            results.(algorithm).(['L', num2str(L)]).Std_abs_magnetization = sqrt(results.(algorithm).(['L', num2str(L)]).Std_abs_magnetization / averaging_runs - (results.(algorithm).(['L', num2str(L)]).Avg_abs_magnetization).^2) / sqrt(averaging_runs);
            results.(algorithm).(['L', num2str(L)]).Std_energy_fluctuation = sqrt(results.(algorithm).(['L', num2str(L)]).Std_energy_fluctuation / averaging_runs - (results.(algorithm).(['L', num2str(L)]).Energy_fluctuation).^2) / sqrt(averaging_runs);
            results.(algorithm).(['L', num2str(L)]).Std_magnetization_fluctuation = sqrt(results.(algorithm).(['L', num2str(L)]).Std_magnetization_fluctuation / averaging_runs - (results.(algorithm).(['L', num2str(L)]).Magnetization_fluctuation).^2) / sqrt(averaging_runs);
        end
    end
    combined_time = toc;
    fprintf('组合算法运行时间: %.2f 秒\n', combined_time);
    disp('模拟结束后 results 结构体的内容：');
    disp(results);

    % 保存数据
    save_combined_data(L_values, T_values, results);

    % 可视化结果
    plot_combined_results(L_values, T_values, results);
end

function [spins, delta_energy_total] = metropolis_update_combined(spins, L, T)
    delta_energy_total = 0;
    for attempt = 1:L*L
        i = randi(L);
        j = randi(L);
        s = spins(i, j);
        neighbors = get_neighbors_combined(i, j, L);
        neighbor_sum = 0;
        for n = 1:4
            ni = neighbors(n, 1);
            nj = neighbors(n, 2);
            neighbor_sum = neighbor_sum + spins(ni, nj);
        end
        delta_E = 2 * s * neighbor_sum;
        if delta_E <= 0 || rand() < exp(-delta_E / T)
            spins(i, j) = -s;
            delta_energy_total = delta_energy_total + delta_E;
        end
    end
end

function [spins, delta_energy] = wolff_cluster_update_combined(spins, L, T)
    seed_i = randi(L);
    seed_j = randi(L);

    cluster = false(L, L);
    cluster(seed_i, seed_j) = true;
    to_check = [seed_i, seed_j];

    seed_spin = spins(seed_i, seed_j);
    P_add = 1 - exp(-2 / T);

    while ~isempty(to_check)
        i = to_check(1, 1);
        j = to_check(1, 2);
        to_check(1, :) = [];

        neighbors = get_neighbors_combined(i, j, L);

        for n = 1:size(neighbors, 1)
            ni = neighbors(n, 1);
            nj = neighbors(n, 2);
            if ~cluster(ni, nj) && spins(ni, nj) == seed_spin && rand() < P_add
                cluster(ni, nj) = true;
                to_check = [to_check; ni, nj];
            end
        end
    end
    delta_energy = 0;
    for i = 1:L
        for j = 1:L
            if cluster(i, j)
                neighbors = get_neighbors_combined(i, j, L);
                for n = 1:size(neighbors, 1)
                    ni = neighbors(n, 1);
                    nj = neighbors(n, 2);
                    if ~cluster(ni, nj)
                        delta_energy = delta_energy + 2 * spins(i, j) * spins(ni, nj);
                    end
                end
            end
        end
    end
    spins(cluster) = -spins(cluster);
end

function neighbors = get_neighbors_combined(i, j, L)
    i_prev = mod(i - 2, L) + 1;
    i_next = mod(i, L) + 1;
    j_prev = mod(j - 2, L) + 1;
    j_next = mod(j, L) + 1;
    neighbors = [i_prev, j; i_next, j; i, j_prev; i, j_next];
end

function E = calculate_total_energy_combined(spins, L)
    E = 0;
    for i = 1:L
        for j = 1:L
            s = spins(i, j);
            neighbors = get_neighbors_combined(i, j, L);
            for n = 1:4
                ni = neighbors(n, 1);
                nj = neighbors(n, 2);
                E = E - s * spins(ni, nj);
            end
        end
    end
    E = E / 2;
end

function save_combined_data(L_values, T_values, results)
    % 保存为 .mat 格式
    mat_filename = 'ising_data_combined.mat';
    save(mat_filename, 'L_values', 'T_values', 'results');
    fprintf('组合算法数据已保存到 %s\n', mat_filename);

    % 保存为 .csv 格式
    csv_filename = 'ising_data_combined.csv';
    file_id = fopen(csv_filename, 'w');
    if file_id == -1
        fprintf('无法打开文件 %s 进行写入。\n', csv_filename);
        return;
    end

    % 写入 CSV 文件头
    fprintf(file_id, 'Algorithm,L,Temperature,Avg_energy,Std_energy,Avg_abs_magnetization,Std_abs_magnetization,Energy_fluctuation,Std_energy_fluctuation,Magnetization_fluctuation,Std_magnetization_fluctuation\n');

    algorithm_names = fieldnames(results);
    for i = 1:length(algorithm_names)
        algorithm = algorithm_names{i};
        for j = 1:length(L_values)
            L = L_values(j);
            L_str = ['L', num2str(L)];
            if isfield(results.(algorithm), L_str)
                for t_idx = 1:length(T_values)
                    T = T_values(t_idx);
                    data = results.(algorithm).(L_str);
                    fprintf(file_id, '%s,%d,%.4f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f\n', ...
                        algorithm, L, T, ...
                        data.Avg_energy(t_idx), data.Std_energy(t_idx), ...
                        data.Avg_abs_magnetization(t_idx), data.Std_abs_magnetization(t_idx), ...
                        data.Energy_fluctuation(t_idx), data.Std_energy_fluctuation(t_idx), ...
                        data.Magnetization_fluctuation(t_idx), data.Std_magnetization_fluctuation(t_idx));
                end
            end
        end
    end

    fclose(file_id);
    fprintf('组合算法数据已保存到 %s\n', csv_filename);
end

function plot_combined_results(L_values, T_values, results)
    figure('Position', [100, 100, 1200, 900]);
    colors = [ % 更友好的颜色
        0.1, 0.4, 0.8; % 蓝色
        0.9, 0.6, 0.2; % 橙色
        0.5, 0.2, 0.7; % 紫色
        0.2, 0.6, 0.3; % 绿色 (稍暗)
        0.8, 0.3, 0.3; % 红色 (稍淡)
        0.4, 0.8, 0.8  % 青色 (稍亮)
    ];
    markers = {'-', '--'}; % 线型用于区分算法
    line_width = 1.5;
    physical_quantities = {'Avg_energy', 'Avg_abs_magnetization', 'Energy_fluctuation', 'Magnetization_fluctuation'};
    y_labels = {'平均能量 per site', '平均绝对磁化强度', '能量波动', '磁化强度波动'};
    titles = {'平均能量 per site vs 温度', '平均绝对磁化强度 vs 温度', '能量波动 vs 温度', '磁化强度波动 vs 温度'};
    algorithm_names = fieldnames(results);
    for alg_idx = 1:length(algorithm_names)
        algorithm = algorithm_names{alg_idx};
        for l_idx = 1:length(L_values)
            L = L_values(l_idx);
            if isfield(results.(algorithm), ['L', num2str(L)])
                for p_idx = 1:length(physical_quantities) % 将 p_idx 的循环放在这里
                    subplot(2, 2, p_idx);
                    hold on;
                    data = results.(algorithm).(['L', num2str(L)]).(physical_quantities{p_idx});
                    std_field_name = strrep(['Std_', lower(physical_quantities{p_idx})], 'avg_', ''); % 移除 "avg_" 并转换为小写
                    std_data = results.(algorithm).(['L', num2str(L)]).(std_field_name);
                    if strcmp(physical_quantities{p_idx}, 'Avg_energy')
                        data = data / (L * L);
                        std_data = std_data / (L * L);
                    end
                    plot(T_values, data, 'DisplayName', sprintf('%s (L=%d)', algorithm, L), 'Color', colors(l_idx,:), 'LineStyle', markers{alg_idx}, 'LineWidth', line_width);
                    errorbar(T_values, data, std_data, 'Color', colors(l_idx,:), 'LineStyle', markers{alg_idx}, 'LineStyle', 'none', 'HandleVisibility', 'off'); % 添加 'HandleVisibility', 'off'
                    hold off;
                    title(titles{p_idx});
                    xlabel('温度 T');
                    ylabel(y_labels{p_idx});
                    legend('Location', 'best'); % 确保每个子图都有图例
                    grid on;
                end
            else
                warning('在 results.%s 中找不到 L=%d 的数据。', algorithm, L);
            end
        end
    end
    sgtitle('2D Ising Model - Metropolis vs Wolff');
end