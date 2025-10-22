%% 修改参数以实现5km最大探测距离的RD谱
% 在你的主代码后面添加这段代码

%% 重新配置参数以实现5km最大探测距离
fprintf('\n========== Reconfiguring for 5km Max Range ==========\n');

% 要实现5km最大探测距离，需要：
% R_max = c * Tc / 2 = 5000 m
% 因此 Tc = 2 * R_max / c
R_max_target = 5000;  % 目标最大距离 5km
Tc_new = 2 * R_max_target / c;
fprintf('Required Tc for 5km range: %.3f μs\n', Tc_new*1e6);

% 但这会导致距离分辨率太差，我们需要权衡
% 建议方案：使用较大的Tc以获得合理的距离分辨率
% 例如：Tc = 3.33 ms 可以给出 1 MHz 有效带宽，150m 距离分辨率
Tc_optimal = 3.33e-3;  % 3.33 ms
B_eff_optimal = k * Tc_optimal;
range_res_optimal = c / (2 * B_eff_optimal);
R_max_optimal = c * Tc_optimal / 2;

fprintf('\nOptimal configuration:\n');
fprintf('Tc = %.3f ms\n', Tc_optimal*1e3);
fprintf('Effective bandwidth = %.2f MHz\n', B_eff_optimal/1e6);
fprintf('Range resolution = %.2f m\n', range_res_optimal);
fprintf('Max range = %.2f km\n', R_max_optimal/1e3);

% 使用优化参数重新处理
Tc = Tc_optimal;
N_chirp = floor(T/Tc);
N_per_chirp = round(Tc*Fs);

fprintf('Number of chirps: %d\n', N_chirp);
fprintf('Samples per chirp: %d\n', N_per_chirp);

%% 重新进行RD处理
N_total_samples = N_chirp * N_per_chirp;
s_if_rd = s_if_total(1:N_total_samples);
s_rd_matrix = reshape(s_if_rd, N_per_chirp, N_chirp);

% Range FFT
Nfft_range = 2^nextpow2(N_per_chirp);
win_range = hamming(N_per_chirp);
S_range = zeros(Nfft_range, N_chirp);
for i = 1:N_chirp
    S_range(:,i) = fft(s_rd_matrix(:,i) .* win_range, Nfft_range);
end

% 计算有效频率范围
f_beat_max = k * Tc;
freq_res_fft = Fs / Nfft_range;
N_range_bins = min(ceil(f_beat_max / freq_res_fft), Nfft_range/2);
S_range = S_range(1:N_range_bins, :);

% Doppler FFT
Nfft_doppler = 2^nextpow2(N_chirp);
win_doppler = hamming(N_chirp);
S_rd = zeros(N_range_bins, Nfft_doppler);
for i = 1:N_range_bins
    S_rd(i,:) = fftshift(fft(S_range(i,:) .* win_doppler.', Nfft_doppler));
end

%% 计算坐标轴
freq_res_range = Fs / Nfft_range;
freq_axis_range = (0:N_range_bins-1) * freq_res_range;
range_axis = freq_axis_range * c / (2*k);

doppler_res = 1 / (N_chirp * Tc);
doppler_axis = (-Nfft_doppler/2:Nfft_doppler/2-1) * doppler_res;

% 速度轴（从多普勒频率转换）
velocity_axis = doppler_axis * lambda / 2;  % v = fd * λ/2

% 计算理论位置
R_body_theory = R0 + v*(T/2);
fd_body_theory = 2*v/lambda;
v_body_theory = v;  % 真实速度

% 速度模糊检查
PRF = 1/Tc;
v_max_unamb = (PRF/2) * lambda/2;
fd_max_unamb = PRF/2;
fd_body_apparent = mod(fd_body_theory + fd_max_unamb, PRF) - fd_max_unamb;
v_body_apparent = fd_body_apparent * lambda / 2;

fprintf('\n--- Velocity Parameters ---\n');
fprintf('Max unambiguous velocity: %.2f m/s\n', v_max_unamb);
fprintf('Target velocity: %.2f m/s', v);
if abs(v) > v_max_unamb
    fprintf(' [AMBIGUOUS]\n');
    fprintf('Apparent velocity: %.2f m/s\n', v_body_apparent);
else
    fprintf('\n');
end

%% 绘制距离-多普勒频率图（限制到5km）
S_rd_dB = 20*log10(abs(S_rd) + 1e-12);
S_rd_dB = S_rd_dB - max(S_rd_dB(:));

% 显示范围设置
range_display = [0, 5000];  % 0-5km
doppler_display = [fd_body_apparent-3e3, fd_body_apparent+3e3];  % ±3kHz around target

figure('Name','RD Spectrum (5km Range, Doppler)','Color','w','Position',[200 100 1000 700]);
imagesc(doppler_axis/1e3, range_axis, S_rd_dB);
axis xy;
colormap('jet');
caxis([-60 0]);
colorbar;
xlabel('Doppler Frequency [kHz]', 'FontSize', 12);
ylabel('Range [m]', 'FontSize', 12);
title('Range-Doppler Spectrum (Max Range: 5 km)', 'FontSize', 14);
grid on;
xlim(doppler_display/1e3);
ylim(range_display);

hold on;
if abs(v) > v_max_unamb
    plot(fd_body_apparent/1e3, R_body_theory, 'ro', 'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', 'r');
    text(fd_body_apparent/1e3 + 0.3, R_body_theory + 200, ...
        sprintf('Target\n(%.0fm, %.1f Hz)', R_body_theory, fd_body_apparent), ...
        'Color', 'white', 'FontSize', 10, 'FontWeight', 'bold');
else
    plot(fd_body_theory/1e3, R_body_theory, 'ro', 'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', 'r');
    text(fd_body_theory/1e3 + 0.3, R_body_theory + 200, ...
        sprintf('Target\n(%.0fm, %.1f Hz)', R_body_theory, fd_body_theory), ...
        'Color', 'white', 'FontSize', 10, 'FontWeight', 'bold');
end

%% 绘制距离-速度图（推荐，更直观）
velocity_display = doppler_display * lambda / 2;  % 转换为速度范围

figure('Name','Range-Velocity Map (5km Range)','Color','w','Position',[220 120 1000 700]);
imagesc(velocity_axis, range_axis, S_rd_dB);
axis xy;
colormap('jet');
caxis([-60 0]);
h = colorbar;
ylabel(h, 'Normalized Power [dB]', 'FontSize', 11);
xlabel('Radial Velocity [m/s]', 'FontSize', 12);
ylabel('Range [m]', 'FontSize', 12);
title('Range-Velocity Map (Max Range: 5 km)', 'FontSize', 14);
grid on;
xlim(velocity_display);
ylim(range_display);

hold on;
if abs(v) > v_max_unamb
    plot(v_body_apparent, R_body_theory, 'ro', 'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', 'r');
    plot(v_body_theory, R_body_theory, 'bx', 'MarkerSize', 12, 'LineWidth', 2);
    text(v_body_apparent + 0.5, R_body_theory + 200, ...
        sprintf('Target (aliased)\n(%.0fm, %.1f m/s)', R_body_theory, v_body_apparent), ...
        'Color', 'white', 'FontSize', 10, 'FontWeight', 'bold');
    legend('RD Map','Target (apparent)','Target (true)', 'Location', 'northeast', ...
        'TextColor', 'white', 'Color', [0.2 0.2 0.2]);
else
    plot(v_body_theory, R_body_theory, 'ro', 'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', 'r');
    text(v_body_theory + 0.5, R_body_theory + 200, ...
        sprintf('Target\n(%.0fm, %.1f m/s)', R_body_theory, v_body_theory), ...
        'Color', 'white', 'FontSize', 10, 'FontWeight', 'bold');
    legend('RD Map','Target', 'Location', 'northeast', ...
        'TextColor', 'white', 'Color', [0.2 0.2 0.2]);
end

%% 打印最终参数
fprintf('\n========== Final Configuration ==========\n');
fprintf('Chirp duration (Tc):        %.3f ms\n', Tc*1e3);
fprintf('Effective bandwidth:        %.2f MHz\n', B_eff_optimal/1e6);
fprintf('Number of chirps:           %d\n', N_chirp);
fprintf('Range resolution:           %.2f m\n', range_res_optimal);
fprintf('Doppler resolution:         %.3f Hz\n', doppler_res);
fprintf('Velocity resolution:        %.3f m/s\n', doppler_res*lambda/2);
fprintf('Max range:                  %.2f km\n', R_max_optimal/1e3);
fprintf('Max unambiguous velocity:   %.2f m/s\n', v_max_unamb);
fprintf('Display range:              [%.0f, %.0f] m\n', range_display(1), range_display(2));
fprintf('Display velocity:           [%.1f, %.1f] m/s\n', velocity_display(1), velocity_display(2));
fprintf('==========================================\n\n');
