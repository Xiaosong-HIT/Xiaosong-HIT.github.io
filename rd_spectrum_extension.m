%% -------------------- RD谱生成（Range-Doppler Spectrum） --------------------
% 在上述代码基础上，继续生成 RD 谱
% 使用总的 IF 信号 s_if_total 进行 RD 处理

%% 参数设置
Tc = 1e-3;                      % 单个chirp的时长 [s]，可根据需要调整
N_chirp = floor(T/Tc);          % chirp个数
N_per_chirp = round(Tc*Fs);     % 每个chirp的采样点数

% 截取完整的chirp数据
N_total_samples = N_chirp * N_per_chirp;
s_if_rd = s_if_total(1:N_total_samples);

% 重排为矩阵：[距离维度(快时间) × 多普勒维度(慢时间)]
% 每一列对应一个chirp
s_rd_matrix = reshape(s_if_rd, N_per_chirp, N_chirp);

%% Range FFT（距离向FFT）
Nfft_range = 2^nextpow2(N_per_chirp);
% 对每一列（每个chirp）加窗并做FFT
win_range = hamming(N_per_chirp);
S_range = zeros(Nfft_range, N_chirp);
for i = 1:N_chirp
    S_range(:,i) = fft(s_rd_matrix(:,i) .* win_range, Nfft_range);
end

%% Doppler FFT（多普勒向FFT）
Nfft_doppler = 2^nextpow2(N_chirp);
% 对每一行（同一距离单元的所有chirp）加窗并做FFT
win_doppler = hamming(N_chirp);
S_rd = zeros(Nfft_range, Nfft_doppler);
for i = 1:Nfft_range
    S_rd(i,:) = fftshift(fft(S_range(i,:) .* win_doppler.', Nfft_doppler));
end

%% 计算坐标轴
% 距离轴：基于Range FFT后的频率bin转换为实际距离
% 在FMCW雷达中：f_beat = 2*k*R/c，因此 R = f_beat * c / (2*k)
freq_res_range = Fs / Nfft_range;                               % 频率分辨率 [Hz]
freq_axis_range = (0:Nfft_range-1) * freq_res_range;           % 频率轴 [Hz]
range_axis = freq_axis_range * c / (2*k);                       % 距离轴 [m]
range_res = c / (2*B);                                          % 理论距离分辨率 [m]

% 多普勒轴
doppler_res = 1 / (N_chirp * Tc);                               % 多普勒分辨率 [Hz]
doppler_axis = (-Nfft_doppler/2:Nfft_doppler/2-1) * doppler_res;  % 多普勒轴 [Hz]

%% 绘制 RD 谱（dB形式）
S_rd_dB = 20*log10(abs(S_rd) + 1e-12);
S_rd_dB = S_rd_dB - max(S_rd_dB(:));  % 归一化到0 dB

% 计算理论位置（用于确定显示范围）
R_body_theory = R0 + v*(T/2);          % 机体理论距离
fd_body_theory = 2*v/lambda;           % 机体理论多普勒

% 设置显示范围：以理论位置为中心
range_margin = 50;   % 距离方向显示余量 [m]
doppler_margin = 10e3;  % 多普勒方向显示余量 [Hz]
range_display = [max(0, R_body_theory-range_margin), R_body_theory+range_margin];
doppler_display = [fd_body_theory-doppler_margin, fd_body_theory+doppler_margin];

figure('Name','Range-Doppler Spectrum','Color','w','Position',[100 100 900 600]);
imagesc(doppler_axis/1e3, range_axis, S_rd_dB);
axis xy;
colormap('jet');
caxis([-60 0]);  % 动态范围 60 dB
colorbar;
xlabel('Doppler Frequency [kHz]');
ylabel('Range [m]');
title('Range-Doppler Spectrum (Total Signal: Body + Rotors)');
grid on;
xlim(doppler_display/1e3);
ylim(range_display);

% 添加参考线：理论机体位置
hold on;
plot(fd_body_theory/1e3, R_body_theory, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
legend('RD Spectrum','Body (theory)','Location','best');

%% 仅绘制旋翼的RD谱（用于对比）
s_if_rotor_rd = s_if_rotor(1:N_total_samples);
s_rd_matrix_rotor = reshape(s_if_rotor_rd, N_per_chirp, N_chirp);

S_range_rotor = zeros(Nfft_range, N_chirp);
for i = 1:N_chirp
    S_range_rotor(:,i) = fft(s_rd_matrix_rotor(:,i) .* win_range, Nfft_range);
end

S_rd_rotor = zeros(Nfft_range, Nfft_doppler);
for i = 1:Nfft_range
    S_rd_rotor(i,:) = fftshift(fft(S_range_rotor(i,:) .* win_doppler.', Nfft_doppler));
end

S_rd_rotor_dB = 20*log10(abs(S_rd_rotor) + 1e-12);
S_rd_rotor_dB = S_rd_rotor_dB - max(S_rd_rotor_dB(:));

figure('Name','Range-Doppler Spectrum (Rotors Only)','Color','w','Position',[120 120 900 600]);
imagesc(doppler_axis/1e3, range_axis, S_rd_rotor_dB);
axis xy;
colormap('jet');
caxis([-60 0]);
colorbar;
xlabel('Doppler Frequency [kHz]');
ylabel('Range [m]');
title('Range-Doppler Spectrum (Rotors Only, 4 rotors × 2 blades)');
grid on;
xlim(doppler_display/1e3);
ylim(range_display);

%% 仅绘制机体的RD谱（用于对比）
s_if_body_rd = s_if_body(1:N_total_samples);
s_rd_matrix_body = reshape(s_if_body_rd, N_per_chirp, N_chirp);

S_range_body = zeros(Nfft_range, N_chirp);
for i = 1:N_chirp
    S_range_body(:,i) = fft(s_rd_matrix_body(:,i) .* win_range, Nfft_range);
end

S_rd_body = zeros(Nfft_range, Nfft_doppler);
for i = 1:Nfft_range
    S_rd_body(i,:) = fftshift(fft(S_range_body(i,:) .* win_doppler.', Nfft_doppler));
end

S_rd_body_dB = 20*log10(abs(S_rd_body) + 1e-12);
S_rd_body_dB = S_rd_body_dB - max(S_rd_body_dB(:));

figure('Name','Range-Doppler Spectrum (Body Only)','Color','w','Position',[140 140 900 600]);
imagesc(doppler_axis/1e3, range_axis, S_rd_body_dB);
axis xy;
colormap('jet');
caxis([-60 0]);
colorbar;
xlabel('Doppler Frequency [kHz]');
ylabel('Range [m]');
title('Range-Doppler Spectrum (Body Only)');
grid on;
xlim(doppler_display/1e3);
ylim(range_display);
hold on;
plot(fd_body_theory/1e3, R_body_theory, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
legend('RD Spectrum','Body (theory)','Location','best');

%% 打印参数信息
fprintf('\n========== RD Spectrum Parameters ==========\n');
fprintf('Chirp duration (Tc):        %.3f ms\n', Tc*1e3);
fprintf('Number of chirps:           %d\n', N_chirp);
fprintf('Samples per chirp:          %d\n', N_per_chirp);
fprintf('Sampling rate:              %.2f MHz\n', Fs/1e6);
fprintf('Chirp bandwidth:            %.2f MHz\n', B/1e6);
fprintf('Chirp rate (k):             %.2e Hz/s\n', k);
fprintf('\n--- Resolution ---\n');
fprintf('Range resolution:           %.3f m\n', range_res);
fprintf('Doppler resolution:         %.3f Hz\n', doppler_res);
fprintf('Frequency resolution:       %.2f Hz\n', freq_res_range);
fprintf('\n--- Max Unambiguous ---\n');
fprintf('Max unambiguous range:      %.2f m\n', range_axis(end));
fprintf('Max unambiguous velocity:   %.2f m/s\n', doppler_axis(end)*lambda/2);
fprintf('\n--- Body Theoretical Values ---\n');
fprintf('Body position (R):          %.2f m\n', R_body_theory);
fprintf('Body velocity (v):          %.2f m/s\n', v);
fprintf('Body Doppler (fd):          %.2f Hz (%.3f kHz)\n', fd_body_theory, fd_body_theory/1e3);
fprintf('Expected beat freq:         %.2f Hz (%.3f kHz)\n', 2*k*R_body_theory/c, 2*k*R_body_theory/c/1e3);
fprintf('\n--- Display Range ---\n');
fprintf('Range display:              [%.1f, %.1f] m\n', range_display(1), range_display(2));
fprintf('Doppler display:            [%.1f, %.1f] kHz\n', doppler_display(1)/1e3, doppler_display(2)/1e3);
fprintf('==========================================\n\n');
