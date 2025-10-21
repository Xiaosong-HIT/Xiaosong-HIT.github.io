# FMCW雷达参数配置指南

## 问题诊断

你当前的参数设置：
- `Tc = 27e-3` (27 ms)
- `T = 500e-3` (500 ms)
- 结果：`v_max = 0.14 m/s` ❌ **太小了！**

## 根本原因

**关键公式：**
```
v_max = λ / (4 × Tc)
```

其中：
- λ = 0.015 m (波长)
- Tc = chirp时长

**重要规律：Tc 与 v_max 成反比！**

| Tc (ms) | PRF (Hz) | v_max (m/s) | 说明 |
|---------|----------|-------------|------|
| 1.0     | 1000     | 3.75        | 你最初的设置 |
| 27.0    | 37       | 0.14        | 你当前的设置（错误！） |
| 0.15    | 6667     | 25.0        | 推荐设置 ✓ |
| 0.1     | 10000    | 37.5        | 更高速度测量 |

---

## 正确的参数配置

### 方案 1：推荐配置（平衡性能）

```matlab
T = 500e-3;      % 总时长 500 ms
Tc = 0.15e-3;    % Chirp时长 0.15 ms (150 μs)
```

**性能指标：**
- ✅ PRF = 6667 Hz
- ✅ v_max = ±25 m/s (足够测量 -20 m/s)
- ✅ R_max_chirp = 22.5 km
- ✅ N_chirp = 3333 (充足的多普勒采样)
- ✅ 多普勒分辨率 ≈ 2 Hz

### 方案 2：高速度测量

```matlab
T = 500e-3;      % 总时长 500 ms
Tc = 0.1e-3;     % Chirp时长 0.1 ms (100 μs)
```

**性能指标：**
- ✅ PRF = 10000 Hz
- ✅ v_max = ±37.5 m/s (更高速度范围)
- ✅ R_max_chirp = 15 km
- ✅ N_chirp = 5000

### 方案 3：长距离测量（如果需要）

如果你的应用需要测量更远距离（>100 km），可以适当增大 Tc，但会牺牲速度测量范围：

```matlab
T = 500e-3;      % 总时长 500 ms
Tc = 1.0e-3;     % Chirp时长 1 ms
```

**性能指标：**
- ⚠️ PRF = 1000 Hz
- ⚠️ v_max = ±3.75 m/s (会产生速度模糊)
- ✅ R_max_chirp = 150 km
- ✅ N_chirp = 500

---

## 设计权衡

FMCW雷达存在固有的权衡关系：

```
短 Tc (< 0.2 ms):
  ✅ 高 PRF → 高 v_max (适合高速目标)
  ❌ 低 R_max_chirp (但通常仍足够)
  ✅ 多 N_chirp → 高多普勒分辨率

长 Tc (> 1 ms):
  ❌ 低 PRF → 低 v_max (速度模糊严重)
  ✅ 高 R_max_chirp (适合超远距离)
  ❌ 少 N_chirp → 低多普勒分辨率
```

---

## 针对你的应用

**目标参数：**
- 距离：R0 = 200 m
- 速度：v = -20 m/s
- 需求：无速度模糊

**最佳选择：方案 1**

```matlab
Tc = 0.15e-3;    % 150 μs
```

理由：
1. v_max = 25 m/s > 20 m/s ✓ 无速度模糊
2. R_max_chirp = 22.5 km >> 200 m ✓ 距离范围充足
3. N_chirp = 3333 ✓ 充足的多普勒分辨率
4. 计算量适中

---

## 修改代码

在你的代码中，找到这一行：

```matlab
Tc = 27e-3;      % 单个chirp的时长 [s]
```

改为：

```matlab
Tc = 0.15e-3;    % 单个chirp的时长 [s] - 优化用于 v_max ≈ 25 m/s
```

---

## 预期结果

修改后，你应该看到：

```
========== RD Spectrum Parameters ==========
Chirp duration (Tc):        0.150 ms
Number of chirps:           3333
PRF:                        6666.67 Hz

--- Max Unambiguous ---
Max unambiguous velocity:   25.00 m/s  ✓

--- Body Theoretical Values ---
Body velocity (v):          -20.00 m/s  [No ambiguity!]
Body Doppler (fd):          -2666.67 Hz (-2.667 kHz)
```

**RD谱上的红色标记将精确对齐峰值！** 🎯

---

## 总结

- ❌ **错误**：增大 Tc (27 ms) → v_max 变小 (0.14 m/s)
- ✅ **正确**：减小 Tc (0.15 ms) → v_max 变大 (25 m/s)

**记住：Tc 越小，v_max 越大！**
