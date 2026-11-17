# GCRS/TT 核对算例

从仓库根目录运行专用入口；Windows 构建和 GM 说明见[交付文档](../../doc/18_gcrs-tt-orbit-check.md)。所有配置项均在 `setting.txt` 中显式给出。

- `reference_cartesian.csv`：60 s 步长的完整十天参考轨迹，14401 行，含初态；位置 km、速度 km/s。
- `reference_summary.txt`：该次运行的实际 GM、模型、初末状态和构建信息。
- `reference_evidence.json`：源码提交、平台、输入及输出 SHA-256、步长复算和检查记录。

双方按 `elapsed_TT_s` 和绝对历元对齐，再计算 `sqrt(dx²+dy²+dz²)`；CSV 位置差乘 1000 后单位为 m。参考是在 macOS 的 SYSTEM 后端生成，Windows 结果可据此做数值比较。

自动检查：以 `-DBUILD_TESTING=ON -DASTER_MATH_BACKEND=SYSTEM` 配置并构建后，执行：

```sh
ctest --test-dir build/client -C Release -R "^(core_gcrs_tt_check|gcrs_tt_check_runtime)$" --output-on-failure
```

`eop_two_days.txt` 是检查插值边界用的合成资料；默认十天算例使用 `setting.txt` 指定的公共 EOP 文件。
