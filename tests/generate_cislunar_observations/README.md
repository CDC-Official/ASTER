# 地月单程测距最小算例

本目录包含两颗星的状态、两条接收事件、DE440 历表和运行输出。12 为导航星，901 为用户星，两个方向各生成一条观测。

## 输入

`aster.settings` 中的路径相对于配置文件。状态与事件历元均为真实 TDB，位置和速度为 E-BCRS。

| 文件 | 内容 |
|---|---|
| `input/navigation_state.txt` | 导航星的十条状态记录 |
| `input/user_state.txt` | 用户星的十条状态记录 |
| `input/receive_events.txt` | 真实接收时刻与发收编号 |
| `input/JPLEPH_C.440` | DE440 二进制历表 |

状态行依次为：`年 月 日 时 分 秒 x y z vx vy vz 钟差`，单位为 m、m/s、s。钟差为卫星钟读数减 TDB。接收事件行依次为：`年 月 日 时 分 秒 发射编号 接收编号`。

位置、速度和钟差逐分量进行最多十点的 Neville 插值。两端位置分别加各自事件的地球 BCRS 位置，求解光行时。引力时延使用 DE440 GM 和接收时刻的天体位置，参与天体为水星、金星、火星、木星、土星、天王星、海王星、冥王星、月球和太阳；太阳项采用 2GM/c² 修正。

`gravity_delay` 控制引力时延；`systematic_bias_m`、`noise_sigma_m` 控制伪距系统差和白噪声。该案例二者均为零。`noise_seed` 指定随机种子。

## 输出

- `output/observations.txt`：接收钟的年月日时分秒、发射编号、接收编号、伪距（m）。
- `output/components.csv`：真实发收 TDB、光行时、几何距离、引力时延、两端钟差、理论伪距、加入的误差和迭代情况。
- `output/run.log`：输入路径、观测计算记录和运行结果。
- `output/verification.json`：两条样例的数值核查结果。

接收钟读数包含接收钟差；分项文件中的接收时刻为 TDB。`reference.csv` 保存两条已验证基准的距离分项，`source.json` 记录状态来源、处理方式和输入校验值。

## 运行

从 ASTER 根目录执行：

```sh
cmake -S . -B build/cislunar_minimal -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build/cislunar_minimal --target aster_generate_cislunar_observations aster_test_cislunar_observation -j 3
build/cislunar_minimal/bin/aster_generate_cislunar_observations -I tests/generate_cislunar_observations/aster.settings
ctest --test-dir build/cislunar_minimal -R '^(core_cislunar_observation|cislunar_observation_runtime)$' --output-on-failure
```

本目录已保存一次运行结果。再次运行时，把 `directory_output` 改成新的目录，例如 `output_repeat`。CLion 选择 `aster_generate_cislunar_observations`，运行参数填写 `-I tests/generate_cislunar_observations/aster.settings`，工作目录为仓库根目录。

核心验证覆盖解析插值、亚微秒历元、引力时延、跨日发收事件和钟差符号；主程序验证覆盖双向观测、系统差、噪声复现、引力开关及错误追踪。

## 待查事项

DE440 常数统一后，600 秒近地预报算例的位置变化约 3.65 cm，留待单独检查。
