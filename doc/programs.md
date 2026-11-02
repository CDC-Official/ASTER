# 程序使用

## 初轨、预报与批处理

八个配置型程序使用 `程序名 -I setting.txt`。配置中的输入、输出相对路径均从配置所在目录解释。省略配置参数时，程序显示用法。最小配置位于 `tests/` 中与程序同名的目录，去掉 `aster_` 前缀即可找到。

GL 与 GLML 使用地基 RA/Dec 观测。观测行依次为 UTC MJD、观测类型、目标编号、测站编号、赤经和赤纬；输入角度为度。测站文件保留已有大地坐标或 ITRS 坐标格式。计算中的状态采用库内归一化单位，文件输出单位由 `output_state_units` 决定。

GLML 先求基准 GL 解，再生成固定偏心率候选，经过边界处理后按最小 J 选取结果。`random_seed` 控制失败回退时的随机抽样，省略时为 0，可取 0 至 2147483647 的整数。种子写入控制台、报告和轨道解。`source_solution_flag` 保存来源求解状态：正数表示收敛轮次，-1 表示达到迭代上限，-2 表示数值失败。当前选择规则允许有限 J 的迭代上限候选参与比较，判断结果时同时查看该字段。

Cartesian 和 nonsingular2 程序处理地球目标，cislunar 程序处理地月预报。模型、积分器、步长、时段和输出文件按各自最小配置中的关键词设置。地月程序沿用当前 E-BCRS 动力学约定及地球非球形项近似。

批处理程序读取观测与初轨，迭代估计状态及所选附加参数。初轨文件接口和配置关键词保持现有格式，BODN 专用的关键词后缀也保持不变。已有输出的覆盖行为沿用当前批处理程序；日常实验应为每次任务准备独立目录。

控制台用于查看进度与摘要，report 保存科学过程，program log 保存事件和诊断。三者的开关和路径按配置分别设置。

## 观测生成

`aster_generate_observations` 读取星历和地面测站，输出理论 RA/Dec 弧段及可选真值；`tools/add_angular_noise.py` 生成可复现的东西/南北等标准差带噪样本。配置位于 `tests/generate_observations/setting.txt`，输入、筛选、单位和完整运行命令见 [观测生成](observations.md)。

## 状态转换

```sh
aster_convert_keplerian_to_cartesian input.txt output.txt
aster_convert_cartesian_to_keplerian input.txt output.txt m
aster_convert_osculating_to_mean -I setting.txt
```

根数输入每行是 `MJD_UTC a[km] e i[deg] OMEGA[deg] omega[deg] M[deg]`，M 为平近点角，适用椭圆轨道。转为 Cartesian 后，每行是 `MJD_UTC x[m] y[m] z[m] vx[m/s] vy[m/s] vz[m/s]`。反向转换的最后一个参数可用 `m` 或 `km` 指定输入单位，输出根数采用 km 和度。省略该参数会沿用交互选择，批量运行时应显式填写。

两项转换跳过空行，要求每个非空行恰好包含七个有限数值，并保留输入历元和参考系。尾随空白允许，额外列和非数值内容会报错。输出文件按现有行为覆盖。

瞬转平程序逐历元把密切 Kepler 根数转换为均根数，输入输出均为七列，采用 UTC MJD、km 和度。原有重力、日月、光压及迭代设置见 `tests/convert_osculating_to_mean/setting.txt`。

## 运行自己的资料

先用 `tests/run_case.py` 跑通对应入口，再复制该次生成的配置作为实验起点。将输入指向自己的资料，为输出填写独立路径，并检查资料覆盖的时段。配置参数的读取说明保留在各主程序中，涉及单位、参考系和归一化的计算说明位于对应头文件与实现。
