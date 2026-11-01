# 最小算例与检查

九个程序各有一个最小算例目录，保存输入、配置和迁移前参考科学结果。`cases.json` 记录入口与参考文件。`run_case.py` 为每次运行创建独立目录，解析配置中的相对输入路径，并逐字节比较输出；结果写入构建目录。

```sh
ctest --test-dir build/release --output-on-failure
python3 tests/run_case.py --case opn_nonsingular2 --binary build/release/bin/aster_opn_nonsingular2 --output build/manual/nonsingular2
```

配置型算例共七个，另外两项转换直接使用输入表。九个算例合计比较 16 份科学文件，覆盖轨道解、残差、候选表、预报星历和转换结果。参考结果来自迁移前当前工作区，已经包含 GLML 角度与随机种子修复。

`core/` 保存坐标、瞬转平、批处理、附加参数、RKF78、Adams 和偏心率公式检查。偏心率检查使用已保存的独立参考表；Adams 检查保留导师参考实现及来源说明。只用于旧、新 IOD 求解器相互比较的测试留在原工程。

`check_failures.py` 检查七列输入、有限值、空行和尾随空白、额外列拒绝、缺少配置、随机种子范围、程序版本以及后期报告写出失败。输入输出格式和错误行为随这些检查一起维护。
