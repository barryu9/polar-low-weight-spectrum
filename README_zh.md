# Polar / PAC 低重量谱：论文 Algorithm 1 开源实现

本仓库包含独立的 **Python 原型、纯 MATLAB 原型和 C++17 kernel**。
三者均实现论文 Algorithm 1 的完整流程，并精确计算截至 D 的低重量谱。

## 直接运行

在仓库根目录执行：

~~~shell
python examples/python/run_spectrum.py --n 64 --kind polar -D 8
python examples/python/run_spectrum.py --n 128 --kind pac --mode puncture -D 8
python tools/build_cpp.py
python examples/python/run_spectrum.py --n 256 --kind pac --backend cpp --mode shorten -D 8 --repeat 7
~~~

MATLAB 中：

~~~matlab
addpath('examples/matlab');
r = run_spectrum(64,'polar','none',8);
r = run_spectrum(128,'pac','puncture',8);
r = run_spectrum(256,'pac','shorten',8);
~~~

MATLAB 示例直接调用 MATLAB 实现，使用自带的精确大整数运算，
不依赖 Python、C++、Java 或符号工具箱。Python 示例可选择 python/cpp。

## 只求 dmin 和 A_dmin

Python 示例增加 **--dmin-only**，不需要指定 D，可选 Python 或 C++：

~~~shell
python examples/python/run_spectrum.py --n 64 --kind pac --mode mixed --dmin-only
python examples/python/run_spectrum.py --n 128 --kind polar --backend cpp --dmin-only
~~~

MATLAB 把 D 参数设为 **'auto'**：

~~~matlab
r = run_spectrum(64,'pac','mixed','auto');
disp(r.dmin);
disp(r.A_dmin);
~~~

算法从可证明的重量下界开始，逐层寻找第一个非零重量，
完整统计该层重数后立即停止，不计算更高重量谱。
预处理只做一次，搜索过程共用 max_work 总预算，计时包含整个搜索。
缩短、凿孔后会按实际发送码重新识别 dmin；不假定 dmin=8。
只有全零码字时，dmin 无定义，A_dmin=0。
Python 的 --dmin-only 与 -D 互斥；原有数值 D 模式仍可使用。
旧版 C++ 库需要重新执行 python tools/build_cpp.py。

## 已固定的约定

- 六组母码：Polar 和 PAC 各有 (64,32)、(128,64)、(256,128)。
- GA 的 2.0 dB 指母码率 R=1/2 下的 Eb/N0；两类码共用 GA 信息位。
- PAC taps 为 [1,0,1,1,0,1,1]，按 u_i = Σ g_l v_(i-l) 编码，初始状态全零。
- 信息位、缩短和凿孔坐标均已硬编码；另提供 JSON 和二进制选择矩阵 CSV。
- D 默认 8，可自由指定。增大 D 可能显著增加搜索量。
- 矩阵针对输出码字坐标。缩短先要求这些位为零，再删除；凿孔直接删除。
- 结果统计不同的发送码字，凿孔产生的重复映射会被消除。

默认 (64,32) Polar、D=8 的结果是 A_0=1、A_8=664。
对应 PAC 的结果是 A_0=1、A_8=472。

## 与论文的关系

实现保留无噪声阈值 coset 枚举及截断 WEF 递归，并加入仿射族计数。
内部需要继续递归的路径仍以实际 (u,x) 对流式产生；可直接结算的族
以精确乘数一次计入。局部正代价偏差可选枚举或 syndrome DP。

[算法、公式和复杂度对应](docs/algorithm.md)说明了论文与代码的映射。
[构造约定及 API](docs/conventions.md)给出坐标、GA、PAC 和矩阵定义。
[核验记录](docs/validation.md)说明实测范围。

每次返回的谱均为截止 D 内的精确结果。超过 max_work 会明确报错，
不会把部分结果标成完整谱。计时包括预处理和谱计算，不包括固定 profile
加载、编译或 MATLAB 启动；C++ 另报告包含绑定开销的时间。

源码采用 MIT 许可证。仓库可独立构建，不依赖原有 SCL 工程。
