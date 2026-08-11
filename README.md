# Orbis

Orbis 是一个面向高精度 GNSS 定位与导航的 C++ 软件项目，提供传统滤波（FLT）和因子图优化（FGO）两条解算流程。项目面向 PPP、RTK 及 PPP/INS 紧耦合等应用，重点支持多系统、多频观测、模糊度固定和可重复的结果评估。

> 当前仓库处于持续开发阶段。配置项、输出格式和实验脚本可能随版本演进，请以当前源码和配置文件为准。

## 主要能力

- 基于 Ceres 的滑动窗口因子图优化、重线性化与边缘化。
- 基于滤波器的定位流程，可作为 FGO 的结果基线。
- PPP 原始观测（`RAW_ALL`）和无电离层组合（`IONO_FREE`）处理。
- GPS、Galileo、BDS-2/3、GLONASS 等卫星系统。
- 多频载波相位与码观测，包含五频模糊度、IFB 等参数路径。
- UPD 和 OSB 产品支持，以及浮点、参数反馈和约束反馈模式。
- PPP/INS 紧耦合和多源融合应用程序。
- 结果文件、日志和统计报告生成，以及 FGO/FLT 对照分析工具。

FGO 的 `RAW_MIX` 目前仅完成配置解析，原始观测因子图暂不支持该模式；请使用 `RAW_ALL` 或 `IONO_FREE`。

## 仓库结构

```text
Orbis/
├── src/
│   ├── LibGREAT/       FGO、GNSS、模糊度和融合算法库
│   ├── LibGnut/        GNSS 数据、产品和通用工具库
│   └── app/            应用程序入口
├── sample_data/        示例观测、导航和精密产品
├── scripts/            实验启动、监控和结果分析脚本
└── doc/                测试记录与设计说明
```

## 编译环境

- CMake 3.21 或更高版本。
- 支持 C++17 的编译器：MSVC、GCC 或 Clang。
- Eigen3、Ceres 2.x、GLFW、pugixml 和 spdlog。
- Windows 构建依赖可通过项目中的 vcpkg 配置准备；Linux/macOS 需要自行安装对应开发包。

## 编译

在仓库根目录执行：

```powershell
cmake -S src -B build -G "Visual Studio 18 2026" -A x64
cmake --build build --config Release --parallel
```

Linux 示例：

```bash
cmake -S src -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

主要应用目标为：

- `GREAT_PVT`：GNSS 定位，可通过 XML 选择 FLT 或 FGO。
- `GREAT_GINSFGO`：PPP/INS 紧耦合 FGO。
- `GREAT_MSF`：多源融合应用。

应用目标名称暂保留原有兼容名称，但项目对外名称为 Orbis。

## 运行

Orbis 使用 XML 配置文件描述输入数据、精密产品、观测组合、解算器、输出路径和日志设置。以 GNSS 定位程序为例：

```powershell
.\build\Bin\Release\GREAT_PVT.exe `
  -x .\path\to\config.xml
```

也可以在配置文件后追加节点或属性覆盖，例如指定测站和日志名：

```powershell
.\build\Bin\Release\GREAT_PVT.exe `
  -x .\path\to\config.xml `
  node:config:gen:rec=GODN `
  attr:config:outputs:log:name=GODN_run
```

FGO 配置中的关键选择包括：

- `estimator`：`FGO` 或 `FLT`。
- `obs_combin`：`RAW_ALL` 或 `IONO_FREE`。
- `fix_mode`：浮点或模糊度搜索固定。
- `ambiguity_fix_factor_enable`：`true` 保留模糊度固定约束因子并参与后续图优化；`false` 仅条件化当前图参数，不保留显式约束因子。
- `upd_mode`：按配置选择 UPD 或 OSB 产品。

启用 FGO 模糊度固定反馈时，应使用 `RAW_ALL` 并同时启用模糊度固定；`ambiguity_fix_factor_enable` 仅控制约束因子是否保留。窗口大小、GNSS 窗口大小和 INS 窗口大小不能超过源码编译时的容量上限。

## 实验与结果分析

`scripts/run_fgo_experiments.py` 可启动单个实验或并行实验矩阵，并持续记录状态、日志、性能和精度报告。示例：

```powershell
.\.venv\Scripts\python.exe scripts\run_fgo_experiments.py run `
  --manifest scripts\fgo_experiment_matrix.example.json `
  --run-id orbis_ff_30m
```

查看运行状态：

```powershell
.\.venv\Scripts\python.exe scripts\run_fgo_experiments.py status `
  --run-dir build\fgo_runs\orbis_ff_30m
```

实验目录通常包含：

- `state.json`：当前进度和最终状态。
- `logs/`：各进程的标准输出与错误输出。
- `analysis/`：单站 Markdown/JSON 统计报告。
- `report.md`：实验矩阵汇总。

分析脚本支持收敛时间、首次正确固定时间、固定率、三维 RMS、平均绝对误差、最大误差以及 P95/P99 误差分位数等指标。并行实验必须为每个进程配置独立的输出文件，避免结果相互覆盖。

## 输出文件

输出格式由 XML 配置决定，常见文件包括：

- `.fgo`：因子图优化结果。
- `.flt`：滤波结果及模糊度解算结果。
- `.log`：运行日志。

FGO 与 FLT 的结果文件应使用独立路径；分析脚本会检查文件是否由当前进程实际更新、历元是否连续以及坐标是否有限。

## 开发说明

- 修改观测模型、状态维度或边缘化接口后，应进行全量重新编译。
- RAW 因子涉及轨道、钟差、潮汐、光行时、姿态和偏导计算；修改并行化或缓存逻辑时必须与串行结果进行逐项对照。
- 五频测试应同时检查 IFB 参数是否由实际观测因子激活，避免不可观参数进入协方差和模糊度固定流程。
- 运行记录和异常历元分析统一保存在 `doc/` 或实验运行目录中。

## 第三方组件

Orbis 使用以下开源组件：

- [Eigen](https://eigen.tuxfamily.org)
- [Ceres Solver](https://ceres-solver.org)
- [GLFW](https://www.glfw.org)
- [pugixml](https://pugixml.org)
- [spdlog](https://github.com/gabime/spdlog)
- [Newmat](http://www.robertnz.net/nm_intro.htm)
- [G-Nut](http://www.pecny.cz)
- [PSINS](https://psins.org.cn)

各组件的版权和许可证以其随附文件及官方发布说明为准。

## 许可证与贡献

如果需要提交问题、改进代码或补充实验，请提供：

1. 使用的提交版本和操作系统；
2. XML 配置及输入产品的基本信息；
3. 可复现的命令行、日志和异常历元；
4. 对应的 `.fgo`/`.flt` 结果或统计报告。

## 项目链接

[https://github.com/thetengda/Orbis](https://github.com/thetengda/Orbis)

## 致谢

感谢 GREAT 项目及其开发团队提供的开源基础与支持。
