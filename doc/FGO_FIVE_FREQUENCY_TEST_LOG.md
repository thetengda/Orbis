# FGO 五频 RAW/OSB 开发与测试记录

## 任务范围

- UPD 算例：`sample_data/PPPFLT_2023305`
- OSB 算例：`sample_data/PPPFLT_2023305_OSB`
- 首轮测试时段：2023-11-01 00:00:00--02:00:00（30 s）
- 最终回归时段：2023-11-01 全天
- 对比基线：相同频率、系统、观测和固定设置下的 PPP 滤波解

## 阶段 0：既有修改基线

- 基线分支：`feature/osb-support`
- 基线提交：`a9be762 feat(fgo): checkpoint RAW PPP graph implementation`
- 提交范围：23 个文件，新增 3,213 行、删除 177 行；未包含 `.vscode`、样例数据和既有结果。
- 编译命令：`cmake --build build --config RelWithDebInfo --target GREAT_PVTFGO -- /m`
- 编译结果：通过；生成 `build_Windows/Bin/RelWithDebInfo/GREAT_PVTFGO.exe`。

## 现有观测模型与数据流审计

### 观测模型

`RAW_ALL` 对每颗卫星的每个已配置且有效的频点分别建立码和相位观测。公共精密模型 `t_gprecisebiasFGO::cmb_equ()` 负责卫星轨钟、天线、潮汐、相位缠绕、对流层、偏差等改正和线性化；RAW 因子只把 Ceres 参数块映射为 `t_gallpar`，读取模型返回的残差、权和雅可比，避免在因子内重复施加观测改正。

主要状态与符号如下：

- 公共历元状态：`CRD_X/Y/Z`、`CLK`、`TRP`、各系统 `ISB`。
- 卫星历元状态：参考频率上的斜向电离层 `SION`。
- 卫星频点弧段状态：`AMB_L1`--`AMB_L5`，以卫星、频率和弧段 ID 管理。
- 码方程的电离层项为正；相位方程的电离层项为负；非参考频率按 `(f1/fi)^2` 缩放。
- 相位模糊度以米为单位进入因子，与滤波 RAW 参数类型一致。

### 单历元到滑窗图的数据流

1. `processBatch()` 逐历元调用 `processWindow()`，先执行周跳记录、时间更新和公共 PPP 数据预处理。
2. `_get_initial_value()` 从滤波式 PPP 参数初始化当前 FGO 节点，并由 `t_gambRAW_manager` 创建/复用逐频模糊度弧段；单频周跳只重建受影响频点，卫星级 ID 和其他频点保持连续。
3. `_combine_RAW()` 对配置的 `FREQ_1`--`FREQ_5` 选择精确 `GOBS`，生成一条观测对应一个 `RAWEquMsg`；消息保留历元、卫星、频率、波段、信号、SION ID 和模糊度弧段 ID。
4. `_optimization_PPP_RAW()` 注册公共状态、被观测引用的 SION 和模糊度参数，添加码/相位因子、SION/ISB 随机游走、初值先验和上一窗口边缘化先验，然后由 Ceres 重优化。
5. `_posteriori_test_PPP_RAW()` 生成与 `t_gallpar` 顺序一致的后验雅可比、残差、权和协方差；异常值检测按历元、卫星、信号、频率反查并迭代剔除。
6. `_pre_amb_resolution()` 将 FGO 后验装入既有滤波/模糊度固定接口；当前 `_amb_resolution()` 仅计算条件固定解并写 FLT 输出，固定结果尚未反馈到 Ceres 图。
7. `_marginalization_PPP_RAW()` 将最老节点观测、过程因子和已有先验线性化，消去该节点的公共状态与 SION；仍被窗口使用的跨历元模糊度保留，并进行参数地址平移。
8. `_slide_window()` 移除最老观测和节点，移动公共/SION 状态，并同步 RAW 模糊度管理器。

### OSB 数据流

观测读取后由 `t_gsatdata::apply_bias()` 对每个精确信号施加 Bias-SINEX OSB，并记录 `osb_corrected(obs)`。`_combine_RAW()` 在 `upd_mode=OSB` 时只保留已成功改正的码/相位观测；`_amb_resolution()` 还会检查参与固定的相位及配套码 OSB 完整性。RAW 因子本身不再次施加 OSB。

### 审计发现与待验证项

- 已有 FGO 两小时配置实际是双频；五频测试需要使用 GPS 三频、GLO 双频、GAL/BDS 五频的 FF 配置。
- 现有构造日志明确说明固定约束尚未反馈到图，正是阶段 4 要补齐的功能。
- RAW 图目前对每个窗口节点重复添加坐标、钟差、对流层、ISB 和全部活动模糊度初值先验；需通过与滤波基线的数值比较判断这些先验是否导致收敛或精度偏差。
- 固定入口依赖 `_all_para_win` 与后验矩阵列严格同序；五频下需重点检查参数数量、弧段映射和协方差可逆性。
- OSB 模式会按单观测丢弃缺少 OSB 的信号；需统计保留频点、可固定模糊度数和输出连续性。

## 测试记录

后续每次运行按“提交/可执行文件、配置、命令、耗时、退出码、日志诊断、输出完整性、定位统计、固定统计、结论”记录于本节。

### 浮点阶段：诊断运行与 IFCB 修正

不计入验收的运行尝试也保留如下，避免把不完整或输出冲突的结果误作基线：

- 首次按双站顺序运行在 20 分钟工具时限内仅推进约 24%，因预计总耗时过长而停止；随后一次长时运行发现单进程工作集接近 1 GB，改为分站并行。
- 第一次分站并行使用了错误的 XML 命令行覆盖层级，两个进程都处理 GODN 且争用同名输出，结果作废。用生成配置验证后，确认正确站点覆盖形式为 `node:config:gen:rec=...`，日志属性覆盖为 `attr:config:outputs:log:name=...`。
- 命令行对含时分秒的 `<end>` 文本覆盖未生效，一次拟运行 10 分钟的进程仍读取两小时终止时间并写默认输出；发现后立即停止。后续短段测试直接临时修改测试 XML，并在运行后用补丁恢复，未把临时值留入配置。

- 候选基线：`a9be762` 构建的 `GREAT_PVTFGO.exe`。
- 配置：新增 `GREAT_PPPFGO_kinematic_FF_Float_2h.xml`，GPS 3 频、GLO 2 频、GAL/BDS 5 频，`RAW_ALL`、`fix_mode=NO`、`upd_mode=UPD`。
- 诊断运行：GODN、HARB 分站并行；分别完成 61/239 和 72/239 个有效历元后主动停止。停止前无崩溃、NaN/Inf、Ceres failure 或协方差失败；单进程工作集约 0.98--1.00 GB。
- 发现：滤波 `t_gcombALL` 在 `ifcb_model=COR`、GPS `FREQ_3` 且无 phase OSB 时会把 IFCB 加到 RAW 相位残差；FGO 主程序却仅在 `frequency == 3` 时向处理数据注册 IFCB，五频运行未注册，RAW 因子也没有该观测域修正。因此五频 FGO 与滤波观测模型不一致。
- 修正：主程序在 `frequency >= 3` 且 IFCB 数据存在时注册 IFCB；`_combine_RAW()` 计算 GPS `FREQ_3` 相位 IFCB，并随 `RAWEquMsg::additive_correction` 固化，使优化、后验检验和边缘化使用同一修正。存在 phase OSB 时不叠加 IFCB。
- 编译：`cmake --build build --config RelWithDebInfo --target GREAT_PVTFGO -- /m`，通过。
- 修正后 10 分钟 GODN 试跑：19/19 个有效历元，退出码 0，墙钟 120.692 s；末历元坐标 `(1130760.6430, -4831298.6937, 3994155.1425) m`；无 NaN/Inf、Ceres failure、协方差失败或异常中断。
- 性能试验（未采用）：用全图雅可比的稠密法方程代替 Ceres 全协方差，两个 10 分钟对照分别为 129.292 s 和 108.530 s，输出与原方法逐行相同，但相对 120.692 s 的稳定收益有限且增加尺度/秩风险，故已恢复原 Ceres `SPARSE_QR` 协方差路径，不纳入代码 diff。
- 完整两小时验收：GODN、HARB 分站并行，两个进程退出码均为 0；总墙钟 1,784.897 s，程序内处理耗时分别为 1,768.93 s 和 1,548.91 s。可执行文件 SHA-256 为 `B5940C90DD02F1403092EE86BB99BA08546A85B4461A45953D902CD990D9B5D9`。
- 输出完整性：两站 FGO/FLT 各 239 个历元，SOW 均为 259230--266370，全部相邻间隔 30 s；坐标无 NaN/Inf；全部为 Float。每站 FGO 与 FLT 的 239 个同历元坐标完全一致（输出精度下最大三维差 0）。
- 末历元：GODN 为 `(1130760.7021, -4831298.6690, 3994155.2014) m`，相对既有同日五频滤波参考末值三维差 0.0132 m；HARB 为 `(5084657.5947, 2670325.4524, -2768480.8300) m`，三维差 0.0443 m。
- 日志诊断：关键错误匹配 0，RAW 异常值剔除 0；未出现 Ceres failure、协方差不可用、NaN、异常、节点回滚或可用卫星不足。每站存在 24 条主程序既有的误标级别消息 `path is not file (skipped)`（代码实际随后正常打开 `file://` 输入），以及 510 条 ATX 对未支持 R04/R06 频率码的既有解码消息；均发生在输入初始化且未影响本次观测频率、历元连续性或结果。
- 结论：五频 RAW 浮点流程通过两站两小时验收；GPS `FREQ_3` IFCB 与滤波模型已对齐，可提交阶段 1。

### 固定阶段：双差约束诊断

- 输入提交：`5f98ca7 fix(fgo): align five-frequency RAW float IFCB model`。
- 配置：新增 `GREAT_PPPFGO_kinematic_FF_Fixed_2h.xml`，频率/系统与浮点配置一致，改为 `fix_mode=SEARCH`、`upd_mode=UPD`。
- 修正前 GODN 10 分钟：退出码 0、墙钟 112.907 s；FGO 19 个历元，FLT 从 00:01 起输出 18 个历元，其中 Fixed 6、Float 12，首次固定 SOW 259260。以既有五频滤波末值为临时参考，固定历元三维差均值 0.3397 m、最大值 0.5903 m；后段浮点收敛到约 0.06--0.08 m。该短段说明固定流程可执行，但早期 FGO 浮点状态尚未收敛，不能仅凭 ratio 将其视作最终定位验收。
- 代码审计发现：`t_gambiguity::_addFixConstraint()` 把双差系数 `Ba/Bb` 定义在循环外；GLONASS 分支将其改为逐星波长系数后，若迭代顺序再进入其他系统/频率，旧系数可能泄漏。现已改为每条双差方程重置 `1/-1`。同时移除固定成功后写死到 `D:\GREAT-FGO\amb_WL_NL_debug.txt` 的调试输出，并把无条件 stdout ratio 改为 spdlog debug。
- 修正后同段：退出码 0、墙钟 114.081 s；18 个 FLT 历元、6 个 Fixed，与修正前逐行相同，说明本算例当前双差排序未触发系数泄漏；stdout 中孤立 ratio 数字由 2 行/次降为 0。该修正消除了排序相关隐患和外部文件副作用，不改变本次已正确的数值路径。
- 完整两小时验收：GODN、HARB 分站并行，两个进程退出码均为 0；总墙钟 1,887.347 s，程序内处理耗时分别为 1,871.55 s 和 1,633.63 s。动态库 `LibGREATrd.dll` SHA-256 为 `7D541CFECC3675D327CB17B5CDA018770E5D4BE78FA51B371EBCDE2A4DA7B8F8`。
- 输出完整性：两站 FGO 各 239 个历元（SOW 259230--266370），FLT 各 238 个历元（SOW 259260--266370）；全部相邻间隔 30 s，坐标无 NaN/Inf。当前尚未向图反馈固定结果，因此 FGO 文件按设计仍全部标记为 Float。
- 固定统计：GODN 条件固定 56/238（23.5%），首次固定 SOW 259260，共 10 段、最长连续 32 个历元，固定 ratio 范围 2.06--50.58；HARB 条件固定 87/238（36.6%），首次固定 SOW 259410，共 24 段、最长连续 18 个历元，固定 ratio 范围 2.03--32.46。
- 定位统计：以既有同日五频滤波末值作稳定参考，GODN 条件固定历元三维误差均值/RMS/最大值为 0.0804/0.1379/0.5903 m，HARB 为 0.0485/0.0839/0.4517 m；大值集中于初始化阶段。去除前 30 分钟后，不分固定状态的三维误差均值/RMS/最大值分别为 GODN 0.0312/0.0354/0.1188 m、HARB 0.0326/0.0389/0.0947 m。条件固定相对同历元 FGO 浮点坐标的三维改变量均值/最大值分别为 GODN 0.0600/0.5857 m、HARB 0.0377/0.4070 m。
- 日志诊断：主日志、stdout、stderr 中 Ceres/协方差失败、NaN/Inf、异常退出、求解失败和可用卫星不足等关键模式合计命中 0；移除的硬编码文件输出未再产生，无条件 ratio stdout 行为也已消失。
- 当时结论：五频 RAW 的 NL 搜索与条件固定主流程可执行；后续用正确 Bias-SINEX 解码重跑时发现 EWL24/EWL25 协方差和整数回写分支缺失，因此本节数值仅保留为修正前基线，最终五频固定验收以后续回归为准。

### OSB 阶段：Bias-SINEX 解码入口修正

- 输入提交：`298ac7e fix(ambiguity): stabilize five-frequency RAW fixing`。
- 配置：新增 `sample_data/PPPFLT_2023305_OSB/xml/GREAT_PPPFGO_kinematic_FF_Fixed_2h.xml`，使用 WUM 快速轨道、钟差和 `WUM0MGXRAP_20233050000_01D_01D_OSB.BIA`，五频 `RAW_ALL`、`fix_mode=SEARCH`、`upd_mode=OSB`，不配置 UPD 或 IFCB。
- 无效命令记录：首次冒烟命令遗漏 `-x`，程序在 0.181 s 内由参数解析返回 1，未开始数据读取；更正命令后才计入流程诊断。
- 修正前 GODN 10 分钟：正确命令退出码 0，但 19 个候选历元全部报告 `combine PPP observations wrong`，FGO/FLT 都是 0 行。新增载入摘要后确认 `Loaded bias products: 0 analysis center(s) [], phase OSB=false`，说明不是观测缺失，而是 Bias-SINEX 没有进入偏差容器。
- 根因：`GREAT_PVTFGO` 对通用 `<bias>` 输入错误地实例化 `t_biasinex`；同仓库滤波程序 `GREAT_PVT` 对该输入使用支持 WUM/CAS 等 Bias-SINEX OSB 记录的 `t_biabernese`。现已把 FGO 入口与滤波入口对齐；`<biasinex>` 的原有 `t_biasinex` 分支保持不变。另增加一条 INFO 摘要，报告实际载入的分析中心及是否包含绝对相位 OSB，便于防止“退出码 0、结果 0 行”的静默退化。
- 修正后 GODN 10 分钟：退出码 0、墙钟 71.131 s；载入 1 个分析中心 `[WHU_A]` 且 `phase OSB=true`；FGO 19/19、FLT 18/18 个历元，FLT 中 Fixed 15、Float 3；缺 OSB 信号跳过警告 0，关键错误 0。
- 解码修正后的首轮两小时运行（在后续 Win64 地址键和 EWL24/EWL25 修正之前）：GODN、HARB 分站并行，两个进程退出码均为 0；总墙钟 870.967 s，程序内处理耗时分别为 846.608 s 和 860.841 s。可执行文件 SHA-256 为 `2A47E7D6027549DA47A15A306D5025EE9139FA7B7B28B9B0E75DF844CA1C6877`。
- 输出完整性：两站 FGO 各 239 个历元（SOW 259230--266370），FLT 各 238 个历元（SOW 259260--266370）；相邻间隔均为 30 s，坐标无 NaN/Inf。两站日志均且仅有 1 条 `[WHU_A], phase OSB=true` 载入摘要，`OSB mode skipped` 为 0。
- 固定统计：GODN 条件固定 235/238（98.7%），首次固定 SOW 259290，固定 ratio 范围 1.51--18.51；3 个 Float 均在前 3 分钟。HARB 条件固定 232/238（97.5%），首次固定 SOW 259260，固定 ratio 范围 1.50--81.07；6 个 Float 分布在 SOW 260310--260340、262680--262710 和 265650--265680 的短弧段。
- 定位统计：以既有两小时 OSB 滤波末值作稳定参考，去除前 30 分钟后 GODN 条件解三维均值/RMS/最大值为 0.0159/0.0191/0.0540 m，HARB 为 0.0124/0.0145/0.0405 m；相对既有同历元 OSB 滤波结果的后 90 分钟三维 RMS 分别为 0.0206 m 和 0.0133 m。全段差值 RMS 受初始化影响，分别为 0.0856 m 和 0.0481 m；HARB 前两历元虽通过低 ratio 固定但误差约 0.53--0.59 m，阶段 4 反馈测试需重点检查错误固定的传播风险。
- 日志诊断：主日志、stdout、stderr 中 Ceres/协方差失败、NaN/Inf、异常退出、求解失败、可用卫星不足和观测组合失败等关键模式合计命中 0。
- 首轮结论：五频 RAW OSB 读取、逐信号改正和筛选已得到验证；该轮固定结果随后被 EWL24/EWL25 缺失分支的诊断取代，不作为最终固定验收。修复恢复了 FGO 对通用 Bias-SINEX `<bias>` 输入的既有语义，`<biasinex>` 分支保持不变。

### 五频 EWL24/EWL25 与 Win64 参数地址键修正

- 正确的 CAS Bias-SINEX 解码使 UPD 算例实际载入 DCB 后，完整运行曾非确定性触发 Ceres `Map key not found`。根因是 Windows LLP64 上 `long` 为 32 位，而 GNSS 后验和边缘化代码把 64 位参数指针转换为 `long` 作为 map key，地址高位被截断并可能碰撞；短段是否触发取决于进程内存布局。
- 首轮只修改容器类型后，仍有一个局部 `const long address` 使后验描述符查找全部失败；该次进程虽然退出 0，但所有历元均被拒绝，明确作废。补齐局部变量后，2 分钟诊断恢复有效 posterior、协方差以及 FGO/FLT 输出。
- 修正范围随后扩展到 `GNSSInfo`、`MarginalizationInfo` 及其 PVT/GINS 调用链，地址键统一为 `std::uintptr_t`。XML、CLI 和结果文件格式未变；但导出 C++ 类的方法签名及 public/protected 容器类型发生 ABI/源码接口变化，必须同步重编 `LibGREAT` 与所有消费者，不能与旧二进制混用。
- 地址键修正后的 UPD 两站 2 小时兼容回归（尚未包含下述 EWL 修正）：GODN/HARB 均退出 0，FGO 各 239 行、FLT 各 238 行且 30 s 连续，无 NaN/Inf、Ceres/地址键/后验描述符失败。GODN 固定 217/238、首次固定 SOW 259470；HARB 固定 222/238、首次固定 SOW 259350。去除前 30 分钟后，FLT 相对同历元滤波基线的三维 RMS 分别为 0.0157 m 和 0.0186 m。
- 该回归的详细日志每个处理历元都出现 `_prepareCovarianceWL ... prepare Double-Difference covariance is Wrong`：GODN 420 条、HARB 442 条。审计确认这不是可忽略的协方差偶发错误，而是 `_prepareCovarianceWL()` 只实现 WL/EWL，EWL24/EWL25 进入时波长保持 0、值向量为空；上层又忽略中间组合返回值，所以进程仍可完成，但 L4/L5 的 LAMBDA、整数保存和宽巷约束实际缺失。
- 修正为四种宽巷组合的统一显式映射：WL=`L1-L2`、EWL=`L2-L3`、EWL24=`L2-L4`、EWL25=`L2-L5`；补齐 EWL24/EWL25 的协方差、部分固定标志清理和整数回写，并在未知模式、端点不足、无效模糊度、缺波长、非正/非有限波长、非有限协方差及维数不一致时安全返回且输出准确诊断。现有 `_addFixConstraintWL()` 已能把两种组合继续写入滤波虚拟观测。
- 第一轮 EWL 修正的 GODN 10 分钟门禁：UPD/OSB 并行，退出码均为 0，程序内耗时 94.750 s/63.628 s；两者 FGO 19 行、FLT 18 行，分别固定 12/18 与 17/18。旧协方差错误、`_selectAmb Wrong`、未知组合、非正波长、维数不一致、地址键和 Ceres failure 均为 0。
- 随后启动的 UPD/OSB × GODN/HARB 两小时并行对照均退出 0、总墙钟 1,350.479 s，但代码复核发现既有 WL/EWL 协方差只在两条 DD 的首颗卫星 PRN 数字相同时填非对角元：它会漏掉共享第二端点/参考星的相关性、误把跨系统同号卫星关联，并对 GLONASS FDMA 复用第一条 DD 的波长。该轮因此只作为“EWL24/EWL25 已执行、相关性模型尚未修正”的中间诊断，不计入最终验收。
- 协方差最终修正：按每条 DD 的 `[satA@band1, satB@band1, satA@band2, satB@band2]` 四个真实参数构造 `[+1/lambda_A1, -1/lambda_B1, -1/lambda_A2, +1/lambda_B2]`，逐条计算全部上三角 `c_i^T Qx c_j`；GLONASS 按卫星取 FDMA 波长，其他系统按系统取波长，不再按 PRN 猜测相关性。`Qx` 为单位权协因数，最后统一乘 `sigma0^2`；同时要求 `sigma0>0`，校验四端点索引/波长/协方差有限性，并显式拒绝未实现的 GLONASS 非 WL 约束，避免零除。
- 参数 shift 安全性补强：两个 `getParameterBlocks()` 不再用 `operator[]` 为缺失 shift 静默插入空指针，而是把 prior 标为无效并安全返回；IF-PPP 的旧先验 ISB0 在星座刚丢失时也无条件进入 drop set。RAW 路径原有 retained-block fallback 保持不变。
- WL/EWL 完整协方差与地址键修正后的构建：`cmake --build build --config Release --target GREAT_PVTFGO -j 4` 与同参数 `GREAT_GINSFGO` 均通过，确认 PVT 与 GINS 两个主要消费者可针对新的边缘化 API 完整链接。
- 完整协方差版 GODN 10 分钟门禁：UPD/OSB 并行，退出码均为 0，程序内耗时 95.142 s/63.737 s；两者 FGO 19 行、FLT 18 行，分别固定 12/18 与 18/18。末历元分别为 `(1130760.6830, -4831298.6350, 3994155.1548) m` 和 `(1130760.6991, -4831298.6994, 3994155.2197) m`。EWL/波长/索引/协方差/`sigma0`/维数、未知组合、地址键和 Ceres failure 均为 0；测试后两份 XML 的结束时间均恢复为 `01:59:30`。
- 随后一次拟作为最终验收的四组两小时运行，在 585.228 s 时由静态复核发现 NL 协方差阻断问题后主动停止；四个子进程退出码均为 -1，截断结果全部作废。问题是 `_prepareCovariance()` 对 GLONASS 两端仍使用同一 `factor`，与浮点 NL 值和最终约束所用的逐星波长不一致，且 NL 未像 WL/EWL 一样乘 `sigma0^2`。
- NL 最终修正：每条 DD 按两个真实端点构造 `[+1/lambda_A, -1/lambda_B]`，RAW 按 `AMB_L1`--`AMB_L5` 选择波段、IF 使用 `NL`，GLONASS 逐星取波长、其他系统按系统取波长；计算全部上三角 `c_i^T Qx c_j * sigma0^2`，并加入候选端点、波长、参数索引、有限性和维数检查。UPD/OSB 产品误差继续遵循既有“确定性改正”策略；当前产品接口只保留最后一次查询的 sigma，未伪造不完整的跨 DD 产品协方差。
- 同时把 `_selectAmb()` 的 `ndef==0` 从误导性 ERROR 改为 DEBUG“无独立候选”，仍返回 -1 并由上层按 Float 回退，数值语义不变。最终 `GREAT_PVTFGO`、`GREAT_GINSFGO` Release 重编均通过。
- NL 完整协方差版四组 20 分钟门禁：UPD/OSB × GODN/HARB 全部退出 0、总墙钟 240.156 s；每组 FGO 39 行、FLT 38 行且连续。固定数为 UPD GODN 32/38、UPD HARB 25/38、OSB GODN 38/38、OSB HARB 32/38；末历元坐标分别为 `(1130760.6860, -4831298.6496, 3994155.1696)`、`(5084657.6094, 2670325.4740, -2768480.8424)`、`(1130760.6920, -4831298.6671, 3994155.1965)`、`(5084657.6075, 2670325.4800, -2768480.8448) m`。NL/WL/EWL、波长/索引/维数、地址键、未知组合和 Ceres failure 均为 0；门禁后 XML 恢复为两小时。
- 可复现运行命令（分别在两个算例目录执行；`STATION` 取 `GODN` 或 `HARB`，并为并发进程指定不同 `LOG`）：`..\..\build_Windows\Bin\Release\GREAT_PVTFGO.exe -x .\xml\GREAT_PPPFGO_kinematic_FF_Fixed_2h.xml node:config:gen:rec=STATION attr:config:outputs:log:name=LOG`。UPD 算例目录为 `sample_data/PPPFLT_2023305`，OSB 为 `sample_data/PPPFLT_2023305_OSB`；20 分钟门禁只临时把 XML `<end>` 改为 `00:19:30`，运行后立即用补丁恢复。

### NL/WL 完整协方差后的最终两小时回归

- 最终配置均为 `00:00:00--01:59:30`，UPD/OSB × GODN/HARB 四个进程并行运行；全部退出码为 0，总墙钟 1,350.402 s。OSB 两站程序内处理耗时为 884/906 s；UPD 两站分别约在总墙钟 1,125/1,350 s 时结束。
- 输出完整性：四组 FGO 均为 239 历元、SOW 259230--266370；FLT 均为 238 历元、SOW 259260--266370；全部相邻间隔严格为 30 s。所有坐标和 ratio 均为有限数，原文 NaN/Inf 扫描为 0。固定结果在本阶段仍只写入 FLT，因此四组 FGO 的 239 个历元均标记为 Float。
- UPD 固定统计：GODN 为 Fixed 219/238（92.02%）、首次固定 SOW 259260，Fixed ratio 的 min/median/mean/p95/max 为 2.01/6.80/8.9414/21.50/58.00；HARB 为 Fixed 218/238（91.60%）、首次固定 SOW 259320，ratio 为 2.00/4.455/5.3961/11.1005/25.23。
- OSB 固定统计：GODN 为 Fixed 238/238（100%）、首次固定 SOW 259260，ratio 的 min/median/mean/max 为 1.74/7.275/8.0158/28.74；HARB 为 Fixed 228/238（95.80%）、首次固定 SOW 259260，Fixed ratio 为 1.58/19.565/22.3166/80.01。HARB 的 10 个 Float 历元为 SOW 259830、260040、260070、260100、260310、260340、262680、262710、265650、265680。
- 与同历元滤波条件解比较：按 XML 名义起点去除前 30 分钟（SOW >= 261000，180 历元），FGO 与 FLT 的三维 ECEF 差 RMS/均值/最大值分别为 UPD GODN 0.02358/0.02018/0.04485 m、UPD HARB 0.03412/0.03159/0.06415 m、OSB GODN 0.02038/0.01883/0.03972 m、OSB HARB 0.03205/0.02606/0.08478 m。全段 RMS 因初始化及条件固定改变量分别为 0.06606、0.04417、0.12561、0.19134 m。
- 以“从该历元到结束一直不超过阈值”为口径，FGO/FLT 三维差永久进入 0.1/0.05 m 的 SOW 分别为：UPD GODN 259440/259950，UPD HARB 259830/263850，OSB GODN 259890/260190，OSB HARB 260580/261570。该指标衡量因子图浮点解和滤波条件解的一致性，不等同于相对外部真值的定位收敛时间。
- 末历元 FGO/FLT 坐标：UPD GODN 为 `(1130760.6928,-4831298.6698,3994155.1898)` / `(1130760.6953,-4831298.6713,3994155.1980)` m，UPD HARB 为 `(5084657.6104,2670325.4782,-2768480.8330)` / `(5084657.5964,2670325.4731,-2768480.8274)` m；OSB GODN 为 `(1130760.7058,-4831298.6625,3994155.1943)` / `(1130760.6894,-4831298.6613,3994155.1933)` m，OSB HARB 为 `(5084657.6022,2670325.4684,-2768480.8285)` / `(5084657.6146,2670325.4809,-2768480.8414)` m。
- 产品载入摘要：UPD 两站均且仅有 1 条 `1 analysis center(s) [CAS_R], phase OSB=false`；OSB 两站均且仅有 1 条 `1 analysis center(s) [WHU_A], phase OSB=true`。OSB 保留可改正信号并逐历元跳过缺 OSB 信号：GODN 共 7,017 条（每历元 20--36），HARB 共 5,367 条（每历元 21--24）；输出仍保持连续且固定率如上。
- 硬错误扫描：NL/WL/EWL 协方差、未知组合、非正波长、参数索引、非有限协方差、`sigma0`、维数、`Map key not found`、`_selectAmb Wrong`、Ceres failure、观测组合失败均为 0。既有非致命噪声包括每组 239 条空的 `t_gpvtfgo` ERROR 级标记，以及 ATX 对未使用频率码和输入路径跳过的初始化消息；不计为求解失败。
- 最终构建复核：沙箱内首次增量构建仅因 MSBuild `FileTracker` 权限返回 `E_ACCESSDENIED`，按相同命令在允许环境重跑后，`GREAT_PVTFGO` 与 `GREAT_GINSFGO` Release 均成功。SHA-256：`GREAT_PVTFGO.exe=2751E8DF3C8E0E80A74883AD847748A41EA5F8D80ADBECA907F8F873C7194B43`，`LibGREAT.dll=714620EA84BE7BCD039902CF69AC4AF326B25B4F9236F484CF73CBD3C1316C66`，`GREAT_GINSFGO.exe=A9184116BA9A2A7ADF2A79CCE4D3598A90E8223C112859DA89FCC3C4D67C2F71`。
- 结论：五频 RAW 的 UPD 与 OSB 固定流程在最终 NL/WL 完整协方差实现上均通过双站两小时验收；Win64 参数地址不再截断，Bias-SINEX 入口正确，OSB 缺失信号按设计局部跳过。下一阶段以这些结果作为 `NONE` 基线，实现参数反馈和可边缘化固定约束反馈。

### 阶段 4：固定结果反馈到因子图（NONE / PARAMETER / CONSTRAINT 三模式）

- 输入提交：`df67f88 fix(fgo): load Bias-SINEX safely on 64-bit` 之后的待提交工作集。
- 目标：把上一阶段"只写 FLT 条件解、图仍全 Float"的缺口补上，让整周固定结果真正反馈进 Ceres 因子图。新增可配置反馈模式 `<fgo><ambiguity_feedback_mode>`（`NONE` / `PARAMETER` / `CONSTRAINT`），三种模式仅作用于固定阶段，浮点估计完全共用。
- 配置解析：`t_gsetfgo::ambiguity_feedback_mode()` 读取 `<fgo><ambiguity_feedback_mode>` 并做大小写归一；未知值抛出 `std::invalid_argument`。`GREAT_PVTFGO` 与 `GREAT_GINSFGO` 两个入口在启动时立即校验该配置，非法值快速失败（打印配置错误并返回 1），避免"退出码 0、固定模式未生效"的静默退化。
- FLT 侧导出（`gambfix/t_gambiguity`）：`t_gambiguity::_addFixConstraint()` 把每条被接受的 NL 绝对双差方程写入 `FixedAmbiguityConstraint`（两端参数索引、卫星、类型、系数 `Ba/Bb`、目标整数、信息量），放入 `_pending_fixed_constraints`；`processBatch()` 在条件更新前调用 `_validateFixedConstraints()` 复核全部导出方程，任一出错（非有限、系数为零、参数越界或残差超容差 `max(1e-6, 10/sqrt(info))`）则回退到 `fltpreAMB` 并拒绝本次候选。只有通过校验后 `_pending` 才转为 `_fixed_constraints`，并 commit 固定历史（`_DD_previous`、连续固定计数），保证被拒绝的候选不会污染后续参考星选取和连续计数。另新增 `_amb_fixed=false` 时代理恢复与约束向量清空。
- 图侧消费（`t_gpfgo`）：新增 `RawFixedConstraint` 与 `fixed_ambiguity_factor.h` 中 Ceres `FixedAmbiguityFactor`（绝对方程 `ca*a + cb*b = target`，权重为 `sqrt_information`）。
  - `NONE`：维持旧行为，仅写条件 FLT，图不反馈。
  - `PARAMETER`：`_apply_RAW_parameter_feedback()` 调用 `_write_RAW_fixed_solution()` 把固定解直接写回图参数并用 `_double_to_vector()` 落地，设 `_graph_ambiguity_fixed=true`。
  - `CONSTRAINT`：`_apply_RAW_constraint_feedback()` 把固定方程作为 `FixedAmbiguityFactor` 残差块逐条加入 Ceres 图，回写参数后重解；失败时用回滚 lambda（删除已安装残差、恢复被替换约束、恢复参数快照、恢复浮点 `_all_para_win` 与旧的固定标志）整段回退到上一有效图。`_rebuild_RAW_posterior_transactional()` 以事务方式重建后验，失败即整体还原。
- 输出：`.flt` 的 `AmbStatus` 现在由 `_graph_ambiguity_fixed ? "Fixed" : "Float"` 决定，即图真正进入固定状态才标 Fixed；被拒绝/回退的历元保持 Float。
- 三种模式浮点结果一致性验证（逐历元 X/Y/Z 全相等，`max_dxyz=0`）：UPD/OSB × GODN/HARB × 2h 与 20m 的 8 组浮点运行，`none-para`、`none-cons`、`para-cons` 逐历元差异均为 0，证明三模式共用同一浮点滤波、差异只出现在固定解上。
- 外部真值：改用 `sample_data/refsnx/IGS0OPSSNX_20233050000_01D_01D_CRD.SNX` 的 `SOLUTION/ESTIMATE` 块坐标（GODN `1130760.6931,-4831298.6759,3994155.1990`、HARB `5084657.6078,2670325.4787,-2768480.8416` m）作为 GODN/HARB 参考，替代此前的 `APX`/共识均值（旧 APX 在 HARB X 方向偏差约 2.4 m）。验证脚本以 `平面收敛门限 10 cm 且高程 20 cm（连续 5 分钟窗口内全部满足）`统计收敛时间、首次固定、首次正确固定与精度（GODN 用 SINEX 估计值、HARB 亦用 SINEX，不再需要共识）。
- 验证结果要点（DF/FF × 双频/浮点 × FLT/FGO-none/para/cons）：
  - 三模式固定率在多数组合无差别；双频 RAW 固定解精度 FLT 与 FGO 同为 cm 级（水平 0.3--0.8 cm，高程 1.3--2.5 cm），已接近 SINEX 参考本身精度。
  - FLT 收敛最快且首固定即正确；GO 首固定极快（0--3.5 min）但常"先错后对"——首正确固定显著晚于首固定（如 `2023305` GODN FF 固定 FGO 首正确 6.5 min、DF-IF 首正确 13--16 min），仅报首固定时间会系统性高估 GO 的收敛速度。
  - FGO 五频固定率低于 FLT（66.7--84.6%），DF-IF 组合固定率最低（56.7--82.5%）且收敛最慢（4--16 min）。
  - 浮点上 FGO 收敛普遍慢于 FLT（约 1.5--3 倍，例 `OSB` HARB FF 浮点 FGO 全程未达 20 cm 高程门限，vRMS 37.7 cm，标"未收敛"）。

### 阶段 4 现存问题（open items）

- 首正确固定 vs 首固定："先固定但定错、后纠正"在 FGO 上普遍出现，尤其 FF 与 DF-IF；当前 `.flt` 只提供 `AmbStatus`，无法从输出区分"错误固定"与"正确固定"，验证脚本需以 SINEX 真值为准后验标注。建议后续输出固定 ratio 与参与固定弧段信息，便于离线甄别 false fixed。
- 三模式的图中真实固定判定差异：`CONSTRAINT` 在候选顶点不足（空区间、窗口滑动）时语义复杂（依赖 `_raw_prior_contains_fixed_information` 与历史约束图），`PARAMETER` 直接覆盖参数，`NONE` 不回馈；三者对"固定后回落浮点"重现速度（收敛保持）的行为不一致，需要长段一致性回归锁定。
- 固定率对固定策略的敏感性：`2023305` 的 HARB DF 固定 `CONSTRAINT` 固定率仅 32.8%、RMS 15--19 cm，明显劣于 `PARAMETER`/`NONE`；这是当前数据里最突出的策略差异，指向约束模式在该站(UPD 模式、低可见性弧段)的数值稳定性或参数化问题，未归入本次验收。
- FGO 进程墙钟/内存：单进程工作集约 1 GB，两小时双站并行总墙钟约 1.9 k s；对项目级多站批处理仍是瓶颈，未优化。
- 现有待办：`t_gambiguity` 的 `_validateFixedConstraints` 对 GLONASS FDMA 非 NL 约束显式拒绝；RAW 后期初态与后验矩阵列序、早期收敛前的错误固定传播（前两历元误差约 0.53--0.59 m）仍需在更长测段和多站上回归。
