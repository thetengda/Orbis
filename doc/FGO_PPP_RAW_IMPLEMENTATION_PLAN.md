# FGO PPP 非组合模式实现计划

## 1. 目标与范围

在现有 `GREAT_PVTFGO` 的 PPP 因子图处理中增加非组合观测模式，使配置项 `obs_combination` 选择 `RAW_ALL` 时，不再构造无电离层组合（IF）观测，而是按原始频点分别建立伪距和载波相位因子。

首阶段以双频、浮点 PPP 为主，显式估计斜向电离层延迟 `SION`，并保留现有坐标、钟差、对流层、系统间偏差和模糊度状态。随后扩展多系统、滑动窗口、周跳重置、异常值处理、模糊度固定及 OSB 支持。

现有 IF 模式保持原有处理路径，非组合模式通过独立的数据结构、状态容器和因子实现，避免改变已有 IF 方程的行为。

## 2. 非组合观测模型

以参考频率 `f1` 上的斜向电离层延迟为 `SION`，频率 `fi` 的电离层比例系数为：

```text
gamma_i = (f1 / fi)^2
```

伪距观测方程中的主要参数项为：

```text
P_i = rho + c(dt_r - dt_s) + T + gamma_i * SION + ISB + code_bias
```

载波相位观测方程中的主要参数项为：

```text
L_i = rho + c(dt_r - dt_s) + T - gamma_i * SION + AMB_i + ISB + phase_bias
```

其中：

- `CRD`：接收机坐标；
- `CLK`：接收机钟差；
- `TRP`：对流层湿延迟参数；
- `SION`：卫星对应的参考频率斜向电离层延迟；
- `AMB_L1` 至 `AMB_L5`：各频率独立载波模糊度；
- `ISB`：非参考系统相对于参考系统的系统间偏差；
- `code_bias`、`phase_bias`：由当前偏差模型和 OSB 预处理提供。

优先复用 `t_gprecisebiasFGO::cmb_equ()` 生成线性化方程。现有精密模型已经包含 `SION` 的伪距正号、相位负号及频率缩放关系，因此新增因子应从该模型提取设计矩阵、残差和权阵，避免在因子内部重复实现改正模型。

首阶段不引入 `VION`。当前精密模型对 `VION` 尚未提供可直接使用的实现，非组合 PPP 统一使用按卫星维护的 `SION`。

## 3. 原始观测消息结构

在 `src/LibGREAT/gproc/gpvtfgo.h` 中新增原始观测消息结构，例如 `RAWEquMsg`，用于保存单颗卫星、单个频率、单种观测量对应的因子输入。

建议字段包括：

```cpp
t_gtime epoch;
t_gsatdata satdata;
GOBSTYPE obs_type;
GOBS obs;
FREQ_SEQ freq;
GOBSBAND band;
string site;
string sat_id;
string amb_id;
string ion_id;
```

同时增加：

- 当前历元原始消息容器 `_RAW_msg`；
- 窗口内原始消息容器 `_vRAW_msg`；
- 用于稳定标识模糊度弧段和电离层参数的键值生成函数。

消息中必须保存最终选用的精确信号 `GOBS`，不能只保存频率序号。这样才能正确处理同一频率上的码型切换、偏差查询和 OSB 完整性判断。

## 4. 观测组合入口改造

当前 PPP FGO 路径无条件调用 `_combine_IF()`。需要在 PPP 历元处理入口根据 `_observ` 或 `obs_combination` 分流：

```cpp
if (_observ == OBS_COMBIN::IONO_FREE) {
    _combine_IF(...);
} else if (_observ == OBS_COMBIN::RAW_ALL) {
    _combine_RAW(...);
}
```

新增 `_combine_RAW()`，主要职责为：

1. 遍历当前卫星可用频率和观测类型；
2. 依据现有选频规则确定每个频率实际采用的码和相位信号；
3. 对观测执行天线、潮汐、相位缠绕、钟差、轨道和偏差等已有预处理；
4. 为每条有效伪距和载波观测生成一条 `RAWEquMsg`；
5. 建立卫星级 `SION` 标识；
6. 为载波观测建立频率级、弧段级模糊度标识；
7. 在 OSB 模式下，仅保留偏差完整且满足当前策略的观测。

`RAW_ALL` 初期按实际存在的频率生成观测，不强制所有卫星具有相同频率数。`RAW_MIX` 可在 `RAW_ALL` 稳定后复用相同框架，通过额外的选频策略实现。

## 5. 新增非组合因子

在 `src/LibGREAT/gfactor/` 中新增原始观测因子：

- `pseudorange_RAW_factor.h/.cpp`；
- `carrierphase_RAW_factor.h/.cpp`；
- 如现有代码按单系统和多系统拆分，则同步增加对应的多系统版本。

GPS 或参考系统伪距因子的参数块：

```text
CRD, CLK, TRP, SION
```

GPS 或参考系统相位因子的参数块：

```text
CRD, CLK, TRP, SION, AMB_Li
```

非参考系统在上述参数块基础上增加相应 `ISB`。

每个因子需要：

- 从 `RAWEquMsg` 获得卫星、历元、信号、频率和站点信息；
- 将 Ceres 参数块映射为 `t_gallpar` 中的对应参数；
- 调用 `t_gprecisebiasFGO::cmb_equ()` 计算线性化观测方程；
- 从设计矩阵提取各状态的雅可比；
- 使用模型返回的残差及权信息构造加权残差；
- 对不存在的系统间偏差参数采用与现有 IF 因子一致的系统映射规则。

因子中不要自行再次施加 OSB。OSB 应由观测预处理和精密偏差模型统一处理，防止重复改正。

## 6. FGO 状态容器扩展

在 `src/LibGREAT/gproc/gfgo_estimator.h` 及相关实现中增加：

- 非组合模糊度状态容器 `_para_AMB_RAW`；
- 斜向电离层状态容器 `_para_SION`；
- 对应的历元、卫星、频率和弧段索引；
- 必要的先验值、过程噪声和生命周期信息。

Ceres 注册参数块后要求参数地址稳定，因此状态容器不能在优化期间因扩容或重新组织而导致地址变化。可使用节点对象、稳定地址的关联容器，或在建图前完成全部内存分配。

推荐键定义：

```text
SION key = epoch/node + satellite
AMB key  = satellite + frequency + arc_number
```

`SION` 是随历元变化的窗口状态；模糊度在连续弧段内跨历元共享，发生周跳、信号切换或长时间失锁时创建新弧段。

## 7. 原始模糊度管理器

参照现有 `gambIF_manager` 新增 `gambRAW_manager`，负责非组合模糊度状态的创建、复用和释放。

主要功能包括：

- 以“卫星、频率、弧段号”为唯一键管理模糊度；
- 根据周跳标志为受影响频率开启新弧段；
- 在码型或相位信号切换时判断是否能够延续原弧段；
- 处理卫星进入、退出和重新进入窗口；
- 记录模糊度首次出现及最后使用的节点；
- 在窗口滑动和回滚时同步更新引用关系；
- 提供浮点模糊度和固定模糊度之间的索引映射；
- 在 OSB 模式下记录各频率相位 OSB 是否完整。

模糊度初值可参考滤波器 `t_gcombALL` 和 `gpvtflt.cpp` 的 RAW 初始化逻辑，使用同频相位与伪距差值初始化，并结合当前参数单位完成转换。不能沿用 IF 模糊度的初始化公式。

## 8. SION 状态建模

每颗卫星每个窗口节点维护一个 `SION` 状态。首个有效历元可使用当前模型或伪距组合结果初始化，后续节点通过随机游走连接：

```text
SION_k - SION_(k-1) = 0
```

过程噪声应按实际历元间隔 `dt` 缩放，而不是假定固定采样间隔。卫星中断后重新出现时，根据中断时长决定继续旧过程链或重新初始化。

需要增加独立的 `SION` 过程因子，或复用现有随机游走因子并为其提供卫星级状态索引。窗口内没有任何有效观测引用的 `SION` 参数不得继续注册到优化问题中。

## 9. PPP 优化流程改造

修改 `t_gpvtfgo::_optimization_PPP()`，根据观测组合模式分别建立 IF 或 RAW 图。

RAW 路径需要完成：

1. 注册窗口节点的坐标、钟差、对流层和系统间偏差参数；
2. 注册每颗卫星各节点的 `SION` 参数；
3. 注册每个有效频率弧段的模糊度参数；
4. 为每条码观测添加 `PseudorangeRAWFactor`；
5. 为每条相位观测添加 `CarrierphaseRAWFactor`；
6. 添加相邻节点间的运动学、钟差、对流层及 `SION` 过程约束；
7. 仅在新弧段建立时施加模糊度初始先验；
8. 接入上一窗口生成的边缘化先验；
9. 保持现有鲁棒核、参数固定策略和多系统钟差策略。

IF 和 RAW 两条路径应共享公共节点状态的注册代码，但分别管理观测因子、模糊度和电离层状态。

## 10. 状态与向量转换

当前多个函数只处理 `_para_AMB_IF`。以下流程需要增加 RAW 分支：

- `_double_to_vector()`：将 `_para_AMB_RAW` 和 `_para_SION` 写入统一参数向量；
- `_vector_to_double()`：把优化结果同步回 RAW 状态容器；
- 初始值设置与恢复：为新节点和新弧段填充合理初值；
- 参数输出：按卫星、频率和弧段输出 RAW 模糊度，并按卫星输出 `SION`。

参数匹配必须使用稳定标识，不能依赖容器遍历顺序，否则窗口滑动后容易将状态写回错误的卫星或频率。

## 11. 周跳、信号切换与卫星状态变化

将当前针对 IF 模糊度的周跳处理扩展为逐频率处理：

- 某频率发生周跳时，仅重建该频率模糊度弧段；
- 多频组合检测无法确定具体频率时，可重建该卫星全部相位弧段；
- 相位观测码型切换且相位偏差不可连续时，开启新弧段；
- 卫星短时失锁是否保留弧段，由现有数据间断阈值控制；
- 卫星重新进入后重新建立必要的 `SION` 状态和观测消息关联。

原始观测消息、模糊度管理器和窗口节点对周跳结果必须使用同一弧段编号。

## 12. 滑动窗口、回滚与边缘化

以下现有 IF 专用流程需要增加 RAW 状态处理：

- `_remove_outlier_sat()`；
- `_rollback_current_ppp_node()`；
- `_slide_window()`；
- `_marginalization_PPP()`。

窗口滑动时：

- 删除移出节点对应的 RAW 观测消息；
- 删除该节点的 `SION` 状态；
- 对仍被后续节点引用的模糊度弧段保留状态；
- 对不再被任何观测引用的模糊度状态进行清理；
- 将被边缘化节点中与保留变量相关的信息写入边缘化先验。

回滚当前节点时，需要同时恢复：

- RAW 观测消息；
- 新创建的 `SION` 状态；
- 当前节点新增的模糊度弧段；
- 弧段引用计数和周跳状态；
- 优化前的公共参数值。

边缘化参数块列表必须包含 RAW 模糊度和 `SION`。对于跨窗口共享的模糊度，只能保留或重参数化，不能随最老节点一并删除。

## 13. 异常值与验后处理

修改以下 PPP 流程，使其能够识别 RAW 观测因子：

- `_posteriori_test_PPP()`；
- `_gobs_outlier_detection()`；
- 卫星级异常剔除逻辑；
- 因子残差到原始观测的反向索引。

异常值索引至少包含历元、卫星、观测类型、频率和信号。码异常只删除对应码因子；相位异常需要根据现有策略决定仅删除当前观测，还是同时结束该频率模糊度弧段。

不能继续用“每颗卫星固定两条 IF 观测”的位置关系推断残差归属，因为 RAW 模式下每颗卫星的因子数会随可用频率变化。

## 14. 多系统支持

在 GPS 双频 RAW 路径完成后扩展到 BDS、Galileo、GLONASS 和 QZSS 等系统。

需要处理：

- 各系统频率序号到实际载波频率和波段的映射；
- 各系统参考频率及 `gamma_i` 计算；
- 非参考系统 `ISB` 参数块；
- GLONASS 频分多址下按卫星频率计算波长和电离层比例；
- 各系统信号优先级与码型切换；
- 不同系统可用频率数量不一致的情况。

多系统因子应复用现有 IF 多系统因子的钟差和 `ISB` 组织方式，避免另建一套系统时钟定义。

## 15. 非组合模糊度固定

扩展 `_pre_amb_resolution()` 和 `_amb_resolution()`，使模糊度固定模块能够使用 `AMB_L1` 至 `AMB_L5`，而不是只识别 `AMB_IF`。

主要修改包括：

- 从 FGO 状态容器提取逐频率浮点模糊度及协方差；
- 建立卫星、频率、弧段与 `t_gallpar` 参数索引的映射；
- 生成宽巷、窄巷或原始整数模糊度组合；
- 根据固定结果向图中添加约束或更新固定解；
- 仅对连续弧段且偏差条件完整的模糊度参与固定；
- 固定失败时保持浮点状态，不破坏窗口内原始参数。

该部分应尽量复用现有滤波 PPP RAW 模糊度参数类型和 `gambiguity` 的组合关系，FGO 层主要负责状态提取、协方差映射与固定约束回写。

## 16. OSB 模式衔接

现有 OSB 读取和观测改正位于公共观测及偏差模型层，可直接服务于 RAW 观测。FGO RAW 路径还需补齐以下约束：

- `RAWEquMsg` 保存精确 `GOBS`，以便判断对应码和相位 OSB；
- `_combine_RAW()` 在建因子前确认所需 OSB 是否存在；
- 码 OSB 以米为单位进入码观测改正；
- 相位 OSB 按对应频率波长转换后进入相位观测改正；
- 模糊度固定前检查参与组合的各频率相位 OSB 完整性；
- OSB 已包含的偏差项不得在因子模型或固定模块中重复施加；
- 保持 `<upd_mode>OSB</upd_mode>` 和传统 UPD 模式互斥。

RAW 浮点 PPP 可以允许部分频率缺少相位 OSB，并按现有策略排除对应相位或模糊度；RAW 固定解必须要求参与固定的信号偏差完整。

## 17. 配置与模式限制

配置层继续使用现有 `obs_combination`：

```xml
<obs_combination>RAW_ALL</obs_combination>
```

需要在 `GREAT_PVTFGO` 初始化阶段显式检查：

- PPP 模式下 `IONO_FREE` 进入现有 IF 路径；
- PPP 模式下 `RAW_ALL` 进入新增 RAW 路径；
- 尚未实现的组合模式给出明确错误或警告；
- RAW 模式必须启用 `SION` 状态；
- RAW 与 `VION` 的不兼容配置应在建图前被拒绝；
- OSB 模式继续由 `<upd_mode>OSB</upd_mode>` 指定。

配置帮助文本应明确说明 FGO PPP 对 `IONO_FREE`、`RAW_ALL` 和 `RAW_MIX` 的实际支持状态。

## 18. 主要文件修改清单

预计涉及以下文件或目录：

```text
src/LibGREAT/gproc/gpvtfgo.h
src/LibGREAT/gproc/gpvtfgo.cpp
src/LibGREAT/gproc/gfgo_estimator.h
src/LibGREAT/gproc/gfgo_estimator.cpp
src/LibGREAT/gambfix/gambRAW_manager.h
src/LibGREAT/gambfix/gambRAW_manager.cpp
src/LibGREAT/gfactor/pseudorange_RAW_factor.h
src/LibGREAT/gfactor/pseudorange_RAW_factor.cpp
src/LibGREAT/gfactor/carrierphase_RAW_factor.h
src/LibGREAT/gfactor/carrierphase_RAW_factor.cpp
src/LibGREAT/gfactor/（多系统 RAW 因子文件）
src/LibGREAT/gmodels/gprecisebias.cpp
src/LibGREAT/gmodels/gprecisebiasFGO.cpp
src/LibGREAT/gset/（观测组合配置和帮助文本）
src/LibGREAT/CMakeLists.txt 或对应源文件清单
```

`gprecisebias` 相关文件优先只做接口适配；若现有 `cmb_equ()` 已能完整返回 RAW 方程和 `SION` 系数，则不修改核心观测模型。

## 19. 建议实施顺序

### 阶段一：GPS 双频浮点 RAW PPP

- 建立 `RAWEquMsg` 和 `_combine_RAW()`；
- 新增码、相位 RAW 因子；
- 增加 `AMB_L1/AMB_L2` 和 `SION` 状态；
- 在 `_optimization_PPP()` 中建立 GPS 双频 RAW 图；
- 保持 IF 路径不变。

### 阶段二：窗口状态完整管理

- 增加 `gambRAW_manager`；
- 完成周跳、弧段切换和卫星进出；
- 补齐向量转换、回滚和窗口滑动；
- 将 `SION` 与 RAW 模糊度纳入边缘化。

### 阶段三：多系统和多频

- 增加各系统频率映射；
- 接入 `ISB` 参数；
- 扩展 `AMB_L3` 至 `AMB_L5`；
- 处理 GLONASS 实际频率和波长；
- 补齐多系统 RAW 因子。

### 阶段四：异常值与模糊度固定

- 将残差索引改为历元、卫星、频率和信号联合索引；
- 扩展验后检验和异常值剔除；
- 接入逐频率浮点模糊度协方差；
- 实现 RAW 模糊度固定及固定约束回写。

### 阶段五：OSB 和扩展组合模式

- 完成 RAW 观测的 OSB 完整性检查；
- 将 OSB 条件接入非组合模糊度固定；
- 在同一 RAW 框架上扩展 `RAW_MIX`；
- 更新配置帮助和运行期模式提示。

## 20. 实现原则

- IF 和 RAW 使用明确的模式分支，公共状态处理尽量复用；
- 观测偏差只在统一模型层施加一次；
- 所有状态通过稳定键关联，不依赖容器位置；
- 原始频率模糊度按弧段管理，不与 IF 模糊度共用容器；
- `SION` 按卫星、按节点维护，并使用实际历元间隔建模；
- 新增功能按“浮点、窗口管理、多系统、固定、OSB”的顺序逐步接入；
- 保持现有 IF PPP、RTK FGO 和滤波 PPP 的代码路径不变。
