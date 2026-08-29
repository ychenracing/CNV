# ChatGPT Project Brief

> 本文件只保存长期稳定、仓库级的信息。当前任务、临时分支、SHA、测试状态和执行进度应保存在当前 Pull Request 正文中。

## 1. Project

- 项目名称：CNV
- GitHub 仓库：`ychenracing/CNV`
- 默认分支：`master`
- 系统定位：用于实验室 CNV 区域统计、重叠分析、覆盖处理和精确率/召回率计算的 Java 工具集合。
- 项目最终目标：更完整的目标与验收口径未在仓库文档中明确。

## 2. Purpose and Non-Goals

仓库包含针对不同 CNV 工具结果、模拟区域、WES 结果和已知 CNV 区域的读取、规范化、重叠与统计程序。

长期非目标未在仓库文档中明确。仓库没有服务、用户界面、通用工作流平台、发布包或生产运行说明，不应自行假定这些职责。

## 3. Architecture and Module Boundaries

- `src/utils/`：共享 `Region` 和 `Pair` 数据结构。
- `src/cnvstatistic/`：CNV 工具结果与模拟区域的统计和 ROC 准备。
- `src/wes/`：WES CNV 结果、已知区域及不同工具间的重叠和精确率/召回率分析。
- `src/excavator/`：EXCAVATOR 目标覆盖与注释转换工具。
- `src/placenta/`、`src/seqcnv/`：特定数据集或工具的评估辅助程序。

各 Java `main` 类是对应批处理入口；`utils.Region` 是区域重叠语义的共享 Owner。输入文件格式与环境路径由具体程序代码决定，仓库未提供统一配置层。不得在治理文档中发明第二套区域 schema 或统一运行配置。

## 4. Non-Negotiable Constraints

- CNV 区域的染色体、起止坐标和重叠计算语义必须以 `src/utils/Region.java` 及调用代码为准。
- 不得把某个 `main` 方法中的环境专用文件路径描述为可移植默认配置。
- 输入文件格式、表头处理和工具特定转换不得在无样本证据时猜测。
- 仓库未定义医学、临床或诊断用途；不得把研究统计工具描述为临床结论系统。
- 研究数据、样本标识或本地绝对路径不得复制到治理文件或 PR 模板中。

## 5. Authoritative Sources

- 项目简介：`README.md`
- 工程约定：`AGENTS.md`
- 区域模型：`src/utils/Region.java`
- 统计入口：`src/cnvstatistic/Statistic.java`
- WES 重叠分析：`src/wes/OverlapAnalysis.java`
- 其他工具入口：`src/excavator/`、`src/placenta/`、`src/seqcnv/`
- 构建、依赖、测试、版本和发布权威来源：未在仓库中定义

## 6. Standard Commands

- 安装：不适用；仓库未定义安装流程。
- 构建：未在仓库中定义；没有 Maven、Gradle 或 Make 配置。
- 单元测试、集成测试、lint、类型检查和格式检查：未在仓库中定义。
- 本地运行：各 Java 类含独立 `main` 方法，但统一编译 classpath、JDK 版本、参数和输入准备未在仓库中定义。
- 关键验收命令：未在仓库中定义。

## 7. Important Paths

- `src/utils/Region.java`：共享 CNV 区域表示与重叠逻辑。
- `src/cnvstatistic/`：统计和 ROC 数据准备。
- `src/wes/`：WES 重叠、精确率和召回率分析。
- `src/excavator/`：目标覆盖和注释转换。
- `src/placenta/`：胎盘数据相关评估。
- `src/seqcnv/`：SeqCNV 辅助分析。
- `README.md`：简要项目说明。
- `AGENTS.md`：渐进式验证约定。

## 8. CI and Acceptance Entry Points

- 仓库没有 `.github/workflows/`，未定义 GitHub Actions 门。
- 本地验证应遵循 `AGENTS.md` 的影响范围驱动原则。
- 项目没有已定义的统一 Definition of Done；行为改动至少需要用匹配的输入格式验证受影响 `main` 类和共享 `Region` 逻辑。
- 纯文档治理只需验证 Markdown、路径、引用和 diff 范围。

## 9. Prohibited Actions

- 不得把研究统计结果描述为临床诊断或医学建议。
- 不得把环境专用绝对路径、样本数据或隐私信息写入治理文件。
- 不得擅自改写 Git 历史或 force push。
- 不得丢弃未知或未提交工作，也不得覆盖无关改动。
- 不得把计划执行写成已验证完成。
- 不得根据旧聊天猜测当前分支、SHA、PR 或 CI 状态。

## 10. Context Loading Protocol

1. 新开发任务可以直接使用自然语言提出，不要求预先填写固定 Prompt。
2. 开始任务时先读取本文件。
3. 搜索与任务相关的开放 PR、分支和 Issue。
4. 如果存在匹配工作，从现有现场原地继续。
5. 当前动态任务状态默认维护在 Pull Request 正文。
6. 不强制普通单 PR 任务创建 Issue。
7. 优先读取目标代码、直接调用者、相关测试和直接相关配置。
8. 只有证据不足、状态冲突或影响范围扩大时才扩大读取。
9. 不默认加载完整仓库、完整聊天、完整日志或全部 GitHub Actions 历史。
10. 长对话交接使用 `conversation-continuity-guard`，但 GitHub 当前现场仍是状态权威来源。

## 11. References

- `README.md`
- `AGENTS.md`
- `src/utils/Region.java`
- `src/cnvstatistic/Statistic.java`
- `src/wes/OverlapAnalysis.java`
