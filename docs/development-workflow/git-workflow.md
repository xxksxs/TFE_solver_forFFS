# Git 仓库管理约定

仓库根目录为 `bp_filter/`，默认主分支为 `main`。发布前用 `git remote -v` 核对 `origin` 指向本项目的 GitHub 仓库。

## 仓库基础设施

| 文件 | 用途 |
| --- | --- |
| `.gitignore` | 屏蔽 `build*/`、`result/`、IDE 元数据、LaTeX 中间文件、Python 缓存等 |
| `.gitattributes` | 跨平台行尾稳定化（仓库存 LF，工作副本恢复平台原生 EOL）；`*.aedt`、`*.ngmesh`、`*.pdf`、`*.png`、`*.vtu` 标记为 binary |
| [版本日志](../history/CHANGELOG.md) | 用户可见变化的版本日志，每次里程碑发版同步打 tag |

## 提交粒度

- **一次 commit 解决一件事**：例如"`PortMethod::Numerical` 接入"+对应文档+对应 CHANGELOG 条目算一次 commit；若同时要做无关的代码风格清理，分两次 commit。
- **不混合 release artifact**：`result/`、`build*/`、`*.pdf` 中的临时构建输出不入库；只有作为参考的 LaTeX 编译 PDF（如 `docs/optimization/build_affine_system.pdf`）以及 HFSS reference CSV 等"对照基准"明确入库。
- **commit message**：第一行 ≤ 70 字，命令式语气；如果改动会影响 CLI 或数学行为，body 中列一两行影响面。中英文皆可，本工程默认中文。

## 标签 (tags)

- `baseline-v1`：进入 Git 后的首个干净状态。Annotated tag，`git show baseline-v1` 可读到引入说明。
- 未来里程碑（`alps-mor-v1` / `tfe-multimode-v1` 等）建议都用 annotated tag (`git tag -a`) 而非 lightweight，便于回溯说明。

## 分支策略

- `main` 为默认集成分支。并行工作在隔离分支进行。
- 新功能或重构超过一次 commit 的，开 `feature/<short-name>` 分支，merge 时用 `--no-ff` 保留拓扑：

  ```powershell
  git checkout -b feature/preconditioner-ams main
  # ... commits ...
  git checkout main
  git merge --no-ff feature/preconditioner-ams
  git branch -d feature/preconditioner-ams
  ```

- 实验性 / 抛弃概率高的代码用 `wip/<short-name>` 分支，最终 squash merge 或丢弃。

## 还原与排错

- 工作副本误改：`git restore <file>`（HEAD 版本覆盖）。
- 已 stage 误添加：`git restore --staged <file>`。
- 想看上次 commit 改了什么：`git show HEAD --stat`。
- 想看仓库体积：`git count-objects -vH`；超过 50 MB 时检查是否误纳入大件。
- 想看正在被 ignore 但你不确定的路径：`git check-ignore -v <path>`。

## 发布到远端

完成[发布验收](../testing/release-validation.md)后，推送该次发布分支和对应附注标签。

```powershell
git push --atomic -u origin <release-branch> "refs/tags/<release-tag>"
```

推送前检查变更文件是否包含未授权的数据或临时产物。推送后用 `git ls-remote` 核对远端分支及标签指向。
