# CMI Explorer 部署指南（shinyapps.io）

## 1. 架构

程序包只保存计算所需的**模型参数**（树结构、分裂点、叶节点值、预处理常数）和计算五个维度所需的参考文件，
用纯 R 代码完成全部计算（20 个 CRAN 包，无需 Java 或 Python）。
界面进程只加载 shiny/bslib/DT（约 190 MB），每次分析在独立的 R 子进程中运行，结束即释放内存；
同一实例同一时间只运行一个分析，其余按先来后到排队，并显示排队位置。上传的程序包约 30 MB。

## 2. 本地准备（只需一次）

需要 R 4.3 或更高版本（建议 R 4.4）。在 CMI_Explorer 文件夹中运行：

```sh
Rscript tools/install_packages.R
```

**重要：** 请在使用 CRAN 安装包的 R 环境中部署（Windows/macOS 的常规 R 即可）。
shinyapps.io 会按本机包的来源重新安装依赖；Linux 发行版自带的 r-cran-* 包可能缺少 DT 的部分 JS 文件。

## 3. 网站结构

网站只有三个页面：**Analysis**（单细胞 / 空间 / bulk 三个分析，用顶部按钮切换）、**Help**、**About**。
不再包含引用、文献、作者信息和模型说明页面，上线前无需填写任何出版信息。

## 4. 本地试运行

```r
shiny::runApp(".")        # 托管版：与 shinyapps.io 上的限制完全相同
```

在 Analysis 页依次切换 Single-cell / Spatial / Bulk，点击 “Try the example”（立即显示预计算结果），
再点 “Run the example live” 确认实时计算可用。

（可选）运行验证测试：`Rscript tests/test_equivalence.R .`，应显示 “16 of 16 checks passed”
（未安装 spacexr 时 RCTD 一项会显示 skipped，属正常）。

## 5. 部署

1. 登录 shinyapps.io → Account → Tokens → Show → 复制命令，在 R 中执行：
   ```r
   install.packages("rsconnect")
   rsconnect::setAccountInfo(name = "<账号名>", token = "<token>", secret = "<secret>")
   ```
2. 先做一次演练（不需要账号，不上传任何文件）：
   ```sh
   Rscript tools/deploy_shinyapps.R --dry-run
   ```
   它会做下面全部检查，并列出将上传的文件数、大小和服务器要安装的 R 包。
3. 在 CMI_Explorer 文件夹中运行：
   ```sh
   Rscript tools/deploy_shinyapps.R <账号名> cmi-explorer
   ```
   脚本会依次检查：`ref/` 中模型文件的 SHA-256、预计算示例是否为本版本生成、
   `tools/check_integrity.R`（锁定文件未改动、R 源文件为 ASCII、`ref/SHA256SUMS` 为最新、文件中没有令牌）；
   然后只上传网站需要的文件（app.R、R/、engine/、ref/、www/、examples/ 与 .renvignore），
   并把实例大小设为 Large（1 GB）。tools/、tests/、docs/ 不会上传。
4. 首次部署约需 5–10 分钟（服务器安装 R 包）。完成后访问 `https://<账号名>.shinyapps.io/cmi-explorer/`。

也可以在 RStudio 中点 “Publish”：`.rscignore` 已排除 tools/、tests/、docs/ 等（但不会做上述检查）。

## 6. shinyapps.io 推荐设置（Dashboard → Applications → cmi-explorer → Settings）

| 设置 | 推荐值 | 说明 |
|---|---|---|
| Instance Size | Large (1 GB) | 托管限制按 1 GB 设计；付费计划可选 xlarge（2 GB）并放宽限制 |
| **Max Worker Processes** | **1** | 每个实例只保留一个界面进程（约 190 MB）+ 一个分析（≤ 700 MB），总计 < 0.9 GB；设为 2–3 时可能超过 1 GB |
| Max Connections per Worker | 20–50 | 界面本身很轻；分析在排队中依次运行 |
| Connection Timeout / Read Timeout | 900 s / 3600 s | 排队时分析可能需要几分钟 |
| Instance Idle Timeout | 15 分钟或更长 | |
| Instances（付费计划） | 按访问量 | 每个实例独立排队，可并行处理多个分析 |

## 7. 托管版上传限制（可调整）

| 工作台 | 默认上限 | 实测峰值内存 / 时间 |
|---|---|---|
| 单细胞 | 10,000 个细胞，且 ≤ 1,500 万非零计数；上传 ≤ 300 MB | 676 MB / 47 s |
| Bulk | ≤ 2,500 万个表达值（如 20,000 基因 × 1,250 样本） | 601 MB / 12 s |
| 空间 | ≤ 20,000 spots / 切片 | — |
| RCTD 反卷积 | 仅本地版（3,000 spots 约需 2.2 GB）；托管版提供 Mal-1/2/3 截断点与可在本机运行的 R 脚本（www/downloads/rctd_mal_states.R，附测试数据） | — |

如果使用 2 GB 以上实例，可在 CMI_Explorer 文件夹中新建 `.Renviron` 放宽限制，例如：

```
CMI_MAX_CELLS=20000
CMI_MAX_NNZ=30000000
CMI_MAX_BULK_ENTRIES=50000000
```

部署脚本会自动上传存在的 `.Renviron`。（`CMI_ENABLE_RCTD=true` 可在托管版开启 RCTD，仅建议 ≥ 4 GB 实例。）

## 8. 本地版（文件最大 4 GB，不限细胞数与 spot 数，含 RCTD）

```sh
Rscript tools/install_packages.R --rctd
Rscript run_local.R
```

批量处理多个样本（不需要浏览器，每个样本一次）：

```sh
for f in samples/*.zip; do s=$(basename "$f" .zip)
  CMI_PROFILE=local Rscript engine/worker.R --root . --module single_cell --input "$f" --output "results/$s" --state "results/$s/state.json"
done
```

## 9. 修改后需要重新生成的文件

| 修改了什么 | 运行 |
|---|---|
| `engine/`、`R/figures.R` 或示例输入 | `Rscript tools/precompute_examples.R .`（否则部署脚本会拒绝部署） |
| `ref/` 或 `engine/constants.R` | `Rscript tools/export_gene_sets.R .` |
| 任何锁定文件（经作者批准后） | 运行测试，再运行 `Rscript tools/check_integrity.R . --write-locks` 更新 `tools/locked_files.sha256` |

每次修改后都可以运行 `Rscript tools/check_integrity.R .`。它检查以下几项：
- 锁定文件（模型参数、引擎、测试与夹具、下载脚本、工具）和 `engine/constants.R` 的冻结部分是否与本版本一致；
- R 源文件是否只含 ASCII 并能解析；
- `ref/SHA256SUMS` 是否为最新（过期时运行 `Rscript tools/export_gene_sets.R .` 即可更新）；
- 文件中是否误写了令牌或密钥，软件包里的 `.Renviron` 是否只含 `CMI_*` 设置（它会随部署上传；GITHUB_PAT 等请放在主目录的 `~/.Renviron`）。

部署脚本在上传前还会检查服务器要安装的 R 包：出现分析用包（GEOquery、Seurat、spacexr、h2o 等）或非 CRAN 来源时停止。

## 10. 常见问题

- **部署时提示缺包：** 运行 `Rscript tools/install_packages.R` 后重试。
- **页面显示 “ran out of memory”：** 上传超过了实例可承受的规模，或 Max Worker Processes 大于 1；请减少细胞数或使用本地版。
- **分析显示 “Queued”：** 同一实例上有其他用户的分析在运行，结束后会自动开始（页面显示排队位置）。
- **出现 “Error ID”：** 在 shinyapps.io 的 Logs 中搜索该 ID 可看到详细错误信息。
- **更新网站：** 修改后重新运行部署命令即可（同名应用会覆盖更新）。
