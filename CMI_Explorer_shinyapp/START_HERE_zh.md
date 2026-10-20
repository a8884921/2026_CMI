# CMI Explorer 2.0.0 —— 先读这里

网站只有三个页面：**Analysis**（单细胞 CMI、空间 CMI 分级、bulk CMI-riskscore）、**Help**、**About**。
程序包只含计算所需的模型参数，用户上传数据后直接按这些参数算出每个细胞、spot 或样本的 CMI 与分级。

1. **安装与本地试运行**（R ≥ 4.3）：
   ```sh
   Rscript tools/install_packages.R
   Rscript run_local.R        # 本地版；在任何目录下运行都可以（脚本会自己找到本文件夹）
   ```
   macOS 也可以在 Finder 里双击 **`启动CMI_Explorer.command`**，它会自动切到本文件夹再启动并打开浏览器。
   只想跑线上同款的限制模式：`Rscript -e 'shiny::runApp(".")'`（需在本文件夹内执行）。
2. **部署到 shinyapps.io**：
   ```sh
   Rscript tools/deploy_shinyapps.R <账号名> cmi-explorer
   ```
   部署后在 Dashboard → Settings → Advanced 中把 **Max Worker Processes 设为 1**。
3. **修改后检查**：`Rscript tests/test_equivalence.R .` 与 `Rscript tools/check_integrity.R .`。

| 文件 | 内容 |
|---|---|
| `docs/DEPLOYMENT_zh.md` | 部署指南：内存预算、推荐设置、上传限制、常见问题 |
| `README.md`、`CHANGELOG.md` | 英文说明与版本记录 |
