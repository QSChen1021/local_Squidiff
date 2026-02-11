# Squidiff LabFlow（本仓库说明）

> 这是一个“研究模型 + 内网 Web 工作流”的混合仓库。  
> 你既可以直接用 `train_squidiff.py` / `sample_squidiff.py` 做模型训练推理，也可以用前后端 Web 流程完成 Seurat 数据上传、500x500 预处理、训练与结果查看。

---

## 1. 项目目标与面向用户

### 项目在解决什么问题？

Squidiff 是一个**用扩散模型预测单细胞转录组变化**的工具（见 [Nature Methods 论文](https://doi.org/10.1038/s41592-025-02877-y)）。通俗说：

- **训练前**：你提供单细胞数据（如 Seurat / h5ad）以及“条件”（例如不同时间点、不同药物/剂量），系统学会“在这种条件下，细胞转录组会变成什么样”。
- **训练后**：用训练好的模型对**新样本、新条件**做**预测**，得到“模型认为的转录组”以及 UMAP、热图等结果，用来做**计算机里的虚拟实验（in silico）**，减少重复湿实验、快速筛条件、辅助发文章。

典型用途包括：**细胞分化轨迹**、**基因扰动**、**药物响应预测**等（与论文中的验证场景一致）。  
更细的“模型能做什么”和设计理念见：`docs/模型能做什么与前端设计理念.md`。

### 产品最终给谁用？

**LabFlow 网页端优先面向不会写代码、不搞生信的生命科学研究者**（硕博、博后、PI 等）。目标是：

- 在**不敲命令、不配环境**的前提下，通过上传数据、点选参数、查看结果，完成“数据 → 训练 → 用模型做预测 → 看图/下载”的全流程。
- 界面和文案要让人一眼看懂：自己在做的是“用 Squidiff 预测转录组”，而不是抽象的“跑一个 AI 模型”。

命令行脚本（`train_squidiff.py` / `sample_squidiff.py`）仍保留给会编程或生信的同学做复现与扩展；但**产品形态与文档以“外行友好”为第一目标**。

---

## 1.5 Quick Start（用户视角）

> 目标：从“刚拿到项目”到“在网页里跑通一次训练并看到结果”。  
> 下面所有文档链接均可直接点击打开。

### Step 0：先选一条部署路径

- **本地开发（推荐第一次使用）**：按 `README` 的“5. 开发与部署（三种方式）”启动前后端，最快看到完整流程。
- **Docker 一键部署（内网环境）**：按 [`docs/部署文档.md`](docs/部署文档.md) 执行，适合团队共享或服务器长期运行。

### Step 1：完成部署并启动服务

1. 按你选择的部署方式完成环境准备。  
2. 启动后端（默认 `http://localhost:8000`）。  
3. 启动前端（默认 `http://localhost:5173`）。  
4. 在浏览器打开前端地址，看到 LabFlow 页面即表示部署成功。

如果你需要从 `.rds/.h5seurat` 转成 `.h5ad`，先看：[`docs/seurat转换指南.md`](docs/seurat转换指南.md)。

### Step 2：按页面完成一次“上传 → 预处理 → 训练 → 查看结果”

建议严格按这个顺序操作（零基础用户可直接照做）：

1. **上传数据并校验**：确认文件可读、字段完整。  
2. **Seurat inspect**：查看 metadata 列与 UMAP 预览。  
3. **500x500 预处理**：筛选 cluster，生成训练用 prepared dataset。  
4. **提交训练任务并轮询状态**：等待任务完成。  
5. **查看结果资产**：模型信息、日志、图像与下载项。

详细“每个按钮点哪里、每个参数怎么填”请看：  
[`docs/LabFlow前端用户操作说明.md`](docs/LabFlow前端用户操作说明.md)

### Step 3：需要快速上手或排查问题时看这些文档

- **10 分钟快速走通一遍**：[`docs/实验室10分钟上手.md`](docs/实验室10分钟上手.md)
- **Seurat 接口与数据流程**：[`docs/api/seurat.md`](docs/api/seurat.md)
- **任务接口（train/predict/log/cancel）**：[`docs/api/jobs.md`](docs/api/jobs.md)
- **数据集接口（上传/校验）**：[`docs/api/datasets.md`](docs/api/datasets.md)
- **UAT 验收清单**：[`docs/UAT_Seurat_V2_检查清单.md`](docs/UAT_Seurat_V2_检查清单.md)
- **产品目标与模型能力说明**：[`docs/模型能做什么与前端设计理念.md`](docs/模型能做什么与前端设计理念.md)

---

## 2. 项目功能总览

### 2.1 研究模型能力（根目录脚本）
- 基于扩散模型的单细胞转录组预测。
- 支持基础模式与药物结构模式（`SMILES + dose`）。
- 入口脚本：
  - `train_squidiff.py`
  - `sample_squidiff.py`

### 2.2 LabFlow Web 能力（`backend/` + `frontend/`）
- 数据上传与格式校验（支持 `.h5ad/.rds/.h5seurat`）。
- Seurat 检查（metadata 字段 + UMAP 预览）。
- 训练前预处理（Phase 2）：
  - cluster 过滤
  - 最多 500 cells 分层抽样
  - 最多 500 genes 筛选（Wilcoxon + fallback）
- 训练任务提交与轮询（Phase 3）：
  - 默认优先使用 `prepared_dataset_id`
  - 训练来源可追溯（`source_dataset_id`、`train_dataset_id`、`prepared_dataset_id`）
- 结果资产查看（模型信息、预测图像、日志）。

---

## 3. 当前架构（前后端 + 任务执行）

```text
Frontend (React/Vite)
   |
   | HTTP / JSON
   v
FastAPI backend
   ├─ /api/datasets   (上传/校验/转换)
   ├─ /api/seurat     (inspect + prepare-training)
   ├─ /api/jobs       (train/predict + poll + log)
   └─ /api/results    (模型/结果/资产)
   |
   v
JsonStateStore (backend/state/*.json)
   |
   v
JobQueue worker
   |
   v
SquidiffRunner -> train_squidiff.py / sample_squidiff.py
```

### 3.1 后端核心目录
- `backend/app/api/`：REST API 路由。
- `backend/app/services/`：业务服务层（转换、检查、预处理、任务执行）。
- `backend/app/storage/state_manager.py`：文件型状态存储。
- `backend/state/`：状态 JSON（`datasets/jobs/seurat_prepare_jobs/models/results`）。
- `backend/uploads/`：上传与预处理输出。
- `backend/artifacts/`：训练/预测任务产物与日志。

### 3.2 前端核心目录
- `frontend/src/App.tsx`：单页流程 UI（上传 -> 校验 -> inspect -> prepare -> train -> 结果）。
- `frontend/src/services/api.ts`：API 类型与请求封装。
- `frontend/src/styles/tokens.css`：样式 token 与页面样式。

---

## 4. API 概览

### 4.1 健康检查
- `GET /api/health`

### 4.2 数据集
- `GET /api/datasets`
- `POST /api/datasets/upload`
- `POST /api/datasets/{dataset_id}/validate`

### 4.3 Seurat（V2）
- `POST /api/seurat/inspect`
- `POST /api/seurat/prepare-training`
- `GET /api/seurat/prepare-training/{job_id}`

### 4.4 任务
- `GET /api/jobs`
- `GET /api/jobs/{job_id}`
- `GET /api/jobs/{job_id}/log`
- `POST /api/jobs/train`
- `POST /api/jobs/predict`

### 4.5 结果
- `GET /api/results`
- `GET /api/results/{result_id}`
- `GET /api/results/job/{job_id}`
- `GET /api/results/models/list`
- `GET /api/results/models/{model_id}`
- `GET /api/results/{result_id}/assets/{asset_name}`

详细接口请看：`docs/api/seurat.md`（Seurat 部分），其余接口可参考 `backend/app/api/*.py`。

---

## 5. 开发与部署（三种方式）

环境三选一即可：**uv**、**conda**、或**本机 Python**（venv + pip）。下面命令按 **Windows** 和 **Linux / macOS** 分开写，复制时只复制你当前系统对应的那一行或一段即可。R + SeuratDisk 仅在需要 `.rds/.h5seurat → h5ad` 转换时安装；前端需 Node.js 20+。**前端每个选项、参数怎么填**见 `docs/LabFlow前端用户操作说明.md`。

### 5.1 环境准备（任选其一）

**uv（推荐）**

```bash
# Windows（在 PowerShell 或 CMD 中执行）
uv venv
.venv\Scripts\activate
uv pip install -r requirements.txt -r backend/requirements.txt

# Linux / macOS（在终端中执行）
uv venv
source .venv/bin/activate
uv pip install -r requirements.txt -r backend/requirements.txt
```

**conda**

```bash
# Windows 与 Linux / macOS 相同
conda create -n labflow python=3.11
conda activate labflow
pip install -r requirements.txt -r backend/requirements.txt
```

**本机 Python（venv + pip）**

```bash
# Windows
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt -r backend/requirements.txt

# Linux / macOS
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt -r backend/requirements.txt
```

### 5.2 启动后端

在项目根目录、已激活环境中执行（Windows 与 Linux / macOS 命令相同）：

```bash
python -m uvicorn backend.app.main:app --reload --host 0.0.0.0 --port 8000
```

### 5.3 启动前端

```bash
# Windows（PowerShell 或 CMD）
cd frontend
npm install
npm run dev

# Linux / macOS
cd frontend
npm install
npm run dev
```

（若习惯一行执行：Linux/mac 与 Windows CMD 用 `cd frontend && npm install && npm run dev`；Windows PowerShell 用 `cd frontend; npm install; npm run dev`。）

- 前端：`http://localhost:5173`
- 后端：`http://localhost:8000`

### 5.4 Docker（内网一键）

```bash
cd infra && cp .env.example .env && docker compose up --build
```

环境变量说明见 `infra/.env.example`，部署细节见 `docs/部署文档.md`。

---

## 7. 开发/运维常用命令

### 7.1 后端静态检查
```bash
ruff check backend/app backend/tests
ruff format --check backend/app backend/tests
```

### 7.2 前端检查与构建
```bash
cd frontend
npm run lint
npm run build
```

### 7.3 训练输入形状检查
```bash
python scripts/check_shape.py --data_path path/to/data.h5ad
```

### 7.4 Phase 4 UAT（至少 2 数据集）
```bash
python scripts/uat_phase4_seurat_v2.py \
  --base-url http://localhost:8000 \
  --dataset-id <A> \
  --dataset-id <B> \
  --group-column sample \
  --cluster-column celltype \
  --selected-clusters T,B,NK \
  --seed 42
```

---

## 8. 目录结构（简版）

```text
.
├─ backend/
│  ├─ app/
│  │  ├─ api/
│  │  ├─ services/
│  │  └─ storage/
│  ├─ state/
│  ├─ uploads/
│  ├─ artifacts/
│  └─ scripts/seurat_to_h5ad.R
├─ frontend/
│  └─ src/
├─ infra/
├─ docs/
├─ scripts/
├─ train_squidiff.py
└─ sample_squidiff.py
```

---

## 9. 文档导航

- **前端用户操作说明（推荐先看）**：`docs/LabFlow前端用户操作说明.md` — 按页面步骤说明每个选项、参数如何填写（含 Windows Conda R 配置、校验/预处理/训练各步）
- **模型能做什么与前端设计理念**：`docs/模型能做什么与前端设计理念.md`（论文依据、训练前后能力、前端根本目标与设计原则）
- 部署与环境：`docs/部署文档.md`
- Seurat 转换（含 V2 补充）：`docs/seurat转换指南.md`
- Seurat API：`docs/api/seurat.md`
- 10 分钟上手：`docs/实验室10分钟上手.md`
- UAT 清单：`docs/UAT_Seurat_V2_检查清单.md`
- 设计/需求：`docs/PRD_Seurat交互筛选与500x500训练管线.md`

---

## 10. 当前状态与注意事项

- V2 Phase 1~4 代码和交付文档已落地（inspect / prepare-training / train 默认 prepared dataset / UAT 资产）。
- 状态存储当前是 JSON 文件方案（MVP 取舍），不等同于数据库事务一致性。
- Seurat 转换依赖本机/容器内 R 运行时与 SeuratDisk 可用。
- 若你在本地开发，建议优先 `LABFLOW_DRY_RUN=true` 验证链路，再切真实训练。

---

## 11. 许可证与引用

- License: MIT（见 `LICENSE`）
- 论文引用见仓库根目录历史信息（`README` 旧版与论文条目）。

## 5.5 ��ҳע��/��¼������������֤��

LabFlow ��֧�������˺�ϵͳ��SQLite ���ؿ� + ��ȫ��ϣ���룩��
- ��ҳ��ע��/��¼���ٵ㡰��ʼ��������
- ��¼��ǰ�˻��Զ�Я�� Bearer Token ���� API��
- �û�˵������ڣ���ҳ���û�˵���顱��ť����� `GET /api/auth/user-guide`����

��� API �ĵ���`docs/api/auth.md`
