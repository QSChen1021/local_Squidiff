# Results API

用于模型与结果资产的查询、下载与删除。

## GET `/api/results`

列出所有结果记录。

Response:
- `200 OK` -> `{ "items": ResultRecord[] }`

说明：
- `items[].summary.assets[]` 会带上 `name` 与可下载 `url` 字段。

## GET `/api/results/{result_id}`

查询单个结果。

Response:
- `200 OK` -> `{ "result": ResultRecord }`
- `404 Not Found` -> `Result not found`

## GET `/api/results/job/{job_id}`

按任务查询结果。

Response:
- `200 OK` -> `{ "result": ResultRecord }`
- `404 Not Found` -> `Result not found for job`

## GET `/api/results/{result_id}/assets/{asset_name}`

下载结果资产文件（图像、summary 等）。

Response:
- `200 OK` -> file stream
- `404 Not Found` -> `Result not found` / `Asset not found`

## DELETE `/api/results/{result_id}`

删除结果记录（可选同时删除文件）。

Query:
- `purge_files` (bool, optional, default `true`)

Response:
- `200 OK` -> `{ "deleted": true }`
- `404 Not Found` -> `Result not found`

---

## GET `/api/results/models/list`

列出模型记录。

Response:
- `200 OK` -> `{ "items": ModelRecord[] }`

## GET `/api/results/models/{model_id}`

查询单个模型记录。

Response:
- `200 OK` -> `{ "model": ModelRecord }`
- `404 Not Found` -> `Model not found`

## GET `/api/results/models/{model_id}/download`

下载模型 checkpoint 文件。

Response:
- `200 OK` -> file stream
- `404 Not Found` -> `Model not found` / `Model checkpoint not found`

## DELETE `/api/results/models/{model_id}`

删除模型记录（可选同时删除 checkpoint 文件）。

Query:
- `purge_files` (bool, optional, default `true`)

Response:
- `200 OK` -> `{ "deleted": true }`
- `404 Not Found` -> `Model not found`
