import { useRef, useState } from "react";

const HOVER_DELAY_MS = 1000;

type ParamTooltipProps = { text: string };

export function ParamTooltip({ text }: ParamTooltipProps) {
  const [show, setShow] = useState(false);
  const timeoutRef = useRef<ReturnType<typeof setTimeout> | null>(null);

  const onEnter = () => {
    timeoutRef.current = window.setTimeout(() => setShow(true), HOVER_DELAY_MS);
  };

  const onLeave = () => {
    if (timeoutRef.current) {
      window.clearTimeout(timeoutRef.current);
      timeoutRef.current = null;
    }
    setShow(false);
  };

  return (
    <span className="param-tooltip-wrap">
      <span
        className="param-tooltip-trigger"
        onMouseEnter={onEnter}
        onMouseLeave={onLeave}
        onFocus={onEnter}
        onBlur={onLeave}
        role="img"
        aria-label="参数说明"
      >
        ?
      </span>
      {show && (
        <span
          className="param-tooltip-bubble"
          role="tooltip"
          onMouseEnter={onEnter}
          onMouseLeave={onLeave}
        >
          {text}
        </span>
      )}
    </span>
  );
}

/** 各参数名称对应的说明文案（单细胞/LabFlow 语境，面向非专业用户） */
export const PARAM_TOOLTIPS: Record<string, string> = {
  "数据名称":
    "给这个数据集起一个名字，方便后面区分多个数据。例如：colon_TC、筋膜_TC。不填也可以。",

  "主数据文件":
    "单细胞测序数据文件。支持 .h5ad（AnnData）、.rds（R/Seurat）、.h5seurat。选你电脑上的文件即可；若是 rds/h5seurat，下一步「校验」会自动转成 h5ad。",

  "SMILES CSV":
    "仅在勾选「使用药物结构模式」时需要。上传一个 CSV，里面包含药物的 SMILES 和剂量等信息，用于药物结构相关的训练。",

  "使用药物结构模式":
    "勾选后会用「药物结构 + 剂量」参与训练，需要同时上传 SMILES CSV。一般做转录组预测时不勾选。",

  "R 执行方式":
    "后端怎么调用 R：Windows 且 R 装在 Conda 里选「cmd_conda」；Linux/mac 且 R 已在系统 PATH 里可选「direct」。",

  "Rscript 命令":
    "在 direct 模式下，后端执行的命令，一般是 Rscript。若 R 在 Conda 里且选了 cmd_conda，这项可忽略。",

  "conda.bat 路径":
    "Windows 下 Conda 的激活脚本完整路径，例如 F:\\software\\Miniconda3\\condabin\\conda.bat。选下拉里自动探测到的，或自己填路径。",

  "Conda R 环境名":
    "Conda 里安装 R/Seurat 的环境名称，例如 r-4.3、r-seurat。选下拉里的环境即可；没有的话先在本机用 conda 建好 R 环境。",

  "分组字段（group_column）":
    "metadata 里用来分组的列名，例如样本、处理条件（sample、condition、Group）。可参考上一步「Seurat 解析」里列出的 metadata 字段。",

  "类群字段（cluster_column）":
    "metadata 里表示细胞类型或聚类结果的列名，例如 celltype、seurat_clusters、cell_type。训练时会按这一列筛选要用的细胞。",

  "随机种子（seed）":
    "抽样和随机过程的种子，固定后结果可复现。默认 42 即可，一般不用改。",

  "待筛选 clusters（逗号分隔）":
    "要从类群字段里保留哪些取值参与训练，用逗号分隔。例如填 T,B,NK 表示只训练 T、B、NK 这几类细胞。可先点「Seurat 解析」里的类群列名查看有哪些取值。",

  gene_size:
    "输入基因数量，要和预处理后的数据一致（通常预处理会给出建议）。单细胞矩阵的基因维数，模型输入层大小。",

  output_dim:
    "模型隐层/输出维度，一般和 gene_size 设成一样。影响表达向量的维度。",

  batch_size:
    "训练时每批喂给模型的样本数。默认 64；显存不够可以改小（如 32）。",

  lr: "学习率，控制参数更新步长。默认 1e-4，一般不用改；调参时可适当减小或增大。",
};
