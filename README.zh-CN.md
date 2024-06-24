
<h1><p align='center'>dms-pipeline: a Pipeline for High-throughput Antibody DMS Data Processing</p></h1>

<p align='center'>
  <a href='./README.md'>English</a>  | 简体中文
</p>

本程序用于高通量深度突变扫描数据的处理和分析。上游数据处理部分遵循[J. Bloom实验室的思路](https://github.com/jbloomlab/SARS-CoV-2-RBD_DMS)，并进行了一些修改以适应高通量数据。

您可以参考我们实验室发表的文章[这里](#Citation)获取更多信息。

值得注意的是，由于已发表文章中呈现的一些数据可能并未通过最新的方法处理，文章中显示的详细结果，包括其对应的Github仓库中的结果，可能与本仓库方法生成的结果不完全相同，但这些结果不应存在结论上的差异。

建议使用此仓库中的最新版本方法。

## 安装
### 准备conda

本方法的大部分脚本是用Python编写的，我们将使用conda来设置环境。为了加速环境解析和安装，推荐使用[Miniforge](https://github.com/conda-forge/miniforge)。

如果使用Miniforge，在以下命令中将`conda`替换为`mamba`。

首先，通过检查当前可用的环境列表，确保conda（或mamba）准备就绪。

```bash
conda env list
```

### 创建环境

我们提供了一个`environment.yml`文件。请将此仓库克隆到您的工作目录并运行

```bash
conda env create -f environment.yml
```

这将构建一个名为`dms-pipeline`的环境，然后激活它

```bash
conda activate dms-pipeline
```
### 安装其他依赖项

我们将使用`pip`安装一些Anaconda未分发的Python依赖。这些包列在`requirements.txt`中。

```bash
pip install -r requirements.txt
```

### 编译C/C++代码（可选）

由于性能原因，这个程序的一小部分是用C/C++编写的。这些代码提供了一些可选的特性，但对于基本的分析不是必需的。如果您需要它们，请进入`scripts`目录并编译它们。目前，这些代码非常简单，因此不应在依赖项上存在复杂问题。

您需要通过`pip`或`conda`安装`pybind11`以将二进制文件导出为Python模块。

```bash
conda install -c conda-forge pybind11

cd scripts
chmod +x compile.sh && ./compile.sh
```

如果不能成功编译，请尝试使用conda安装更新版本的GCC：

```bash
conda install -c conda-forge compilers
```

## 使用
### 配置

请参阅示例配置文件`config.yaml`。您应该将文件路径和其他必要信息根据您的需求进行替换。

**（待完成）**

### 分析PacBio测序数据

使用PacBio测序数据构建barcode-variant查找表。

**（待完成）**

### 分析barcode测序数据

### 下游分析

如果已经为一个抗原或几个密切相关的抗原收集了足够的抗体DMS结果，则可以进行下游分析，以全面了解该抗原的BCR表位分布。
我们将下游分析的核心代码包含在另一个名为[HADI](https://github.com/yunlongcaolab/hadi)的包，HADI代表<ins>**H**</ins>igh-throughput <ins>**A**</ins>ntibody <ins>**D**</ins>MS <ins>**I**</ins>ntegration😄。

因为下游分析并非每次都需要进行，使用HADI进行下游分析的脚本未集成到本仓库的Snakemake流程中。

HADI的运行需要正确编写对应的YAML配置文件并运行`scripts/dms_integration.py` **（待完成）**。

## TODO

- 包括PacBio数据分析脚本。
- 包括Sort-seq和Tite-seq分析脚本。
- 将DMS实验的质量控制（QC）部分从HADI转移到dms-pipeline中。

## 引用

如果这个数据分析流程对您有帮助，请引用以下文章：
- [Cao et al. Nature 2022](https://doi.org/10.1038/s41586-022-04980-y)
- [Cao et al. Nature 2023](https://doi.org/10.1038/s41586-022-05644-7)
- [Yisimayi et al. Nature 2024](https://doi.org/10.1038/s41586-023-06753-7)

您可以在[曹云龙实验室官方网站](https://yunlongcaolab.com)找到更多信息。

---