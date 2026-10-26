# CislunarSuspensionTether
Cislunar Suspension Tether, My Doctral project

这个代码仓是用来储存我的博士论文项目代码的。
博士论文：《月球天梯设计与运行动力学研究》，中山大学，2025 年 5 月。
仓库只放程序、算例参数和推导文件，仿真得到的数值结果和图像结果（约 1.4 GB 的 SVG 等）留在本地，不上传。

This repository stores the code for my doctoral thesis on the Cislunar Suspension Tether (Sun Yat-sen University, 2025). Only source code, case settings and derivation files are kept here; simulation outputs stay on my computer.

## 目录结构

按论文章节组织，每个文件夹都可以单独运行。

| 文件夹 | 论文章节 | 内容 | 主要文件 |
|---|---|---|---|
| `第2章-月球天梯的标称设计/` | 第 2 章 | 基于静力学和材料力学的构型设计：恒定截面、恒定应力（变截面）两种情况，张力上下界，材料对比 | `静力学绘图.py`、`张力上下界绘图.py`、`各类材料图.py`、`3d截面.py`；Maple 推导 `Constant section tethers …mw`、`Constant tensile stress tethers …mw`、`静力学等截面图20221118.mw` |
| `第3章-月球天梯的整体运动特性/` | 第 3 章 | 珠点模型下天梯随月球绕地公转的自由摆动仿真，含末端偏折、突然断裂等工况 | 主程序 `main-月球天梯自由摆动模拟-20230918.py`；工况 `…-末端偏折.py`、`…-突然断裂.py`；后处理 `main-画图和动画.py`、`专门动画.py`；`圆轨道自由摆动.py`；`jacobi.py` |
| `第4章-月球天梯货物运输耦合动力学/` | 第 4 章 | 货物（轿厢）沿天梯爬升的单节点耦合模型，含变速度、天梯被撞断及救援工况 | 主程序 `爬升-主程序.py`；`爬升-变速度-参数化程序.py`、`爬升-突然断裂.py`；后处理 `爬升-画图和动画.py`、`爬升-画图和动画-突然断裂.py` |
| `第5章-引力虹吸自动运行月球天梯/` | 第 5 章 | 引力虹吸月球天梯（绕动滑轮、定滑轮运转的闭环绳索）的多体动力学仿真与参数敏感性分析 | 主程序 `main-虹吸链条-0426.py`（最终版）；`main-虹吸链条-0322-例1.py`（批量运行版）、`main-虹吸链条-60天.py`；后处理 `main-draw.py`、`main-draw-60days.py`、`main-draw-长度图.py`；`批量运行` |
| `第6章-级联式引力虹吸月球天梯/` | 第 6 章 | 级联式虹吸天梯的“集中-连续”一体化（离散-连续复合）动力学模型 | 主程序 `main-离散-连续-复合模型.py`；后处理 `main-draw.py`、`main-draw-长度图.py`；`简化模型推导/`（Maple 推导 `推导.mw`、`推导带离心势.mw` 及 Word 推导笔记） |
| `运动稳定性-平衡态分析/` | 第 3 章相关的早期研究（2021–2023） | 天梯平衡态求解（不动点迭代、约束最小二乘）与雅可比矩阵稳定性分析 | `StableStates*.py`、`Lunar0829~0907.py`、`wlph.py`；`雅克比矩阵.mw`；推导笔记 `平衡方程.docx`、`珠点模型方程（极坐标）.docx`、`寻找稳定态.docx` |
| `早期版本/` | — | 最早上传的两个自由摆动程序（圆轨道 / 椭圆轨道、地心坐标） | `月球公转-天梯自由摆动-无for版-*.py` |

同一程序的多个日期版本（如 `main-虹吸链条-20240322.py`、`-20240324.py`、`-0426.py`）都保留了，文件名日期最新的一般是最终版。

## 运行方法

- 环境：Python 3.9（更新版本一般也可以），需要 `numpy`、`scipy`、`matplotlib`。
  ```
  pip install numpy scipy matplotlib
  ```
- 先进入程序所在的文件夹再运行，程序按相对路径读写文件。
- 算例参数放在同目录的 `Setting-算例X.txt` 里（节点数 `N`、阻尼 `c`、天梯末端位置 `x_end`、初始偏角 `theta_0`、月球轨道偏心率 `epsilon` 等）。程序开头的 `name='例X'` 指定读哪一个算例；部分程序改成了 `name=sys.argv[1]`，用命令行参数传入：
  ```
  python main-虹吸链条-0322-例1.py 例2.1
  ```
- 流程：主程序积分求解，用 `np.savetxt` 把结果存成 txt；画图 / 动画程序再用 `np.loadtxt` 读回来出图。`Earth.png`、`moon*.png`、`moon2.gif` 是画图用的地球、月球贴图。
- `批量运行` 是 Spyder 控制台里的 `runfile(...)` 命令清单，里面还是当年本机的 `D:/米凯尼克斯/...` 绝对路径，用之前要改成自己的路径。
- `.mw` 文件是 Maple 工作表（静力学和简化模型的公式推导）。

## 没有上传的内容

- `图像结果/`、`结果绘图/` 里的 SVG/PNG 图像（约 1.36 GB）和数值结果
- 第 2 章的绘图 PPT、单独的结果截图和动画
- Maple 自动备份（`.bak`）和 `__pycache__`

## 相关成果

- Wang W W, Wu Z G, Liu J F. Conceptual Design and Mechanical Analysis of a Lunar Anchored Cislunar Tether. *Cosmic Research*, 2023, 61(1): 80-89.（第 2 章）
- Wang W W, Wu Z G, Liu J F. Dynamic analysis of cislunar suspension tether swings. *Acta Astronautica*, 2024, 216: 350-369.（第 3 章）
- Wu Z G, Wang W W, Liu L, Liu J F. Dynamic simulation of cargo transport along cislunar suspension tether by single node coupling model. *Acta Astronautica*, 2025, 226(1): 262-274.（第 4 章）
- 发明专利（已受理）：基于引力虹吸效应的自动运行月球天梯建模方法及系统，受理号 202510168560.6（第 5 章）

## 视频

- 第 5 章虹吸天梯仿真动画：https://www.bilibili.com/video/BV1hA2BYEEaH/
- 月球天梯视频：https://www.bilibili.com/video/BV1mJaZ64EtE/
