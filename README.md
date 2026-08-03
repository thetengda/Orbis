# GREAT-PIFGO: PPP/INS Factor Graph Optimization-Based Positioning and Navigation Software by Wuhan University GREAT Group (Version 1.0)

## Overview

  The GREAT (GNSS+ **RE**search, **A**pplication and **T**eaching) software suite is designed and developed by the School of Geodesy and Geomatics, Wuhan University. It is a comprehensive platform for space geodesy data processing, precise positioning and orbit determination, as well as multi-source fusion navigation. <br />
  GREAT-PIFGO is a key module within the GREAT software, dedicated to Factor Graph Optimization (FGO) based navigation solutions. It supports a variety of algorithms including PPP, PPP/INS, and multi-sensor fusion. In the software, core computation modules are implemented in C++, while auxiliary Python 3 scripts are provided for plotting results. GREAT-PIFGO uses CMAKE for build management, allowing users to flexibly choose mainstream C++ compilers such as GCC, Clang, and MSVC. It currently supports build and execution on Windows; for Linux, users are encouraged to compile and test locally. <br />
  GREAT-PIFGO consists of two portable program libraries: LibGREAT and LibGnut. In addition to the GNSS positioning solutions from the original GREAT-PVT and the multi-sensor fusion navigation capabilities provided by GREAT-MSF, GREAT-PIFGO delivers PPP and tightly-coupled PPP/INS algorithms based on factor graph optimization. <br />
  This open-sourced GREAT-PIFGO version 1.0 supports the following capabilities:

1.Supports PPP positioning based on factor graph optimization 

2.Supports tightly-coupled PPP/INS integration based on factor graph optimization

3.Supports satellite navigation systems including GPS, Galileo, BDS-2/3，GLONASS

4.Supports custom IMU data formats and noise models

5.Support for trajectory visualization and Google Maps viewing

6.The software package also provides plotting scripts for positioning results to facilitate user data analysis


## Package Directory Structure

```shell
GREAT-PIFGO
  ./src                   Source code *
    ./app                 Main programs of GREAT-PVTFGO, GREAT-GINSFGO, GREAT-MSF and GREAT-PVT *
    ./LibGREAT            Core algorithm library for Factor Graph Optimization*
    ./LibGnut             Gnut library *
    ./Third-party         Third-party libraries *
  ./sample_data           Example datasets *
  ./doc                   Documentation related to GREAT-PIFGO *
```

## Installation and Usage

See **GREAT-PIFGO Documentation.pdf** included in the **./doc** directory.


## Contributing

Developers:

Wuhan University GREAT Team, Wuhan University.

Third-party libraries:

* GREAT-PIFGO uses the G-Nut library ([http://www.pecny.cz](http://www.pecny.cz))
  Copyright (C) 2011–2016 GOP - Geodetic Observatory Pecny, RIGTC.

* GREAT-PIFGO uses the pugixml library ([http://pugixml.org](http://pugixml.org))
  Copyright (C) 2006–2014 Arseny Kapoulkine.

* GREAT-PIFGO uses the Newmat library ([http://www.robertnz.net/nm_intro.htm](http://www.robertnz.net/nm_intro.htm))
  Copyright (C) 2008: R B Davies.

* GREAT-PIFGO uses the spdlog library ([https://github.com/gabime/spdlog](https://github.com/gabime/spdlog))
  Copyright (C) 2015–present, Gabi Melman & spdlog contributors.

* GREAT-PIFGO uses the GLFW library ([https://www.glfw.org](https://www.glfw.org))
  Copyright (C) 2002–2006 Marcus Geelnard, Copyright (C) 2006–2019 Camilla Löwy

* GREAT-PIFGO uses the Eigen library ([https://eigen.tuxfamily.org](https://eigen.tuxfamily.org))
  Copyright (C) 2008–2011 Gael Guennebaud

* GREAT-PIFGO uses the PSINS library ([https://psins.org.cn](https://psins.org.cn))
  Copyright (c) 2015–2025 Gongmin Yan

* GREAT-PIFGO uses the Ceres library ([http://ceres-solver.org](http://ceres-solver.org))
  Copyright 2023 Google Inc.

## Download

GitHub: [https://github.com/GREAT-WHU/GREAT-PIFGO](https://github.com/GREAT-WHU/GREAT-PIFGO)

## Others

You are welcome to join the QQ group (1009827379) for discussion and exchange.

WeChat Official Account: **GREAT智能导航实验室** — we will continue to share team updates.

bilibili Account: **GREAT智能导航实验室** — we will continue to publish software walkthrough videos.

---





# GREAT-PIFGO: 武汉大学GREAT团队基于因子图优化的PPP/INS定位导航解算软件（1.0版）

## 概述

&emsp;&emsp;GREAT (GNSS+ REsearch, Application and Teaching) 软件由武汉大学测绘学院设计开发，是一个用于空间大地测量数据处理、精密定位和定轨以及多源融合导航的综合性软件平台。<br />
&emsp;&emsp;GREAT-PIFGO是GREAT软件中的一个重要模块，主要用于因子图优化 (Factor Graph Optimization) 导航解算，包括PPP、PPP/INS与多传感器融合等多种算法。软件中，核心计算模块使用C++语言编写，辅助脚本模块使用Python3语言实现结果绘制。GREAT-PIFGO软件使用CMake工具进行编译管理，用户可以灵活选择GCC、Clang、MSVC等主流C++编译器。目前支持在Windows下编译运行，Linux系统需要用户自行编译测试。<br />
&emsp;&emsp;GREAT-PIFGO由2个可移植程序库组成，分别是LibGREAT和LibGnut。除了原GREAT-PVT中的GNSS定位解决方案、GREAT-MSF提供的多传感器融合导航功能外，GREAT-PIFGO提供了基于因子图优化的PPP与PPP/INS紧耦合算法。<br />
&emsp;&emsp;本次开源的GREAT-PIFGO 1.0版本支持以下功能：

	1.支持基于因子图优化的RTK解算
	
	2.支持基于因子图优化的PPP/INS紧耦合解算
	
	3.支持GPS、Galileo、BDS-2/3、GLONASS卫星导航系统
	
	4.支持自定义IMU数据格式、噪声模型
	
	5.支持轨迹动态显示与谷歌地图查看
	
	6.软件包还提供定位结果绘图脚本，便于用户对数据进行结果分析



## 软件包目录结构
```shell
GREAT-PIFGO
  ./src                     源代码 *
    ./app                  GREAT-PVTFGO、GREAT-GINSFGO、GREAT-MSF和GREAT-PVT主程序 *
    ./LibGREAT             因子图优化核心算法库 *
    ./LibGnut              Gnut库 *
    ./Third-party          第三方库 *
  ./sample_data          算例数据 *
  ./doc                  GREAT-PIFGO相关文档 *
```

## 安装和使用

参见**doc**文件夹中的《GREAT-PIFGO说明文档 1.0.pdf》或关注我们团队后续在视频网站bilibili发布的讲解视频。



## 参与贡献

开发人员：

武汉大学GREAT团队, Wuhan University.

三方库：

* GREAT-PIFGO使用G-Nut库(http://www.pecny.cz)
  Copyright (C) 2011-2016 GOP - Geodetic Observatory Pecny, RIGTC.
  
* GREAT-PIFGO使用pugixml库(http://pugixml.org)
  Copyright (C) 2006-2014 Arseny Kapoulkine.

* GREAT-PIFGO使用Newmat库(http://www.robertnz.net/nm_intro.htm)
  Copyright (C) 2008: R B Davies.

* GREAT-PIFGO使用spdlog库(https://github.com/gabime/spdlog)
  Copyright(C) 2015-present, Gabi Melman & spdlog contributors.

* GREAT-PIFGO使用GLFW库(https://www.glfw.org)
  Copyright (C) 2002-2006 Marcus Geelnard, Copyright (C) 2006-2019 Camilla Löwy

* GREAT-PIFGO使用Eigen库(https://eigen.tuxfamily.org)
  Copyright (C) 2008-2011 Gael Guennebaud

* GREAT-PIFGO使用PSINS库(https://psins.org.cn)
  Copyright(c) 2015-2025 Gongmin Yan

* GREAT-PIFGO使用Ceres 库(http://ceres-solver.org)
  Copyright 2023 Google Inc.
## 下载地址

GitHub：https://github.com/GREAT-WHU/GREAT-PIFGO

## 其它

欢迎加入QQ群(1009827379)参与讨论与交流。

微信公众号：GREAT智能导航实验室，我们将持续推送团队成果。

bilibili账号：GREAT智能导航实验室，我们将持续发布软件讲解视频。
