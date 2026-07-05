# Non-Cartesian-SSDU-MRI

<p align="center">
  <img src="./Data/traj_figure.png" width="900">
</p>

## Abstract

Physics-driven deep learning (PD-DL) methods have shown strong promise for accelerated MRI reconstruction by combining learned image priors with the physical MRI forward model. However, conventional PD-DL pipelines can struggle when the target images contain large phase variations, since these phase inconsistencies may degrade data consistency and reconstruction quality. This issue is especially relevant in multi-echo fMRI, where substantial phase differences can appear across echo images.

To address this challenge, we propose a low-passed filter phase-corrected physics-driven deep learning method, referred to as LPF-based PC-PDDL, that explicitly models and corrects phase inconsistencies within the unrolled reconstruction process. 

Experiments on prospectively undersampled non-Cartesian multi-echo fMRI data show that the low-pass filtering strategy improves upon standard non-phase-corrected PD-DL despite its simplicity. These results suggest that explicit phase correction can improve non-Cartesian multi-echo fMRI reconstruction, while also highlighting the importance of robust and generalizable phase estimation within unrolled PD-DL frameworks.

<p align="center">
  <img src="./Data/tsnr_figure.png" width="900">
</p>

```text
Non-Cartesian-SSDU-MRI/
└── Data/
   ├── traj_figure.png
   └── tsnr_figure.png

└── Preprocess/
   ├── Espirit/
   ├── nufft/
   ├── A0_DCF.py
   ├── A0_MaskGen.m
   ├── A1_FastOperator.m
   ├── A2_LabelCoilGen.m
   ├── A3_EHyGen.m
   ├── A4_EHyGen_TestfMRI.m
   └── README.md

└── Training/
   ├── DC.py
   ├── DataLoader.py
   ├── DataLoader_Test.py
   ├── ResNet.py
   ├── Unrolled_Network.py
   ├── train.py
   ├── test.py
   ├── requirements.txt

└── README.md
```


## Quick Start
Note: This code was tested with `torch==2.2.1+cu121`. 

## Installation

**Note:** This code was tested with `torch==2.2.1+cu121`.

### 1. Clone this repository

```bash
git clone https://github.com/MahdiSaberii/Non-Cartesian-SSDU-MRI.git
cd Non-Cartesian-SSDU-MRI
```

### 2. Create and activate conda environment

```bash
conda create -n non_cartesian_ssdu_mri python=3.10 -y
conda activate non_cartesian_ssdu_mri
```

### 3. Install PyTorch

```bash
pip install torch==2.2.1+cu121 --index-url https://download.pytorch.org/whl/cu121
```

### 4. Install remaining requirements

```bash
pip install -r requirements.txt
```
