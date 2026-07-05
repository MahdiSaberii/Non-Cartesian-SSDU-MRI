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
   ├── Tuke_01.mat
   ├── Tukey_32.mat
   ├── Tukey_64.mat
   ├── all_ones_mask.mat
   ├── trj_kxy_full.mat
   └── trj_kxy_single.mat

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
   └── test.py
   
└── requirements.txt
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

## Training and Testing

After completing the preprocessing pipeline, the prepared data can be used for multi-mask SSDU training. The preprocessing step should generate the required NUFFT operators, coil sensitivity maps, density compensation files, and SSDU training/validation masks.

Before starting training, set the desired number of training masks in `train.py`. This controls how many SSDU mask splits are used during multi-mask self-supervised learning.

To train the reconstruction network, run:

```bash
python train.py
```

After training is completed, the trained model can be tested on the full fMRI scan. Specify the path to the saved model weights in `test.py`, then run:

```bash
python test.py
```

The testing script loads the trained model weights and reconstructs the whole fMRI scan using the selected checkpoint.

## 📝 BibTeX

If you find this repository useful in your research, please consider citing our work:

```bibtex
@inproceedings{saberi2026phase,
  title={Phase-Correction Strategies for Physics-Driven Deep Learning Reconstruction of Accelerated Non-Cartesian Multi-Echo fMRI},
  author={Saberi, Mahdi and Yu, Zidan and Rettenmeier, Christoph and Stenger, Andrew and Ak{\c{c}}akaya, Mehmet},
  booktitle={2026 IEEE 23rd International Symposium on Biomedical Imaging (ISBI)},
  pages={1--4},
  year={2026},
  organization={IEEE}
}
