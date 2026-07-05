## Preprocessing for Multi-Mask SSDU Learning

This folder contains the preprocessing pipeline used to prepare prospectively undersampled non-Cartesian multi-echo fMRI data for multi-mask SSDU learning. The preprocessing includes trajectory handling, density compensation, NUFFT operator preparation, coil sensitivity map estimation, mask generation, and creation of the data terms used for self-supervised reconstruction.

The goal of this preprocessing stage is to convert raw non-Cartesian k-space measurements into structured inputs that can be used by the downstream physics-driven deep learning reconstruction pipeline. In the multi-mask SSDU setting, the acquired measurements are split into multiple training and validation subsets, allowing the model to learn from prospectively undersampled data without requiring fully sampled reference images.

### Overview

Physics-driven reconstruction of non-Cartesian MRI requires several acquisition-specific components before training can begin. In this pipeline, the measured k-space data are combined with the spiral trajectory, density compensation function, coil sensitivity maps, and NUFFT operators to generate the inputs needed for SSDU-style learning.

The preprocessing workflow performs the following steps:

1. Prepare non-Cartesian k-space trajectories and density compensation functions.
2. Generate multiple SSDU training and validation masks.
3. Build NUFFT-based forward and adjoint operators.
4. Estimate coil sensitivity maps using ESPIRiT.
5. Generate coil-combined label images or adjoint reconstructions.
6. Save the processed files for multi-mask SSDU training and evaluation.

## Multi-Mask SSDU Preprocessing

In SSDU learning, the acquired k-space samples are divided into two disjoint sets:

- **Training set**: used inside the data-fidelity term during reconstruction.
- **Validation/loss set**: used to compute the self-supervised loss.

Instead of using only one split, this preprocessing pipeline generates multiple train/validation masks. These multiple masks improve the diversity of the self-supervised learning signal and reduce the dependence of training on a single random k-space partition. For each mask, the code prepares the corresponding measurements and operators needed by the reconstruction network.

## Folder Structure

```text
Preprocess/
│
├── Espirit/
│   └── ESPIRiT coil sensitivity estimation utilities
│
├── nufft/
│   └── NUFFT implementation and related utilities
│
├── A0_DCF.py
│   └── Prepares density compensation functions and trajectory-related files
│
├── A0_MaskGen.m
│   └── Generates SSDU training and validation masks
│
├── A1_FastOperator.m
│   └── Builds NUFFT-based forward and adjoint operators
│
├── A2_LabelCoilGen.m
│   └── Generates coil sensitivity maps and coil-combined label images
│
├── A3_EHyGen.m
│   └── Computes adjoint NUFFT reconstructions from measured k-space data
│
├── A4_EHyGen_TestfMRI.m
│   └── Applies the preprocessing pipeline to test fMRI data
│
└── README.md
