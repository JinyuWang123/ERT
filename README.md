# Statistical Inference on Grayscale Images via the Euler-Radon Transform
Tools from topological data analysis have been widely used to represent binary images in many scientific applications. Methods that aim to represent grayscale images (i.e., where pixels instead take on continuous values) have been relatively underdeveloped. 
In this paper, we introduce the Euler-Radon transform, which generalizes Euler characteristic transform to grayscale images by using o-minimal structures and Euler integration over definable functions. Coupling the Karhunen–Loève expansion with our proposed topological representation, we offer hypothesis-testing algorithms based on the $\chi^2$ distribution for detecting significant differences between two groups of grayscale images. We illustrate our framework via extensive numerical experiments and simulations.

A Matlab implementation accompanying the paper "Statistical Inference on Grayscale Images via the Euler-Radon Transform".

The function [test_final.m](https://github.com/JinyuWang123/ERT/blob/main/test_final.m) implement the matching algorithm described in section 8 of "Statistical Inference on Grayscale Images via the Euler-Radon Transform".

## Compatibility and dependencies:
Developed and tested with Matlab R2023a on a Windows (x64) Machine.

## Disclaimer: 
The code for implementing our Algorithm is built upon the materials in the GitHub repository of authors of Crawford et al. (2020). We express our gratitude to Henry Kirveslahti for generously sharing the computer code essential for the implementation of the LECT and SELECT framework presented in ["Representing fields without correspondences: the lifted euler characteristic transform"](https://link.springer.com/article/10.1007/s41468-023-00133-w). We modified the code to suit the analysis in the section 8.

The code is provided as-is for academic use only and without any guarantees. Please contact the authors to report any bugs. Written by Jinyu Wang <jinyu_wang@brown.edu>.
