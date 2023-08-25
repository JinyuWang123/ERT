# Statistical Inference on Grayscale Images via the Euler-Radon Transform
Tools from topological data analysis have been widely used to represent binary images in many scientific applications. Methods that aim to represent grayscale images (i.e., where pixels instead take on continuous values) have been relatively underdeveloped. 
In this paper, we introduce the Euler-Radon transform, which generalizes Euler characteristic transform to grayscale images by using o-minimal structures and Euler integration over definable functions. Coupling the Karhunen–Loève expansion with our proposed topological representation, we offer hypothesis-testing algorithms based on the $\chi^2$ distribution for detecting significant differences between two groups of grayscale images. We illustrate our framework via extensive numerical experiments and simulations.

A Matlab implementation accompanying the paper "Statistical Inference on Grayscale Images via the Euler-Radon Transform".

The function [test_final.m](https://github.com/JinyuWang123/ERT/blob/main/test_final.m) implement the matching algorithm described in section 8 of "Statistical Inference on Grayscale Images via the Euler-Radon Transform".

## Compatibility and dependencies:
Developed and tested with Matlab R2023a on a Windows (x64) Machine.

## Disclaimer: 
The code is provided as-is for academic use only and without any guarantees. Please contact the authors to report any bugs. Written by Jinyu Wang <jinyu_wang@brown.edu>.
