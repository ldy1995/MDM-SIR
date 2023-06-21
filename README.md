# MDM-SIR (SIR-E3DTV)
This is the code of the iterative reconstruction method at Section III-E in [paper](https://doi.org/10.1109/TMI.2020.3001616)  [MDM-PCCT](https://github.com/ldy1995/MDM-PCCT) 

## Requirements
* Linux Platform
* MATLAB
* [IRT](https://web.eecs.umich.edu/~fessler/code/)

## Describe
This model utilize an advanced prior term (enhanced 3DTV, E3DTV) for spectral and photon-counting data reconstruction. To solve the model, we use the alternating direction method of multipliers algorithm (ADMM).

## Results
The following image shows the pig leg results of the compared methods, i.e., FBP, 3DTV, BM3D and the presented MDM-SIR methods, respectively. And zoomed ROIs indicated by the red box illustrate the details of the phantom.

![Fig. 1. Results of the compared methods.](./data/img_compared_method_real_data.png)

## Citation
If you find this work or code is helpful in your research, please cite:
```
@ARTICLE{9115052,
  author={Li, Danyang and Zeng, Dong and Li, Sui and Ge, Yongshuai and Bian, Zhaoying and Huang, Jing and Ma, Jianhua},
  journal={IEEE Transactions on Medical Imaging}, 
  title={MDM-PCCT: Multiple Dynamic Modulations for High-Performance Spectral PCCT Imaging}, 
  year={2020},
  volume={39},
  number={11},
  pages={3630-3642},
  doi={10.1109/TMI.2020.3001616}}
```

## Contact
If you have any question, please feel free to contact Danyang Li (Email: lidanyang1995@smu.edu.cn or dyli0730@gmail.com).
