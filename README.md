# Signature-of-Geometric-Centroids
Signature of Geometric Centroids for 3D Local Shape Description and Partial Shape Matching


### Abstract
Depth scans acquired from different views may contain nuisances such as noise, occlusion, and varying point density. We propose a novel Signature of Geometric Centroids descriptor, supporting direct shape matching on the scans, without requiring any preprocessing such as scan denoising or converting into a mesh. First, we construct the descriptor by voxelizing the local shape within a uniquely defined local reference frame and concatenating geometric centroid and point density features extracted from each voxel. Second, we compare two descriptors by employing only corresponding voxels that are both non-empty, thus supporting matching incomplete local shape such as those close to scan boundary. Third, we propose a descriptor saliency measure and compute it from a descriptor-graph to improve shape matching performance. We demonstrate the descriptor's robustness and effectiveness for shape matching by comparing it with three state-of-the-art descriptors, and applying it to object/scene reconstruction and 3D object recognition.


![prediction example](https://tangbohu.github.io/SGC/images/desc_constru.png)
*Figure 1*: Constructing an SGC descriptor. (a) Construct a unique LRF from a spherical support centered at a feature point (in pink); (b) segment a cubical support centered at the feature point and aligned with the LRF; (c) voxelize the support and extract centroid features from non-empty voxels; the centroid color indicates point density in the voxel, where small and large densities are colored in blue and red respectively.



![prediction example1](https://tangbohu.github.io/SGC/images/match_pipeline.png)
*Figure 2*: Matching two scans using SGC descriptors: (a) sampled feature points (in purple) on two input scans (only part of samples are shown for clarity); (b) calculated LRFs and descriptors; (c) a pair of matched descriptors; (d) match the two scans based on aligning the associated LRFs; and (e) refine the scan alignment using ICP.



![prediction example2](https://tangbohu.github.io/SGC/images/app_reconstru.png)
*Figure 3*: Our reconstruction results. (a) Super Mario; (b) Frog; and (c) Stage scene.


### Code
please refer to https://tangbohu.github.io/SGC/index.html

### Citation
If you find our work useful in your research, please consider citing:

    @article {Tang-2016-SGC,
      title = {Signature of Geometric Centroids for 3D Local Shape Description and Partial Shape Matching},
      author = {Keke Tang and Peng Song and Xiaoping Chen},
      booktitle={Asian Conference on Computer Vision},
     pages={311--326},
     year = {2016}
     }

    @article {Tang-2017-SGCRec,
     title = {3D Object Recognition in Cluttered Scenes With Robust Shape Description and Correspondence Selection},
     author = {Keke Tang and Peng Song and Xiaoping Chen},
     journal={IEEE Access},
     volume={5},
     pages={1833--1845},
     year={2017}
    }
    
### License
Our code is released under MIT License (see LICENSE file for details).
