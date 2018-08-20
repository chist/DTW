# DTW
Matlab mex-function to calculate DTW distance between two 2D-curves and to find corresponding points.

# Installation
Open the repo directory in MATLAB and type

```mex dtw_c.c```

# Description

```[dist, ix, iy] = dtw_c(x, y, window);```

## Input arguments

* `x` is Mx2 array of points (first curve)
* `y` is Nx2 array of points (second curve)
* `window` (optional) is a parameter used to constrain the maximum time difference between corresponding points of different curves (setting this parameter reduces algorithm running time). By default `window` is set to the `min(M, N)`.

## Output

* `dist` is normalized DTW distance between curves
* `ix` is a row-vector containing ordering numbers of points from the first curves, that correspond to `iy`-points from the second curve
* `iy` is a row-vector containing ordering numbers of points from the second curves that correspond to `ix`-points from the first curve 

# Usage example

```
x = [1:1000; cos(1:1000)]';
y = [1:5:1000; sin(1:5:1000)]';

dist = dtw_c(x, y, 10);
fprintf('DTW distance is %0.8g\n', dist);
dist = dtw_c(x, y);
fprintf('DTW distance is %0.8g\n', dist);
```

# Reference

[1] Wikipedia: [Dynamic Time Warping](https://en.wikipedia.org/wiki/Dynamic_time_warping)

[2] [Understanding DTW algorithm](https://www.psb.ugent.be/cbd/papers/gentxwarper/DTWalgorithm.htm)
