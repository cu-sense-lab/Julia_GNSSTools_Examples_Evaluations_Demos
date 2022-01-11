# SeNSe Lab Video Script

This script is generates videos and plots for the SeNSe lab video. For the constellation video. THere are two functions which produce three separate videos in the `make_plots.jl` script. These videos are later put together using OpenShot Video Editor for Linux. Julia 1.6.0-rc1 (unstable release was used). Note that the following was performed when starting Julia to run this script. This script also uses the GNSSTools master version from commit ID [2f233d13](http://192.168.3.66/bilardis/GNSSTools.jl/tree/2f233d1330bf7e47b1e74154e79395c02d6b1ebf)

```bash
export LD_PRELOAD=/usr/lib64/libstdc++.so.6
export JULIA_NUM_THREADS=8
/path/to/julia/binary
```

The two functions that you can use to regenerate these animations are:

- `make_plots`: makes the constellation animation
- `zoom_earth`: makes animation of zooming from the original view of the 1st constellation to a view of the second which is closer to the Earth

Note that the global variable `constellation` must be defined prior to running `make_plots`. Three sets of parameters, each succeeding the banners:

```
###########IRIDIUM############, 
###########STARLINK###########, and
#############GPS##############
```

Comment the parameters for constellations you do not want to generate and uncomment the parameters for the constellation that you want to generate.
