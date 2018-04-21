# Triangle pursuit: About some recurrent sequences

<center>
<img src="https://ahstat.github.io/images/2017-6-11-Triangle-pursuit/rotation_homothety/rot_onenorm_700.png" alt="Rotation and one-norm" width="31%"/>

<img src="https://ahstat.github.io/images/2017-6-11-Triangle-pursuit/rotation_homothety/rot_eucnorm_700.png" alt="Rotation and Euclidian norm" width="31%"/>

<img src="https://ahstat.github.io/images/2017-6-11-Triangle-pursuit/rotation_homothety/rot_maxnorm_700.png" alt="Rotation and maximum norm" width="31%"/>
</center>

**Triangle pursuit** is a program computing recurrent sequences as described in this blog post: https://ahstat.github.io/Triangle-pursuit/ It also offers generalization for different rules, different norms, larger number of initial points and higher dimensions.

**Documentation**

The `main.R` file is self-documented, so you can follow it line by line.
You can keep option `verbose = TRUE` to understand each function through multiple examples.

**Some outputs**

<center>
<img src="https://ahstat.github.io/images/2017-6-11-Triangle-pursuit/map/begin_2pi.png" alt="Initial elements" width="49%"/>
<img src="https://ahstat.github.io/images/2017-6-11-Triangle-pursuit/map/end_2pi.png" alt="Resulting elements" width="49%"/>
</center>

*Mapping described in the blog post and plotted through this program.*

We restrict the mapping on the interval $$(0, \pi)$$ and show a more detailed plot in Fig. 10. Notice that triangle corresponding to $$t = \pi / 3 \approx 1.05$$ remains unchanged by the mapping.

<center>
<img src="https://ahstat.github.io/images/2017-6-11-Triangle-pursuit/map/begin_pi.png" alt="Initial elements" width="49%"/>
<img src="https://ahstat.github.io/images/2017-6-11-Triangle-pursuit/map/end_pi.png" alt="Resulting elements" width="49%"/>
</center>

*Mapping described in the blog post and plotted through this program with additional indications.*

