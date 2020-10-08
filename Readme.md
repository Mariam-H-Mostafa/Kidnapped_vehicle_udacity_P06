
**Kidnapped Vehicle Project**

In this project, a 2-D particle filter implemented in C++ is used to localize the vehicle. The particle filter will be given a map and some initial localization information (analogous to what a GPS would provide). At each time step your filter will also get observation and control data.
At first, we should initialize the filter…

![enter image description here](https://i.ibb.co/10sgmWw/Picture1.png)

Then we should predict the new particle location and add random Gaussian noise…

![enter image description here](https://i.ibb.co/cTcq1fN/Picture3.png)

In data association, I used nearest neighbor data association, and assign each sensor observation the map landmark ID associated with it.

![enter image description here](https://i.ibb.co/S3mDpCq/Picture4.png)

Next we will update the particle weight, first convert all observations from vehicle coordinate to map coordinate, then check the map landmarks that should be considered based on sensor range with respect to each particle location. In this function I will call SetAssociations and dataAssociation then update the weights of each particle using a mult-variate Gaussian distribution.

![enter image description here](https://i.ibb.co/gTy48H2/Picture5.png)

![enter image description here](https://i.ibb.co/270HVCX/Picture6.png)

Then Resample particles with replacement with probability proportional to weight.

![enter image description here](https://i.ibb.co/Wn5r94d/Picture7.png)
