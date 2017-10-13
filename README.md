# Kidnapped Vehicle Project 
---

This is my Project 3 in term 2 of the Ucatity Self-Driving Car Engineer Nanodegree Program, which aims to utilize the final particle filter project for the Localization course to estimate the state of a kidnapped, moving vehicle based on a CTRV model, with noisy measurements using C++. In this project I chose Ubuntu 14.04 LTS for implementation. 

My robot has been kidnapped and transported to a new location! Luckily it has a map of this location, a (noisy) GPS estimate of its initial location, and lots of (noisy) sensor and control data. In this project I will implement a 2 dimensional particle filter in C++. My particle filter will be given a map and some initial localization information (analogous to what a GPS would provide). At each time step my filter will also get sensor observation(p_x,p_y,p_theta) and control data(velocity and yaw_rate). 

[image1]: ./PFresult.png

---

## Content of my project
---
This project origins from the [Udacity Starter Code repository](https://github.com/udacity/CarND-Kidnapped-Vehicle-Project),which includes [uWebSocketIO](https://github.com/uWebSockets/uWebSockets) installation files __install-ubuntu.sh__ and the directory structure of this repository is as follows:

```
root
|   build.sh
|   clean.sh
|   CMakeLists.txt
|   README.md
|   run.sh
|
|___data|  
|   |   map_data.txt|   
|   
|___src
    |   helper_functions.h
    |   main.cpp
    |   map.h
    |   particle_filter.cpp
    |   particle_filter.h
```

Data: `map_data.txt` includes the position of landmarks (in meters) on an arbitrary Cartesian coordinate system. Each row has three columns:(x position, y position, landmark id). All other data the simulator provides, such as observations and controls.: Map data provided by 3D Mapping Solutions GmbH.

Note the only file I should modify is **`particle_filter.cpp`** in the **`src`** directory, or basic requirements in [project rubrics](https://review.udacity.com/#!/rubrics/747/view). And the main protcol that main.cpp uses is the web socket server **uWebSocketIO** acting as a host, connecting the C++ programs to the Unity simulator, where the input values is provided by the simulator to the c++ program, and output values provided by the c++ program to the simulator. The file contains the scaffolding of a **`ParticleFilter`** class and some associated methods. Read through the code, the comments, and the header file **`particle_filter.h`** to get a sense for what this code is expected to do. And I did take a look at **`src/main.cpp`** as well, which contains the code that will actually be running your particle filter and calling the associated methods.

The flow of Particle Filter is Initialization -> Prediction -> Update ( update weight with data association) -> Resampling on this synthetic dataset, so I can narrow down the source of error when run into some problems.

---

## How to run this project

This project involves the Term 2 Simulator which can be downloaded [here](https://github.com/udacity/self-driving-car-sim/releases)

This repository includes two files that can be used to set up and intall uWebSocketIO for either Linux or Mac systems. For windows you can use either Docker, VMware, or even Windows 10 Bash on Ubuntu to install uWebSocketIO.

Once the install for uWebSocketIO is complete, the main program can be built and ran by doing the following from the project top directory.

1. mkdir build
2. cd build
3. cmake ..
4. make
5. ./particle_filter

Alternatively some scripts have been included to streamline this process, these can be leveraged by executing the following in the top directory of the project:

1. ./clean.sh
2. ./build.sh
3. ./run.sh

Tips for setting up your environment can be found [here](https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/0949fca6-b379-42af-a919-ee50aa304e6a/lessons/f758c44c-5e40-4e01-93b5-1a82aa4e044f/concepts/23d376c7-0195-4276-bdf0-e02f1f3c665d)

Note that the programs that need to be written to accomplish the project are src/particle_filter.cpp. The program main.cpp has already been filled out, but feel free to modify it.

---

### Here is the main protcol that main.cpp uses for uWebSocketIO in communicating with the simulator.

**INPUT: values provided by the simulator to the c++ program**
You can find the inputs to the particle filter in the `data` directory. 

- // sense noisy position data from the simulator:   ["sense_x"] , ["sense_y"] , ["sense_theta"] 

- // get the previous velocity and yaw rate to predict the particle's transitioned state:   ["previous_velocity"], ["previous_yawrate"]

- // receive noisy observation data from the simulator, in a respective list of x/y values:   ["sense_observations_x"] , ["sense_observations_y"] 


**OUTPUT: values provided by the c++ program to the simulator**

- // best particle values used for calculating the error evaluation: ["best_particle_x"], ["best_particle_y"],["best_particle_theta"] 

- //Optional message data used for debugging particle's sensing and associations

- // for respective (x,y) sensed positions ID label: ["best_particle_associations"]

- // for respective (x,y) sensed positions:  ["best_particle_sense_x"] <= list of sensed x positions, ["best_particle_sense_y"] <= list of sensed y positions

My job is to build out the methods in `particle_filter.cpp` until the simulator output says:

```
Success! Your particle filter passed!

```
---


## Simulation results
---

![alt text] [image1]

---


## Code Style
---
Here I stick to [Google's C++ style guide](https://google.github.io/styleguide/cppguide.html). And I use the ubuntu file editor gedit with the similar user interface with emacs. 


## Project Instructions and Rubric
---
### Submission
All you will submit is your completed version of `particle_filter.cpp`, which is located in the `src` directory. You should probably do a `git pull` before submitting to verify that your project passes the most up-to-date version of the grading code (there are some parameters in `src/main.cpp` which govern the requirements on accuracy and run time.)

### Success Criteria
If your particle filter passes the current grading code in the simulator (you can make sure you have the current version at any time by doing a `git pull`), then you should pass! 

The things the grading code is looking for are:

1. **Accuracy**: your particle filter should localize vehicle position and yaw to within the values specified in the parameters `max_translation_error` and `max_yaw_error` in `src/main.cpp`.

2. **Performance**: your particle filter should complete execution within the time of 100 seconds.
---
## How to write a README
A well written README file can enhance your project and portfolio.  Develop your abilities to create professional README files by completing [this free course](https://www.udacity.com/course/writing-readmes--ud777).

