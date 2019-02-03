# Run Away Robot with Unscented Kalman Filter Bonus Challenge Starter Code
Self-Driving Car Engineer Nanodegree Program

---

### Overview

This repository contains all the code needed to complete the Bonus Challenge: Catch the Run Away Car with Unscented Kalman Filter.

### Project Introduction

In this project, not only do you implement an UKF, but also use it to catch an escaped car driving in a circular path. 
The run away car will be being sensed by a stationary sensor, that is able to measure both noisy lidar and radar data. The capture vehicle will need to use these measurements to close in on the run away car. To capture the run away car the capture vehicle needs to come within .1 unit distance of its position. However the capture car and the run away car have the same max velocity, so if the capture vehicle wants to catch the car, it will need to predict where the car will be ahead of time.

### Running the Code

This project involves the Term 2 Simulator which can be downloaded [here](https://github.com/udacity/self-driving-car-sim/releases)


## Dependencies

* cmake >= v3.5
* make >= v4.1
* gcc/g++ >= v5.4
* uWebSocketIO

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./UnscentedKF` 


## Implementation
1. For UKF, See `ukf.cpp`
2. For target position estimation, see `main.cpp`
