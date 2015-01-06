Block-Modeling
==============

Creates block models out of social tie data from various algorithms, including CONCOR, random generation, and an algorithm for optimizing from a previous block model.

==============

This program uses the numpy, math, copy, and wxPython modules, so you will need to install them in order to use the program. 

Thorough precautions against improper usage have not been made, so there may not be elegant responses to mistaken inputs. 

FUNCTION: This program is to be used for the purposes of social network analysis, specifically blockmodeling. It is capable of using a number of algorithms to generate block models from social network data. This data should take the form of matrices corresponding to ties between individuals in the network. These data matrices form what is referred to as the stack. The T function developed by Boorman and Levitt is used to compare the effectiveness of various block models in describing the social structure.

USAGE: To launch the program, the main.py file should be called from terminal. The user will first be prompted to specify how many blocks and how many matrices are desired. After inputting this information, a new frame will emerge in which the user can enter the matrix data (one matrix at a time). Once all of this information is inputted, the user can then choose from the button panel which algorithm they would like to use to generate the block model. 

Notes on optimization: 1. Optimization should only be used after a previous block model has been created. 2. Based on current block model, the optimzation function considers moving each agent one by one to a new block, making changes that would increase the optimization function value (in this case T). This process is repeated until a local maximum is reached. 3.Right now, this function works faster than the included "optimization by parts" function, which aims to improve on the number of computations by only considering the parts of the optimization function affected by the agent under consideration. However, it is not clear that this approach would actually work faster, since even calculating part of the optimization function still requires parsing through an asymptotically equivalent number of steps (i.e., it appears that the approaches would be Big-O of each other). The function still works slower than desired, though,so in the future working on speeding up this process would be beneficial.
