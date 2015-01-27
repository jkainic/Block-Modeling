Block-Modeling
==============

Creates block models out of social tie data from various algorithms, including CONCOR, random generation, and an algorithm for optimizing from a previous block model.

==============

This program uses the numpy, math, copy, and wxPython modules, so you will need to install them in order to use the program. 

Thorough precautions against improper usage have not been made, so there may not be elegant responses to mistaken inputs. 

FUNCTION: This program is to be used for the purposes of social network analysis, specifically blockmodeling. It is capable of using a number of algorithms to generate block models from social network data. This data should take the form of matrices corresponding to ties between individuals in the network. These data matrices form what is referred to as the stack. The T function developed by Boorman and Levitt is used to compare the effectiveness of various block models in describing the social structure.

USAGE: To launch the program, the main.py file should be called from terminal. The user will first be prompted to specify how many blocks and how many matrices are desired. After inputting this information, a new frame will emerge in which the user can enter the matrix data (one matrix at a time). Once all of this information is inputted, the user can then choose from the button panel which algorithm they would like to use to generate the block model. 

Notes on optimization: 1. Optimization should only be used after a previous block model has been created. 2. Based on current block model, the optimzation function considers moving each agent one by one to a new block, making changes that would increase the optimization function value (in this case T). This process is repeated until a local maximum is reached. 3.This is a greedy algorithm, so in the future working on speeding up this process would be beneficial. In particular, since the binary matrices used are typically sparse, it would be beneficial to use a different method of encoding the matrices. One attempt is in the process of being made, but has not been completed. This attempt involves encoding each matrix as a dictionary, where the keys are column numbers, and the associated values are lists of the indexes at which there are non-zero entries in the matrix. In this instantiation, there has also been attention paid to avoid expensive function calls where possible. This code is still in the process of being debugged, so speed improvement is not yet known.
