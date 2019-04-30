# CRYSTAL MET

Crystal MET (Crystal Maker Effective Tool) is a python3 programs that helps you building a system using only the informations taken from a CIF file.

Use it with care, always check the output, it may work differently from what you expect.

I'm not responsible of a missuse of the results given, neither am I of any collateral such as your computer bursting into flame or any physical incident that could result of an improper use of this tool.

## Dependencies

Crystal MET uses :

1. Scipy
2. Numpy

## Use

`python3 crystal_met.py INPUT_FILE`


## Performances 

The calculation increases exponentially as you increase the number of atoms in the system you want.
Here are the results of some tests I made :

![GitHub Logo](/images/fig1.png) ![GitHub Logo](/images/fig2.png)
