import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


def remove_file(filedir):
    if os.path.exists(filedir):
        os.system("rm " + filedir)

def create_directories(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
