import matplotlib.pyplot as plt
import statsmodels.graphics.gofplots as st
import numpy as np



def main():
    origin = "Z://Matlab_Data//training_data//chr_1_training.csv"

    x = np.random.normal(loc=8.5, scale=2.5, size=37)
    y = np.random.normal(loc=8.0, scale=3.0, size=37)
    pp_x = sm.ProbPlot(x)
    pp_y = sm.ProbPlot(y)
    st.qqplot_2samples(pp_x, pp_y)



main()