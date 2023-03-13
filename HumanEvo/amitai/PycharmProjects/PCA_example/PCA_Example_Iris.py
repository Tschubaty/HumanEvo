import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

NAN = "nan"


def main():
    # origin = "C://Users//amita//Documents//University//Year_C//SemesterA//project//Examples//meth_window_example.csv"
    origin = "Z://Matlab_Data//chr_1_united.csv"

    df = pd.read_csv(origin, names=['p_0708', 'p_1053', 'p_1116', 'p_1496', 'p_1507', 'p_1631', 'p_2520'])

    print(df)

    df = remove_NANS(df)

    print("nans removed")
    print(df)
    df = best_Var(df)

    print("best var")

    print(df)



    x = df.values

    # amputate.fit_transform(x)

    print(x)

    pca = PCA(n_components=2)

    principalComponents = pca.fit_transform(x)
    principalDf = pd.DataFrame(data=principalComponents
                               , columns=['principal component 1', 'principal component 2'])

    print(principalDf)

    headers = pd.DataFrame({ 'humans': ['p_0708', 'p_1053', 'p_1116', 'p_1496', 'p_1507', 'p_1631', 'p_2520'] })

    finalDf = pd.concat([principalDf, headers], axis=1)

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlabel('Principal Component 1', fontsize=15)
    ax.set_ylabel('Principal Component 2', fontsize=15)
    ax.set_title('2 component PCA', fontsize=20)
    targets = ['p_0708', 'p_1053', 'p_1116', 'p_1496', 'p_1507', 'p_1631', 'p_2520']
    colors = ['r', 'g', 'b', 'c', 'y', 'm', 'k']
    for target, color in zip(targets, colors):
        indicesToKeep = finalDf['humans'] == target
        ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1']
                   , finalDf.loc[indicesToKeep, 'principal component 2']
                   , c=color
                   , s=50)
    ax.legend(targets)
    ax.grid()

    plt.show()


def remove_NANS( df ):
    return df.dropna(axis=0)


def best_Var(df):
    temp = df.reindex(df.var(axis=1).sort_values().index, axis=0)[::-1]
    return temp[:1500]

if __name__ == '__main__':
    main()
