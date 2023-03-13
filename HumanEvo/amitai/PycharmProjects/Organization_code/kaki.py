import pandas as pd


"""This file is to sort and order the time vector and labels and to create repetitive matlab code"""







"""Sorted labels and time vector"""

names = ['I1116', 'I5319', 'I4532', 'I2862', 'I4529', 'I3758', 'I5725', 'I2861', 'I3957', 'I3565', 'I4792', 'I4775', 'I4315', 'I1053', 'I7421', 'I3255', 'I5835', 'I2514', 'I5748', 'I1633', 'I5950', 'I5838', 'I5742', 'I5743', 'I2105', 'I2520', 'I2935', 'I2978', 'I2980', 'I3133', 'I4634', 'I1965', 'I1631', 'I1632', 'I1961', 'LBK', 'I1496', 'I5077', 'I1962', 'I4438', 'I15', 'I2134', 'I1507', 'I4878', 'I4432', 'I4873', 'Los', 'LaB', 'I4596', 'I5233', 'I0708', 'I1960', 'I4914', 'I4875', 'I4877', 'I2139', 'SF1', 'I1734', 'I5236', 'I5235', 'M45', 'Ust']
time_vec = [1, 1, 1.2, 1.6, 1.8, 2.2, 2.5, 2.8, 3, 3.1, 3.3, 3.4, 3.5, 3.8, 3.8, 4, 4.2, 4.3, 4.4, 4.5, 4.5, 4.6, 4.7, 4.7, 4.9, 5.1, 5.1, 5.1, 5.2, 5.4, 5.4, 5.8, 6.1, 6.1, 6.1, 7, 7, 7, 7.3, 7.3, 7.5, 7.6, 7.7, 7.8, 7.9, 7.9, 8, 8, 8, 8, 8.1, 8.1, 8.1, 8.4, 8.5, 8.9, 9, 9.2, 10, 10.8, 12, 45]

def main():


    samples = [
        "Ust",
        "Los",
        "LBK",
        "LaB",
        "M45",
        "I15",
        "SF1",
        "I0708",
        "I1053",
        "I1116",
        "I1496",
        "I1507",
        "I1631",
        "I1632",
        "I1633",
        "I1734",
        "I1960",
        "I1961",
        "I1962",
        "I1965",
        "I2105",
        "I2134",
        "I2139",
        "I2514",
        "I2520",
        "I2861",
        "I2862",
        "I2935",
        "I2978",
        "I2980",
        "I3133",
        "I3255",
        "I3565",
        "I3758",
        "I3957",
        "I4315",
        "I4432",
        "I4438",
        "I4529",
        "I4532",
        "I4596",
        "I4634",
        "I4775",
        "I4792",
        "I4873",
        "I4875",
        "I4877",
        "I4878",
        "I4914",
        "I5077",
        "I5233",
        "I5235",
        "I5236",
        "I5319",
        "I5725",
        "I5742",
        "I5743",
        "I5748",
        "I5835",
        "I5838",
        "I5950",
        "I7421"]

    years = [
        45,
        8,
        7,
        8,
        12,
        7.5,
        9,

        8.1,
        3.8,
        1,
        7,
        7.7,
        6.1,
        6.1,
        4.5,
        9.2,
        8.1,
        6.1,
        7.3,
        5.8,
        4.9,
        7.6,
        8.9,
        4.3,
        5.1,
        2.8,
        1.6,
        5.1,
        5.1,
        5.2,
        5.4,
        4,
        3.1,
        2.2,
        3,
        3.5,
        7.9,
        7.3,
        1.8,
        1.2,
        8,
        5.4,
        3.4,
        3.3,
        7.9,
        8.4,
        8.5,
        7.8,
        8.1,
        7,
        8,
        10.8,
        10,
        1,
        2.5,
        4.7,
        4.7,
        4.4,
        4.2,
        4.6,
        4.5,
        3.8,
    ]

    print(len(years))
    print(len(samples))

    tuple_list = []

    # create a tuple list (sample_label, years)
    # then sort them
    for i in range(len(years)):
        tuple_list.append((samples[i], years[i]))

    sorted_tuple_list = sorted(tuple_list, key=lambda x: x[1])
    print(sorted_tuple_list)



    # technical matlab stuff
    pref = "load(\'C:\\Users\\amitai\\Google Drive\\p_CpGs\\f_"
    suf = ".mat\');"


    for elem in sorted_tuple_list:
        print(pref + elem[0][1:] + suf)

    for i, elem in enumerate(sorted_tuple_list):
        print(pref + elem[0][1:] + suf)
        print("samp" + str(i) + " = fas.getmethylation(1);")

    for i in range(len(sorted_tuple_list)):
        print("samp" + str(i) + " ", end="")


    # print all the labels
    for i in range(len(sorted_tuple_list)):
        print("\'"+ sorted_tuple_list[i][0] + "\'" + "," + " ", end="")


    # print all the years
    for i in range(len(sorted_tuple_list)):
        print(str(sorted_tuple_list[i][1])  + "," + " ", end="")


main()
