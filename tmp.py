from asyncio import run_coroutine_threadsafe
import pandas as pd
import numpy as np
import math

def compute_Q1_EQ(alr_list, alr_list_dropna):

    if len(alr_list_dropna) == 0 :
        return []

    if len(alr_list_dropna) == 1:
        # need to return 
        tmp = alr_list.copy()
        for i in range(len(tmp)):
            if not pd.isnull(tmp[i]):
                tmp[i] = 0
        return tmp

    yy = alr_list.copy()
    yy.reverse()
    result = []
    runningsum = 0
    for item in yy:
        result.append(runningsum)
        if pd.isna(item):
            runningsum+=0
        else:
            runningsum+=item
    result.reverse()   
    # bug: the C10 of Q01-EQ
    for i in range(len(result)):
        if pd.isna(alr_list[i]):
            result[i] = float("nan")
    return result

def compute_P0N(suc_press_list, density_list):
    P0N_list = []
    for i in range(len(suc_press_list)):
        suc_val = suc_press_list[i]
        den_val = density_list[i]
        # special handling
        if pd.isna(den_val):
            den_val = 0
        if pd.isna(suc_val):
            P0N_list.append(float("nan"))
        else:
            pon_val = suc_val - den_val*9.81*250*0.001
            if pon_val<=0:
                P0N_list.append(float("nan"))
            else:
                P0N_list.append(pon_val)
    return P0N_list

def compute_P1_EQ(suc_press_list, P0N_final):
    P01_EQ_list = []
    for i in range(len(suc_press_list)):
        suc_val = suc_press_list[i]
        if pd.isna(suc_val):
            P01_EQ_list.append(float("nan"))
        else:
            if i+1 < len(P0N_final):
                P01_EQ_list.append(P0N_final[i+1])
            else:
                P01_EQ_list.append(0)
    return P01_EQ_list

def compute_R01(P02, P01, Q2, Q1):
    ### need to check delta >= 0 
    res=((Q1*P01+(2*Q1+Q2)*P02)- math.sqrt(pow( Q1*P01+(2*Q1+Q2)*P02 ,2) - 8*Q1*(Q1+Q2)*P01*P02 ) )/ (2*Q1*(Q1+Q2))
    return res

def compute_RN_from_R01_EQ(P02_EQ,P01_EQ, R01_EQ,Q01_EQ):
    res = (P02_EQ-P01_EQ)/(P01_EQ/R01_EQ-Q01_EQ)+R01_EQ
    return res 


def compute_resistance(P0N_final, alr_list, alr_list_dropna, P01_EQ_list, q1_eq_list, P02_EQ_list, q2_eq_list):


    # print(len(P0N_final), len(P01_EQ_list))
    P0N_final_dropna = []
    for item in P0N_final:
        if not pd.isnull(item):
            P0N_final_dropna.append(item)
    num_of_holes = len(P0N_final_dropna)

    if num_of_holes == 1:
        # when n == 1, R01 = P01 / Q1
        res = P0N_final.copy()
        for i in range(len(P0N_final)):
            if not pd.isnull(P0N_final[i]):
                res[i] = P0N_final[i] / alr_list[i]
        return res

    if num_of_holes == 2:
        P02 = P0N_final_dropna[0]
        P01 = P0N_final_dropna[1]
        Q2 = alr_list_dropna[0]
        Q1 = alr_list_dropna[1]
        R01 = 0 
        R02 = 0
        # print("P: {} {} Q: {} {}".format( P02, P01, Q2, Q1 ))
        if P02 - P01 >0:
            R01 = compute_R01(P02, P01, Q2, Q1)
            R02= (P02-P01)/(P01/R01-Q1) + R01
        else:
            R01 = P01 / Q1
            R02 = P02 / Q2 
        res_list = [R01, R02]
        res = P0N_final.copy()
        for i in range(len(P0N_final)):
            if not pd.isnull(P0N_final[i]):
                res[i] = res_list.pop()
        return res

    if num_of_holes == 3:
        P03 = P0N_final_dropna[0]
        P02 = P0N_final_dropna[1]
        P01 = P0N_final_dropna[2]
        Q3 = alr_list_dropna[0]
        Q2 = alr_list_dropna[1]
        Q1 = alr_list_dropna[2]

        if P03 - P02 <= 0:
            # in this case R03 is handled separately, and [R02, R02] is computed as in n == 2 (recurse back)
            R03 = P03 / Q3
            P0N_final_cp  = P0N_final.copy()
            alr_list_cp = alr_list.copy()
            for i in range(len(P0N_final)):
                if not pd.isnull(P0N_final[i]):
                    P0N_final_cp[i] = float("nan") 
                    alr_list_cp[i] = float("nan") 
                    break
            alr_list_dropna_cp = alr_list_dropna.copy()
            alr_list_dropna_cp.pop(0)
            # to compute R01 and R02, make a recursive call to itself
            intermediate_result = compute_resistance(P0N_final_cp, alr_list_cp, alr_list_dropna_cp, P01_EQ_list, q1_eq_list, P02_EQ_list, q2_eq_list)
            for i in range(len(intermediate_result)):

                if pd.isnull(intermediate_result[i]) and not pd.isnull(intermediate_result[i+1]):
                    intermediate_result[i] = R03
                    break 
            return intermediate_result

        if P03 - P02 > 0 and P02 - P01 >0 :
            # find the indices of the not null values: 
            not_null_indices = []
            for i in range(len(P0N_final)):
                if not pd.isnull(P0N_final[i]):
                    not_null_indices.append(i)

            pos = not_null_indices[0]
            P01_EQ = P01_EQ_list[pos]
            P02_EQ = P02_EQ_list[pos]
            Q01_EQ = q1_eq_list[pos]
            Q02_EQ = q2_eq_list[pos]

            R01_EQ = compute_R01(P02_EQ, P01_EQ, Q02_EQ, Q01_EQ)

            # R03 = (P02_EQ-P01_EQ)/(P01_EQ/R01_EQ-Q01_EQ)+R01_EQ
            R03 = compute_RN_from_R01_EQ(P02_EQ,P01_EQ, R01_EQ,Q01_EQ)

            R02 = 2 * R01_EQ 
            R01= ( P02+Q1*R02 -math.sqrt( pow(P02+Q1* R02 , 2) - 4*Q1*P01*R02) ) / (2*Q1) 
            to_pop = [R01, R02, R03]

            # print(R01_EQ, R02)
            res = P0N_final.copy()
            for idx in not_null_indices:
                res[idx] = to_pop.pop()
            return res

        if P03 - P02 > 0 and P02 - P01 <=0 :
            # find the indices of the not null values: 
            not_null_indices = []
            for i in range(len(P0N_final)):
                if not pd.isnull(P0N_final[i]):
                    not_null_indices.append(i)
            pos = not_null_indices[0]
            P01_EQ = P01_EQ_list[pos]
            P02_EQ = P02_EQ_list[pos]
            Q01_EQ = q1_eq_list[pos]
            Q02_EQ = q2_eq_list[pos]
            R01_EQ = compute_R01(P02_EQ, P01_EQ, Q02_EQ, Q01_EQ)
            R03=(P02_EQ-P01_EQ)/(P01_EQ/R01_EQ-Q01_EQ)+R01_EQ
            R01 = P01 / Q1 
            R02 = 1/ ((1/R01_EQ)-1/R01)
            to_pop = [R01, R02, R03]
            res = P0N_final.copy()
            for idx in not_null_indices:
                res[idx] = to_pop.pop()
            return res 
    if num_of_holes > 3:
        print(" num_of_holes = {} ".format(num_of_holes))
        # check Pn and Pn-1:
        P_N   = P0N_final_dropna[0]
        P_N_1 = P0N_final_dropna[1]
        Q_N   = alr_list_dropna[0]
        RN = -1
        # now compute RN:
        not_null_indices = []
        for i in range(len(P0N_final)):
            if not pd.isnull(P0N_final[i]):
                not_null_indices.append(i)
        if P_N - P_N_1 >0 : 
            ## maybe this can be optimized to run faster 
            pos = not_null_indices[0]
            P01_EQ = P01_EQ_list[pos]
            P02_EQ = P02_EQ_list[pos]
            Q01_EQ = q1_eq_list[pos]
            Q02_EQ = q2_eq_list[pos]
            R01_EQ = compute_R01(P02_EQ, P01_EQ, Q02_EQ, Q01_EQ)
            # need to compute RO1_EQ
            RN = compute_RN_from_R01_EQ(P02_EQ,P01_EQ, R01_EQ,Q01_EQ)
        else:
            RN = P_N / Q_N

        # make a copy for recursive call
        P0N_final_cp  = P0N_final.copy()
        alr_list_cp = alr_list.copy()
        P0N_final_cp[not_null_indices[0]] = float("nan") 
        alr_list_cp[not_null_indices[0]] = float("nan") 
        alr_list_dropna_cp = alr_list_dropna.copy()
        alr_list_dropna_cp.pop(0)

        res = compute_resistance(P0N_final_cp, alr_list_cp, alr_list_dropna_cp, 
            P01_EQ_list, q1_eq_list, P02_EQ_list, q2_eq_list)

        for i in range(len(res)-1):
            if pd.isnull(res[i]) and not pd.isnull(res[i+1]):
                res[i] = RN
                break
        return res


file_name = 'raw-data.xlsx'

df1 = pd.read_excel(file_name, sheet_name="Air leakage rate")

df2 = pd.read_excel(file_name, sheet_name="Suction pressure")

df3 = pd.read_excel(file_name, sheet_name="Density")

# print(len(df1), len(df2), len(df3))
if len(df1)!=len(df2) | len(df1)!=len(df3) | len(df2)!=len(df3):
    print("In consistent number of rows for tables. ")


## make this sheet a list of lists. 
Air_leakage_rate = []
Q1_EQ = []
for i in range(len(df1)):
    

    num_itr = i+2

    if num_itr !=51:
        continue

    # clean the air leakage rate 
    alr_list  = df1.iloc[i].tolist()[1:]
    alr_list_dropna  = df1.iloc[i].dropna().tolist()[1:]
    # clean the suction pressure
    suc_press_list = df2.iloc[i].tolist()[1:]
    density_list   = df3.iloc[i].tolist()[1:]

    # n = len(alr_list)
    # compute Q1_EQ 
    q1_eq_list = compute_Q1_EQ(alr_list, alr_list_dropna)

    # compute Q2_EQ
    q2_eq_list = alr_list.copy()

    # compute P0N-final
    P0N_final = compute_P0N(suc_press_list, density_list)

    # compute P01-EQ
    P01_EQ_list = compute_P1_EQ(suc_press_list, P0N_final)
    # compute P02-EQ
    P02_EQ_list = P0N_final.copy()

    # compute Pn-1,n
    PN_1_N = []
    for i in range(len(suc_press_list)-1):
        suc_val = suc_press_list[i]
        if pd.isna(suc_val):
            PN_1_N.append(float("nan"))
        else:
            if pd.isna(suc_press_list[i+1]):
                PN_1_N.append(suc_val)
            else:
                PN_1_N.append(suc_val-suc_press_list[i+1])
    PN_1_N.append(suc_press_list[-1])

    # print(P0N_final)

    resistance_list = compute_resistance(P0N_final, alr_list, alr_list_dropna, P01_EQ_list, q1_eq_list, P02_EQ_list, q2_eq_list)


    # print(num_itr)
    # if num_itr<40:
    if resistance_list is not None:

        tmp = []
        for item in resistance_list:
            if not pd.isnull(item):
                tmp.append(item)

        len_2 = 0
        for i in range(len(P0N_final)):
            if not pd.isnull(P0N_final[i]):      
                 len_2 = len_2+1

        # tmp.reverse()

        if len(tmp)!=len_2:
            print(num_itr, len(tmp), len_2)
            print(resistance_list)
            print(P0N_final)
            print("")

        # res = P0N_final.copy()

        # for i in range(len(P0N_final)):
        #     if not pd.isnull(P0N_final[i]):

        #         res[i] = tmp.pop()
        #         print("res[i] = {}".format(res[i]))
                
        # if len(tmp) < len__:
        #     print(len(tmp), len__, num_itr, res)

        # print(len(tmp), len__)

        # print(num_itr, "->", res)







