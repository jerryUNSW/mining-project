import pandas as pd
import numpy as np
import math

# some constants 
DENSITY_CONST = 0.7

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
                ## instead of getting the adjacent, get the next not-null value. 
                curr = i+1
                while pd.isnull(P0N_final[curr]):
                    curr = curr+1
                    if curr == len(P0N_final):
                        break
                if curr < len(P0N_final):
                    P01_EQ_list.append(P0N_final[curr])
                else:
                     P01_EQ_list.append(float("nan"))
            else:
                P01_EQ_list.append(0)
    return P01_EQ_list

def compute_R01(P02, P01, Q2, Q1):

    # print(P02, P01, Q2, Q1)

    # print("Delta = {}".format(pow( Q1*P01+(2*Q1+Q2)*P02 ,2) - 8*Q1*(Q1+Q2)*P01*P02))

    res=((Q1*P01+(2*Q1+Q2)*P02)- math.sqrt( pow( Q1*P01+(2*Q1+Q2)*P02 ,2) - 8*Q1*(Q1+Q2)*P01*P02 ) )/ (2*Q1*(Q1+Q2))

    return res


def compute_RN_from_R01_EQ(P02_EQ,P01_EQ, R01_EQ,Q01_EQ):
    res = (P02_EQ-P01_EQ)/(P01_EQ/R01_EQ-Q01_EQ)+R01_EQ
    return res 

# PONinterval, the PON array with NULL dropped. 
# not_null_index,   the indicex in the original fixed length array.
def compute_resistance(PONinterval, not_null_index, alr_list, P01_EQ_list, q1_eq_list, P02_EQ_list, q2_eq_list):

    num_of_holes = len(PONinterval)

    if num_of_holes == 0:
        return []

    if num_of_holes == 1:
        # when n == 1, R01 = P01 / Q1
        P01 = PONinterval[0]
        Q1 =  alr_list[not_null_index[0]]
        return [P01 / Q1]

    if num_of_holes == 2:
        P02 = PONinterval[0]
        P01 = PONinterval[1]
        Q2 = alr_list[not_null_index[0]]
        Q1 = alr_list[not_null_index[1]]
        R01 = 0 
        R02 = 0
        # print("P: {} {} Q: {} {}".format( P02, P01, Q2, Q1 ))
        if P02 - P01 >0:
            R01 = compute_R01(P02, P01, Q2, Q1)
            R02= (P02-P01)/(P01/R01-Q1) + R01
        else:
            R01 = P01 / Q1
            R02 = P02 / Q2 
        return [R02, R01]

    if num_of_holes == 3:
        P03 = PONinterval[0]
        P02 = PONinterval[1]
        P01 = PONinterval[2]
        Q3 = alr_list[not_null_index[0]]
        Q2 = alr_list[not_null_index[1]]
        Q1 = alr_list[not_null_index[2]]

        if P03 - P02 <= 0:
            # print("case 1")
            # in this case R03 is handled separately, and [R02, R02] is computed as in n == 2 (recurse back)
            R03 = P03 / Q3
            PONinterval_cp  = PONinterval.copy()
            PONinterval_cp.pop(0)
            not_null_index_cp = not_null_index.copy()
            not_null_index_cp.pop(0)
            # to compute R01 and R02, make a recursive call to itselfneed to remove the first element from PONinterval and not_null_index
            res = compute_resistance(PONinterval_cp, not_null_index_cp, alr_list, P01_EQ_list, q1_eq_list, P02_EQ_list, q2_eq_list)
            res.insert(0, R03)
            return res

        if P03 - P02 > 0 and P02 - P01 >0 :
            pos = not_null_index[0]
            # print("case 2", pos)
            P01_EQ = P01_EQ_list[pos]
            P02_EQ = P02_EQ_list[pos]

            Q01_EQ = 0
            for i in range(len(not_null_index)):
                if i>0:
                    Q01_EQ = Q01_EQ+alr_list[not_null_index[i]]
    
            Q02_EQ = q2_eq_list[pos]

            R01_EQ = compute_R01(P02_EQ, P01_EQ, Q02_EQ, Q01_EQ)

            R03 = compute_RN_from_R01_EQ(P02_EQ,P01_EQ, R01_EQ,Q01_EQ)
            R02 = 2 * R01_EQ 
            R01= ( P02+Q1*R02 -math.sqrt( pow(P02+Q1* R02 , 2) - 4*Q1*P01*R02) ) / (2*Q1) 

            # print("P01_EQ~Q02_EQ:" ,P01_EQ, P02_EQ, Q01_EQ, Q02_EQ)
            # print("R01_EQ = {}".format(R01_EQ))
            return [R03, R02, R01]

        if P03 - P02 > 0 and P02 - P01 <=0 :
            # print("case 3")
            pos = not_null_index[0]
            P01_EQ = P01_EQ_list[pos]
            P02_EQ = P02_EQ_list[pos]

            # Q01_EQ = q1_eq_list[pos]
            Q01_EQ = 0
            for i in range(len(not_null_index)):
                if i>0:
                    Q01_EQ = Q01_EQ+alr_list[not_null_index[i]]

            Q02_EQ = q2_eq_list[pos]
            R01_EQ = compute_R01(P02_EQ, P01_EQ, Q02_EQ, Q01_EQ)
            # compute the three resistance
            R03 = (P02_EQ-P01_EQ)/(P01_EQ/R01_EQ-Q01_EQ)+R01_EQ
            R01 = P01 / Q1 
            R02 = 1/ ((1/R01_EQ)-1/R01)
            return [R03, R02, R01]

    if num_of_holes > 3:
        P_N   = PONinterval[0]
        P_N_1 = PONinterval[1]
        Q_N   = alr_list[not_null_index[0]]
        RN = -1
        # now compute RN:
        if P_N - P_N_1 >0 : 
            pos = not_null_index[0]
            P01_EQ = P01_EQ_list[pos]
            P02_EQ = P02_EQ_list[pos]
            # Q01_EQ = q1_eq_list[pos]
            Q01_EQ = 0
            for i in range(len(not_null_index)):
                if i>0:
                    Q01_EQ = Q01_EQ+alr_list[not_null_index[i]]

            Q02_EQ = q2_eq_list[pos]
            R01_EQ = compute_R01(P02_EQ, P01_EQ, Q02_EQ, Q01_EQ)
            # need to compute RO1_EQ
            RN = compute_RN_from_R01_EQ(P02_EQ,P01_EQ, R01_EQ,Q01_EQ)
        else:
            RN = P_N / Q_N

        # make recursive call to solve the n-1 sub problem 
        PONinterval_cp  = PONinterval.copy()
        PONinterval_cp.pop(0)
        not_null_index_cp = not_null_index.copy()
        not_null_index_cp.pop(0)
        # to compute R01 and R02, make a recursive call to itselfneed to remove the first element from PONinterval and not_null_index
        res = compute_resistance(PONinterval_cp, not_null_index_cp, alr_list, P01_EQ_list, q1_eq_list, P02_EQ_list, q2_eq_list)
        res.insert(0, RN)
        return res

file_name = 'raw-data.xlsx'

df1 = pd.read_excel(file_name, sheet_name="Air leakage rate")

df2 = pd.read_excel(file_name, sheet_name="Suction pressure")

df3 = pd.read_excel(file_name, sheet_name="Density")

# print(len(df1), len(df2), len(df3))
if len(df1)!=len(df2) | len(df1)!=len(df3) | len(df2)!=len(df3):
    print("In consistent number of rows for tables. ")


## make this sheet a list of lists. 
resistance  = [ ]

for i in range(len(df1)):
    
    num_itr = i+2

    # if num_itr !=100:
    #     continue

    # clean the air leakage rate 
    alr_list  = df1.iloc[i].tolist()[1:]
    alr_list_dropna  = df1.iloc[i].dropna().tolist()[1:]
    # clean the suction pressure
    suc_press_list = df2.iloc[i].tolist()[1:]
    density_list   = df3.iloc[i].tolist()[1:]

    # special handling. 
    for i in range(len(density_list)):
        density_list[i] = DENSITY_CONST

    # n = len(alr_list)
    # compute Q1_EQ 
    q1_eq_list = compute_Q1_EQ(alr_list, alr_list_dropna)

    # compute Q2_EQ
    q2_eq_list = alr_list.copy()

    # compute P0N-final
    P0N_final = compute_P0N(suc_press_list, density_list)

    not_null_indices = []
    for i in range(len(P0N_final)):
        if not pd.isna(P0N_final[i]):
            not_null_indices.append(i)

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
                curr = i+1
                while pd.isnull(suc_press_list[curr]):
                    curr = curr+1
                    if curr == len(suc_press_list):
                        break
                if curr < len(suc_press_list):
                    PN_1_N.append(suc_val-suc_press_list[curr])
                else:
                    PN_1_N.append(float("nan"))
                # PN_1_N.append(suc_val-suc_press_list[i+1])
    PN_1_N.append(suc_press_list[-1])


    result = P0N_final.copy()

    # print("row number: ", num_itr, P0N_final)

    ## list of decreasing intervals 
    ## list of the original positions of these intervals. 
    curr_PON = []
    curr_PON_idx = []
    for idx in not_null_indices:
        val =P0N_final[idx]
        if len(curr_PON)==0:
            curr_PON.append(val)
            curr_PON_idx.append(idx)
        else:
            if val < curr_PON[-1]:
                curr_PON.append(val)
                curr_PON_idx.append(idx)
            else:
                # the interval is ready. start computing the resistance for it. 
                res =  compute_resistance(curr_PON, curr_PON_idx, alr_list, P01_EQ_list, q1_eq_list, P02_EQ_list, q2_eq_list)
                # print("interval: {} -> {}".format(curr_PON, res))
                # transfer the resistance to container:
                for itr in range(len(curr_PON)):
                    idx__ = curr_PON_idx[itr]
                    val__ = res[itr]
                    result[idx__] = val__

                curr_PON=[val]
                curr_PON_idx=[idx]

    res = compute_resistance(curr_PON, curr_PON_idx, alr_list, P01_EQ_list, q1_eq_list, P02_EQ_list, q2_eq_list)
    # print("interval: {} -> {}".format(curr_PON, res))
    # transfer the resistance to container:
    for itr in range(len(curr_PON)):
        idx__ = curr_PON_idx[itr]
        val__ = res[itr]
        # print(idx__, val)
        result[idx__] = val__

    resistance.append(result)
    # print(result)
    # print("")

df = pd.DataFrame(resistance[0:], columns=df1.columns[1:])

df.to_excel("resistance.xlsx")  
# print(df)




