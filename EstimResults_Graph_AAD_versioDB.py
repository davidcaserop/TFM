import json 
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import numpy as np
import gSAFTmm as gSAFT
plt.style.use('classic')

#EXCELtoPYTHON
file = 'MonovalentSaltsExperimentalData_20200413.xlsx'
salts = 'KBr','KCl','KF','KI','LiBr','LiCl','LiI','NaBr','NaCl','NaF','NaI','RbBr','RbCl','RbF','RbI','HBr','KOH'

# salt = 'KF'
# df1 = pd.ExcelFile(file).parse(salt)
# m_2 = df1["m2"]
# m_1 = df1["m_1"]
# m_1 = m_1[np.logical_not(np.isnan(m_1))]

# compounds = {}
salts_name = 'PotassiumBromide', 'PotassiumChloride','PotassiumFluoride','PotassiumIodide','LithiumBromide','LithiumChloride','LithiumIodide','SodiumBromide','SodiumChloride','SodiumFluoride','SodiumIodide','RubidiumBromide','RubidiumChloride','RubidiumFluoride','RubidiumIodide','HydrogenBromide','PotassiumHydroxide'
salts_compounds = 'Potassium','Bromide','Potassium','Chlorine','Potassium','Fluoride','Potassium','Iodine','Lithium','Bromide','Lithium','Chlorine','Lithium','Iodine','Sodium','Bromide','Sodium','Chlorine','Sodium','Fluoride','Sodium','Iodine','Rubidium','Bromide','Rubidium','Chlorine','Rubidium','Fluoride','Rubidium','Iodine','Hydronium','Bromide','Potassium','Hydroxyl'



AAD_dl = []
AAD_bp = []
AAD_oc = []


results_AAD = open("AAD.txt","w+")



for i in range(0,len(salts)):
    

    df1 = pd.ExcelFile(file).parse(salts[i])
    m_2 = df1["m2"]
    m_2 = m_2[np.logical_not(np.isnan(m_2))]
    m_1 = df1["m_1"]
    m_1 = m_1[np.logical_not(np.isnan(m_1))]
    
    x_A = df1["X_A"]
    x_A = x_A[np.logical_not(np.isnan(x_A))]           
    x_B = df1["X_B"]
    x_B = x_B[np.logical_not(np.isnan(x_B))] 
    x_C = df1["X_C"]
    x_C = x_C[np.logical_not(np.isnan(x_C))]
    
    x_a = df1["X_a"]
    x_a = x_a[np.logical_not(np.isnan(x_a))]       
    x_b = df1["X_b"]
    x_b = x_b[np.logical_not(np.isnan(x_b))] 
    x_c = df1["X_c"]
    x_c = x_c[np.logical_not(np.isnan(x_c))]

    D = df1["D"]
    D = D[np.logical_not(np.isnan(D))]
    # print(len(x_a))
    # print(len(x_b))
    # print(len(x_c))
    # print(len(D))
    
    P_2 = df1["P_2"]
    P_2 = P_2[np.logical_not(np.isnan(P_2))]
    
    j = i
    a = 'new_'
    b = salts_name[i]
    
    c = '_ED_calculated.json'
    e = '_ED.json'
    d = a+b+c
    f = a+b+e
    
    cat = salts_compounds[i*2]
    an = salts_compounds[i*2+1]
    wat = 'GC_Mie_Databank_TL_18022020_NSW_17MX.xml<water,'
    #_NSW
    end = '>'
    gc = wat + cat +','+ an + end
    
    NSW = gSAFT.System(gc)
    # print(NSW)
    # with open(d) as json_data:
        # data = json.load(json_data)
        
        # tlist = list(zip(*data['BubblePressure']))
        # tlist_oc = list(zip(*data['OsmoticCoefficient']))
        # tlist_ld = list(zip(*data['SinglePhaseDensity']))
        
        
        # T = tlist[0]
        # VP_calculated = tlist[2]

        
        # oc_calculated = tlist_oc[3]
        # m_oc = tlist_oc[2]
        
        # x_calculated = tlist_ld[2]
        # x1_calculated =x_calculated[0]
        # T_LD_calculated = tlist_ld[0]
        # LD_calculated = tlist_ld[3] 

        
    with open(f) as json_data:
        Expdata = json.load(json_data)

        BP = Expdata['BubblePressure']['Data']
        OC = Expdata['OsmoticCoefficient']['Data']
        LD = Expdata['SinglePhaseDensity']['Data']
        
        tlist =list(zip(*BP))
        tlist_OC =list(zip(*OC))
        tlist_LD =list(zip(*LD))
        
        T = tlist[0]
        VP = tlist[2] 
        
        m_OC = tlist_OC[2]
        OC = tlist_OC[3]
        
        T_LD = tlist_LD[0]
        LD = tlist_LD[3]
    
    
    #comented from here:        
    #VP CALCULATED GSAFT    
    VP_calculated = []
    l_bp = len(VP)
    for i in range(0,l_bp):
    
        z = [ x_A[i] ,x_B[i] , x_C[i]]
        VP_calculated.append( NSW.BubblePressure(T[i],z) )
        
    #LD CALCULATED GSAFT     
    LD_calculated = []
    l_dl = len(D)
    for i in range(0,l_dl):
        z = [ x_a[i] ,x_b[i] , x_c[i]]
        LD_calculated.append( NSW.SinglePhaseDensity(T_LD[i],P_2[i],z) )
        
    #oc CALCULATED GSAFT 
    oc_calculated = []    
    l_oc = len(OC) 
    st = {salts_compounds[j*2]:1,salts_compounds[j*2+1]:1}
    for i in range(0,l_oc):
        T = 298.15
        P = 101325
        
        oc_calculated.append( NSW.OsmoticCoefficient(T,P,m_OC[i],Stoichiometry=st) )
        
        
    # AAD: 
    # Density
    
    t_dl = []
    for i in range(0,l_dl):
        t_dl.append(abs(D[i] - LD_calculated[i])/LD[i])
    AAD_dl.append(sum(t_dl)*100/l_dl)
    # # print('% AAD LD =',AAD_dl)   

    # BubblePressure
    
    t_bp = []
    for i in range(0,l_bp):
        t_bp.append(abs(VP[i] - VP_calculated[i])/VP[i])
    AAD_bp.append(sum(t_bp)*100/l_bp)
    # print('% AAD BP=',AAD_bp)   

    # OsmoticCoefficient
    
    t_oc = []
    for i in range(0,l_oc):
        t_oc.append(abs(OC[i] - oc_calculated[i])/OC[i])
    AAD_oc.append(sum(t_oc)*100/l_oc)
    # print('% AAD OC=',AAD_oc) 
    

    results_AAD.write("%s\t%f\t%f\t%f\n" %(salts[j],AAD_bp[j],AAD_dl[j],AAD_oc[j]))
    # results_AAD.write("%s\t%f\t%f\n" %(salts[j],AAD_bp[j],AAD_dl[j]))

    
########################################################################## 


# m_2 = []
# LD_calculated = []
# LD = []
#comented from here:
for i in range(0,len(salts)):

    # df1 = pd.ExcelFile(file).parse(salts[i])

    # m_2 = df1["m2"]
    # m_2 = m_2[np.logical_not(np.isnan(m_2))]
    # m_1 = df1["m_1"]
    # m_1 = m_1[np.logical_not(np.isnan(m_1))]
    
    # j = i
    # a = 'new_'
    # b = salts_name[i]
    # c = '_ED_calculated.json'
    # e = '_ED.json'
    # d = a+b+c
    # f = a+b+e

    # with open(d) as json_data:
        # data = json.load(json_data)
        
        # tlist = list(zip(*data['BubblePressure']))
        # tlist_oc = list(zip(*data['OsmoticCoefficient']))
        # tlist_ld = list(zip(*data['SinglePhaseDensity']))
        
        
        # T = tlist[0]
        # VP_calculated = tlist[2]
        
        # oc_calculated = tlist_oc[3]
        # m_oc = tlist_oc[2]
        
        # x_calculated = tlist_ld[2]
        # x1_calculated =x_calculated[0]
        # T_LD_calculated = tlist_ld[0]
        # LD_calculated.append( tlist_ld[3] )
        # LD_calculated = tlist_ld[3] 

        
    # with open(f) as json_data:
        # Expdata = json.load(json_data)

        # BP = Expdata['BubblePressure']['Data']
        # OC = Expdata['OsmoticCoefficient']['Data']
        # LD = Expdata['SinglePhaseDensity']['Data']
        
        # tlist =list(zip(*BP))
        # tlist_OC =list(zip(*OC))
        # tlist_LD =list(zip(*LD))
        
        # T = tlist[0]
        # VP = tlist[2] 
        
        # m_OC = tlist_OC[2]
        # OC = tlist_OC[3]
        
        # T_LD = tlist_LD[0]
        # LD.append(tlist_LD[3])
        # LD = tlist_LD[3]   
    df1 = pd.ExcelFile(file).parse(salts[i])
    m_2 = df1["m2"]
    m_2 = m_2[np.logical_not(np.isnan(m_2))]
    m_1 = df1["m_1"]
    m_1 = m_1[np.logical_not(np.isnan(m_1))]
    
    x_A = df1["X_A"]
    x_A = x_A[np.logical_not(np.isnan(x_A))]           
    x_B = df1["X_B"]
    x_B = x_B[np.logical_not(np.isnan(x_B))] 
    x_C = df1["X_C"]
    x_C = x_C[np.logical_not(np.isnan(x_C))]
    
    x_a = df1["X_a"]
    x_a = x_a[np.logical_not(np.isnan(x_a))]       
    x_b = df1["X_b"]
    x_b = x_b[np.logical_not(np.isnan(x_b))] 
    x_c = df1["X_c"]
    x_c = x_c[np.logical_not(np.isnan(x_c))]

    D = df1["D"]
    D = D[np.logical_not(np.isnan(D))]
    # print(len(x_a))
    # print(len(x_b))
    # print(len(x_c))
    # print(len(D))
    
    P_2 = df1["P_2"]
    P_2 = P_2[np.logical_not(np.isnan(P_2))]
    
    j = i
    a = 'new_'
    b = salts_name[i]
    
    c = '_ED_calculated.json'
    e = '_ED.json'
    d = a+b+c
    f = a+b+e
    
    cat = salts_compounds[i*2]
    an = salts_compounds[i*2+1]
    wat = 'GC_Mie_Databank_TL_18022020_NSW_17MX.xml<water,'
    #_NSW
    end = '>'
    gc = wat + cat +','+ an + end
    
    NSW = gSAFT.System(gc)
    # print(NSW)
    # with open(d) as json_data:
        # data = json.load(json_data)
        
        # tlist = list(zip(*data['BubblePressure']))
        # tlist_oc = list(zip(*data['OsmoticCoefficient']))
        # tlist_ld = list(zip(*data['SinglePhaseDensity']))
        
        
        # T = tlist[0]
        # VP_calculated = tlist[2]

        
        # oc_calculated = tlist_oc[3]
        # m_oc = tlist_oc[2]
        
        # x_calculated = tlist_ld[2]
        # x1_calculated =x_calculated[0]
        # T_LD_calculated = tlist_ld[0]
        # LD_calculated = tlist_ld[3] 

        
    with open(f) as json_data:
        Expdata = json.load(json_data)

        BP = Expdata['BubblePressure']['Data']
        OC = Expdata['OsmoticCoefficient']['Data']
        LD = Expdata['SinglePhaseDensity']['Data']
        
        tlist =list(zip(*BP))
        tlist_OC =list(zip(*OC))
        tlist_LD =list(zip(*LD))
        
        T = tlist[0]
        VP = tlist[2] 
        
        m_OC = tlist_OC[2]
        OC = tlist_OC[3]
        
        T_LD = tlist_LD[0]
        LD = tlist_LD[3]
    
    
    #comented from here:        
    #VP CALCULATED GSAFT    
    VP_calculated = []
    l_bp = len(VP)
    for i in range(0,l_bp):
    
        z = [ x_A[i] ,x_B[i] , x_C[i]]
        VP_calculated.append( NSW.BubblePressure(T[i],z) )
        
    #LD CALCULATED GSAFT     
    LD_calculated = []
    l_dl = len(D)
    for i in range(0,l_dl):
        z = [ x_a[i] ,x_b[i] , x_c[i]]
        LD_calculated.append( NSW.SinglePhaseDensity(T_LD[i],P_2[i],z) )
        
    #oc CALCULATED GSAFT 
    oc_calculated = []    
    l_oc = len(OC) 
    st = {salts_compounds[j*2]:1,salts_compounds[j*2+1]:1}
    for i in range(0,l_oc):
        T = 298.15
        P = 101325
        
        oc_calculated.append( NSW.OsmoticCoefficient(T,P,m_OC[i],Stoichiometry=st) )        

        
    # Graficos

    
    fig = plt.figure()
    plt.subplot(221)
    plt.plot(m_2,LD_calculated,'b',label="calculated")
    plt.plot(m_2,LD,'o',mfc='none',label="ExpData")
    # plt.legend(loc="upper left",scatterpoints=1,)
    plt.ylabel('Calculated')
    plt.ylabel("LD")
    plt.xlabel("m")
    fig.suptitle(salts_name[j])
    # plt.show()  

    plt.subplot(224)
    plt.plot(m_OC,oc_calculated,'b',label="calculated")
    plt.plot(m_OC,OC,'o',mfc='none',label="ExpData")
    # plt.legend(loc="upper left")
    plt.ylabel('Calculated')
    plt.ylabel("osm")
    plt.xlabel("m")
      

    plt.subplot(223)
    plt.plot(m_1,VP_calculated,'b',label="calculated") 
    plt.plot(m_1,VP,'o',mfc='none',label="ExpData")
    # plt.legend(loc="upper right")
    plt.ylabel('Calculated')
    plt.ylabel("P / Pa")
    plt.xlabel("m / molkg-1")

    fig.suptitle(salts_name[j])
    plt.show()      




