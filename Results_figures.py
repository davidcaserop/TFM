import json 
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import numpy as np
import gSAFTmm as gSAFT
plt.style.use('classic')
plt.rcParams["font.family"] = "Times New Roman"



#NSW DATABANK
NaCl_NSW = gSAFT.System("GC_Mie_Databank_TL_18022020_NSW_17MX.xml<water,Sodium,Chlorine>")
LiI_NSW = gSAFT.System("GC_Mie_Databank_TL_18022020_NSW_17MX.xml<water,Lithium,Iodine>")
LiBr_NSW = gSAFT.System("GC_Mie_Databank_TL_18022020_NSW_17MX.xml<water,Lithium,Bromide>")
LiCl_NSW = gSAFT.System("GC_Mie_Databank_TL_18022020_NSW_17MX.xml<water,Lithium,Chlorine>")
KCl_NSW = gSAFT.System("GC_Mie_Databank_TL_18022020_NSW_17MX.xml<water,Potassium,Chlorine>")
RbBr_NSW = gSAFT.System("GC_Mie_Databank_TL_18022020_NSW_17MX.xml<water,Rubidium,Bromide>")

HI_NSW = gSAFT.System("GC_Mie_Databank_TL_18022020_NSW_17MX.xml<water,Hydronium,Iodine>")
HBr_NSW = gSAFT.System("GC_Mie_Databank_TL_18022020_NSW_17MX.xml<water,Hydronium,Bromide>")
KOH_NSW = gSAFT.System("GC_Mie_Databank_TL_18022020_NSW_17MX.xml<water,Potassium,Hydroxyl>")

NaI_NSW = gSAFT.System("GC_Mie_Databank_TL_18022020_NSW_17MX.xml<water,Sodium,Iodine>")
NaBr_NSW = gSAFT.System("GC_Mie_Databank_TL_18022020_NSW_17MX.xml<water,Sodium,Bromide>")
KF_NSW = gSAFT.System("GC_Mie_Databank_TL_18022020_NSW_17MX.xml<water,Potassium,Fluoride>")

#SW DATABANK
NaCl = gSAFT.System("GC_Mie_Databank_TL_18022020_17MX.xml<water,Sodium,Chlorine>")
LiI = gSAFT.System("GC_Mie_Databank_TL_18022020_17MX.xml<water,Lithium,Iodine>")
LiBr = gSAFT.System("GC_Mie_Databank_TL_18022020_17MX.xml<water,Lithium,Bromide>")
LiCl = gSAFT.System("GC_Mie_Databank_TL_18022020_17MX.xml<water,Lithium,Chlorine>")
KCl = gSAFT.System("GC_Mie_Databank_TL_18022020_17MX.xml<water,Potassium,Chlorine>")
RbBr = gSAFT.System("GC_Mie_Databank_TL_18022020_17MX.xml<water,Rubidium,Bromide>")

HI = gSAFT.System("GC_Mie_Databank_TL_18022020_17MX.xml<water,Hydronium,Iodine>")
HBr = gSAFT.System("GC_Mie_Databank_TL_18022020_17MX.xml<water,Hydronium,Bromide>")
KOH = gSAFT.System("GC_Mie_Databank_TL_18022020_17MX.xml<water,Potassium,Hydroxyl>")

#STOICHIOMETRY
st_NaCl = {"Sodium" : 1, "Chlorine" : 1 }
st_LiI = {"Lithium" : 1, "Iodine" : 1 }
st_LiBr = {"Lithium" : 1, "Bromide" : 1 }
st_LiCl = {"Lithium" : 1, "Chlorine" : 1 }
st_KCl = {"Potassium" : 1, "Chlorine" : 1 }
st_RbBr = {"Rubidium" : 1, "Bromide" : 1 }
st_HI = {"Hydronium" : 1, "Iodine" : 1 }
st_HBr = {"Hydronium" : 1, "Bromide" : 1 }
st_KOH = {"Potassium" : 1, "Hydroxyl" : 1 }
st_NaI = {"Sodium" : 1, "Iodine" : 1 }
st_NaBr = {"Sodium" : 1, "Bromide" : 1 }
st_KF = {"Potassium" : 1, "Fluoride" : 1 }



#EXCELtoPYTHON
file = 'Eriksen_ExpData.xlsx'
figure = 'fig_3'
df1 = pd.ExcelFile(file).parse(figure)

m_1 = df1["m_1"]
m_1 = m_1[np.logical_not(np.isnan(m_1))]
m_1_E = df1["m_1_E"]
m_1_E = m_1_E[np.logical_not(np.isnan(m_1_E))]
oc_1 = df1["oc_1"]
oc_1 = oc_1[np.logical_not(np.isnan(oc_1))]
oc_1_E = df1["oc_1_E"]
oc_1_E = oc_1_E[np.logical_not(np.isnan(oc_1_E))]

m_2 = df1["m_2"]
m_2 = m_2[np.logical_not(np.isnan(m_2))]
m_2_E = df1["m_2_E"]
m_2_E = m_2_E[np.logical_not(np.isnan(m_2_E))]
oc_2 = df1["oc_2"]
oc_2 = oc_2[np.logical_not(np.isnan(oc_2))]
oc_2_E = df1["oc_2_E"]
oc_2_E = oc_2_E[np.logical_not(np.isnan(oc_2_E))]

m_3 = df1["m_3"]
m_3 = m_3[np.logical_not(np.isnan(m_3))]
oc_3 = df1["oc_3"]
oc_3 = oc_3[np.logical_not(np.isnan(oc_3))]
m_3_E = df1["m_3_E"]
m_3_E = m_3_E[np.logical_not(np.isnan(m_3_E))]
oc_3_E = df1["oc_3_E"]
oc_3_E = oc_3_E[np.logical_not(np.isnan(oc_3_E))]

m_4 = df1["m_4"]
m_4 = m_4[np.logical_not(np.isnan(m_4))]
oc_4 = df1["oc_4"]
oc_4 = oc_4[np.logical_not(np.isnan(oc_4))]
m_4_E = df1["m_4_E"]
m_4_E = m_4_E[np.logical_not(np.isnan(m_4_E))]
oc_4_E = df1["oc_4_E"]
oc_4_E = oc_4_E[np.logical_not(np.isnan(oc_4_E))]

#Results from gSAFT
T = 298.15
T_323 = 323.15
T_353 = 353.15
T_373 = 373.15
T_288 = 288.15
T_318 = 318.15
T_333 = 333.15
T_310 = 310.15
T_348 = 348.15
T_373 = 373.15

P = 101000
oc_NaCl= []
oc_LiI= []
oc_KCl = []
oc_RbBr = []
oc_NaCl_NSW= []
oc_LiI_NSW= []
oc_KCl_NSW = []
oc_RbBr_NSW = []
oc_HI_NSW= []
oc_KOH_NSW = []
oc_HBr_NSW = []

oc_HI= []
oc_KOH = []
oc_HBr = []

d_LiI_NSW = []
d_LiI_NSW_T = []
d_LiCl_NSW = []
d_LiCl_NSW_T = []
d_LiBr_NSW = []
d_LiBr_NSW_T = []

d_LiI = []
d_LiI_T = []
d_LiCl = []
d_LiCl_T = []
d_LiBr = []
d_LiBr_T = []


pT298_g_N  = []
pT323_g_N  = []
pT353_g_N  = []
pT373_g_N  = []

pT298_g  = []
pT323_g  = []
pT353_g  = []
pT373_g  = []

miac_LiBr_NSW = []
miac_LiCl_NSW = []
miac_NaI_NSW = []
miac_NaBr_NSW = []
miac_NaCl_NSW = []
miac_KF_NSW = []


miac_HI_NSW = []
miac_HBr_NSW = []
miac_KOH_NSW = []

miac_HI = []
miac_HBr = []
miac_KOH = []

miac_T288_g_N = []
miac_T298_g_N = []
miac_T318_g_N = []
miac_T333_g_N = []

oc_T288_g_N = []
oc_T298_g_N = []
oc_T310_g_N = []
oc_T323_g_N = []
oc_T348_g_N = []
oc_T373_g_N = []

miac_T288_g = []
miac_T298_g = []
miac_T318_g = []
miac_T333_g = []

oc_T288_g = []
oc_T298_g = []
oc_T310_g = []
oc_T323_g = []
oc_T348_g = []
oc_T373_g = []


MWw = 18.01528/1000
# m_min = 0.0001
m_min = 0.01
m_max = 6
m_num = (m_max-m_min+2)/0.01
# m_num = (m_max-m_min+2)/0.0001
m_space = np.linspace(m_min,m_max,m_num)
for i in m_space:
    m=i
    xMX = 1*m/(2*m + (1/MWw))
    z = [1-2*xMX, xMX, xMX]
    
    oc_NaCl.append( NaCl.OsmoticCoefficient(T,P,m,Stoichiometry=st_NaCl))
    oc_NaCl_NSW.append( NaCl_NSW.OsmoticCoefficient(T,P,m,Stoichiometry=st_NaCl))

  
    oc_LiI.append( LiI.OsmoticCoefficient(T,P,m,Stoichiometry=st_LiI))
    oc_LiI_NSW.append( LiI_NSW.OsmoticCoefficient(T,P,m,Stoichiometry=st_LiI))

   
    oc_KCl.append( KCl.OsmoticCoefficient(T,P,m,Stoichiometry=st_KCl))
    oc_KCl_NSW.append( KCl_NSW.OsmoticCoefficient(T,P,m,Stoichiometry=st_KCl))

   
    oc_RbBr.append( RbBr.OsmoticCoefficient(T,P,m,Stoichiometry=st_RbBr))
    oc_RbBr_NSW.append( RbBr_NSW.OsmoticCoefficient(T,P,m,Stoichiometry=st_RbBr))
    
    oc_HI_NSW.append( HI_NSW.OsmoticCoefficient(T,P,m,Stoichiometry=st_HI))
    
    oc_HBr_NSW.append( HBr_NSW.OsmoticCoefficient(T,P,m,Stoichiometry=st_HBr))
    
    oc_KOH_NSW.append( KOH_NSW.OsmoticCoefficient(T,P,m,Stoichiometry=st_KOH)) 
    

    oc_HI.append( HI.OsmoticCoefficient(T,P,m,Stoichiometry=st_HI))
    
    oc_HBr.append( HBr.OsmoticCoefficient(T,P,m,Stoichiometry=st_HBr))
    
    oc_KOH.append( KOH.OsmoticCoefficient(T,P,m,Stoichiometry=st_KOH))


    d_LiI_NSW.append( LiI_NSW.SinglePhaseDensity(T,P,z)/1000)
    d_LiI_NSW_T.append( LiI_NSW.SinglePhaseDensity(T_323,P,z)/1000) 
    d_LiCl_NSW.append( LiCl_NSW.SinglePhaseDensity(T,P,z)/1000)
    d_LiCl_NSW_T.append( LiCl_NSW.SinglePhaseDensity(T_323,P,z)/1000) 
    d_LiBr_NSW.append( LiBr_NSW.SinglePhaseDensity(T,P,z)/1000)
    d_LiBr_NSW_T.append( LiBr_NSW.SinglePhaseDensity(T_323,P,z)/1000)

    d_LiI.append( LiI.SinglePhaseDensity(T,P,z)/1000)
    d_LiI_T.append( LiI.SinglePhaseDensity(T_323,P,z)/1000) 
    d_LiCl.append( LiCl.SinglePhaseDensity(T,P,z)/1000)
    d_LiCl_T.append( LiCl.SinglePhaseDensity(T_323,P,z)/1000) 
    d_LiBr.append( LiBr.SinglePhaseDensity(T,P,z)/1000)
    d_LiBr_T.append( LiBr.SinglePhaseDensity(T_323,P,z)/1000)
    
    
    pT298_g_N.append( NaCl_NSW.BubblePressure(T,z) / 1000 )
    pT323_g_N.append(  NaCl_NSW.BubblePressure(T_323,z)/ 1000 )
    pT353_g_N.append( NaCl_NSW.BubblePressure(T_353,z)/ 1000 )
    pT373_g_N.append(  NaCl_NSW.BubblePressure(T_373,z)/ 1000 )

    pT298_g.append( NaCl.BubblePressure(T,z) / 1000 )
    pT323_g.append(  NaCl.BubblePressure(T_323,z)/ 1000 )
    pT353_g.append( NaCl.BubblePressure(T_353,z)/ 1000 )
    pT373_g.append(  NaCl.BubblePressure(T_373,z)/ 1000 )

    
    miac_LiBr_NSW.append( LiBr_NSW.MeanActivity(T,P,m,Stoichiometry=st_LiBr))
    miac_LiCl_NSW.append( LiCl_NSW.MeanActivity(T,P,m,Stoichiometry=st_LiCl))
    miac_NaI_NSW.append( NaI_NSW.MeanActivity(T,P,m,Stoichiometry=st_NaI))
    miac_NaBr_NSW.append( NaBr_NSW.MeanActivity(T,P,m,Stoichiometry=st_NaBr))
    miac_NaCl_NSW.append( NaCl_NSW.MeanActivity(T,P,m,Stoichiometry=st_NaCl))
    miac_KF_NSW.append( KF_NSW.MeanActivity(T,P,m,Stoichiometry=st_KF))
    
    miac_HI_NSW.append( HI_NSW.MeanActivity(T,P,m,Stoichiometry=st_HI))
    miac_HBr_NSW.append( HBr_NSW.MeanActivity(T,P,m,Stoichiometry=st_HBr))
    miac_KOH_NSW.append( KOH_NSW.MeanActivity(T,P,m,Stoichiometry=st_KOH))

    miac_HI.append( HI.MeanActivity(T,P,m,Stoichiometry=st_HI))
    miac_HBr.append( HBr.MeanActivity(T,P,m,Stoichiometry=st_HBr))
    miac_KOH.append( KOH.MeanActivity(T,P,m,Stoichiometry=st_KOH))
    
    miac_T288_g_N.append( NaCl_NSW.MeanActivity(T_288,P,m,Stoichiometry=st_NaCl))
    miac_T298_g_N.append( NaCl_NSW.MeanActivity(T,P,m,Stoichiometry=st_NaCl))
    miac_T318_g_N.append( NaCl_NSW.MeanActivity(T_318,P,m,Stoichiometry=st_NaCl))
    miac_T333_g_N.append( NaCl_NSW.MeanActivity(T_333,P,m,Stoichiometry=st_NaCl))
    
    
    oc_T288_g_N.append( NaCl_NSW.OsmoticCoefficient(T_288,P,m,Stoichiometry=st_NaCl))
    oc_T298_g_N.append( NaCl_NSW.OsmoticCoefficient(T,P,m,Stoichiometry=st_NaCl))
    oc_T310_g_N.append( NaCl_NSW.OsmoticCoefficient(T_310,P,m,Stoichiometry=st_NaCl))
    oc_T323_g_N.append( NaCl_NSW.OsmoticCoefficient(T_323,P,m,Stoichiometry=st_NaCl))
    oc_T348_g_N.append( NaCl_NSW.OsmoticCoefficient(T_348,P,m,Stoichiometry=st_NaCl))
    oc_T373_g_N.append( NaCl_NSW.OsmoticCoefficient(T_373,P,m,Stoichiometry=st_NaCl))
    
    
    miac_T288_g.append( NaCl.MeanActivity(T_288,P,m,Stoichiometry=st_NaCl))
    miac_T298_g.append( NaCl.MeanActivity(T,P,m,Stoichiometry=st_NaCl))
    miac_T318_g.append( NaCl.MeanActivity(T_318,P,m,Stoichiometry=st_NaCl))
    miac_T333_g.append( NaCl.MeanActivity(T_333,P,m,Stoichiometry=st_NaCl))
    
    
    oc_T288_g.append( NaCl.OsmoticCoefficient(T_288,P,m,Stoichiometry=st_NaCl))
    oc_T298_g.append( NaCl.OsmoticCoefficient(T,P,m,Stoichiometry=st_NaCl))
    oc_T310_g.append( NaCl.OsmoticCoefficient(T_310,P,m,Stoichiometry=st_NaCl))
    oc_T323_g.append( NaCl.OsmoticCoefficient(T_323,P,m,Stoichiometry=st_NaCl))
    oc_T348_g.append( NaCl.OsmoticCoefficient(T_348,P,m,Stoichiometry=st_NaCl))
    oc_T373_g.append( NaCl.OsmoticCoefficient(T_373,P,m,Stoichiometry=st_NaCl))
    
    
    

#GRAPH
#fig_3

plt.figure()
plt.plot(m_1,oc_1,'s',color='none',markeredgecolor='blue')
LiI_g, = plt.plot(m_space,oc_LiI,'--b',label="LiI")
# LiI_g, = plt.plot(m_1_E,oc_1_E,'--b',label="LiI")
LiI_g_N, = plt.plot(m_space,oc_LiI_NSW,'b',label="LiI")

plt.plot(m_2,oc_2,'s',color='none',markeredgecolor='red')
NaCl_g, = plt.plot(m_space,oc_NaCl,'--r',label="NaCl")
# NaCl_g, = plt.plot(m_2_E,oc_2_E,'--r',label="NaCl")
NaCl_g_N, = plt.plot(m_space,oc_NaCl_NSW,'r',label="NaCl")

plt.plot(m_3,oc_3,'s',color='none',markeredgecolor='green')
KCl_g, = plt.plot(m_space,oc_KCl,'--g',label="KCl")
# KCl_g, = plt.plot(m_3_E,oc_3_E,'--g',label="KCl")
KCl_g_N, = plt.plot(m_space,oc_KCl_NSW,'g',label="KCl")

plt.plot(m_4,oc_4,'s',color='none',markeredgecolor='black')
RbBr_g, = plt.plot(m_space,oc_RbBr,'--k',label="RbBr")
# RbBr_g, = plt.plot(m_4_E,oc_4_E,'--k',label="RbBr")
RbBr_g_N, = plt.plot(m_space,oc_RbBr_NSW,'k',label="RbBr")

plt.xlabel(r'${\itm}$ / mol kg$^{-1}$',fontsize=20)
plt.ylabel(r'$\Phi$',fontsize=20)
plt.legend([LiI_g_N,NaCl_g_N,KCl_g_N,RbBr_g_N],['LiI','NaCl','KCl','RbBr'],loc="upper left",fontsize=20, frameon=False)
axes = plt.gca()
# axes.set_xlim([xmin,xmax])
axes.set_ylim([0.8,1.8])
size = 15
plt.xticks(fontsize=size)
plt.yticks(fontsize=size)
plt.tight_layout()
# plt.figure()
plt.show()   



#EXCELtoPYTHON
figure = 'fig_5'
df1 = pd.ExcelFile(file).parse(figure)

m_1 = df1["m_1"]
m_1 = m_1[np.logical_not(np.isnan(m_1))]
m_1_E = df1["m_1_E"]
m_1_E = m_1_E[np.logical_not(np.isnan(m_1_E))]
oc_1 = df1["oc_1"]
oc_1 = oc_1[np.logical_not(np.isnan(oc_1))]
oc_1_E = df1["oc_1_E"]
oc_1_E = oc_1_E[np.logical_not(np.isnan(oc_1_E))]

m_2 = df1["m_2"]
m_2 = m_2[np.logical_not(np.isnan(m_2))]
m_2_E = df1["m_2_E"]
m_2_E = m_2_E[np.logical_not(np.isnan(m_2_E))]
oc_2 = df1["oc_2"]
oc_2 = oc_2[np.logical_not(np.isnan(oc_2))]
oc_2_E = df1["oc_2_E"]
oc_2_E = oc_2_E[np.logical_not(np.isnan(oc_2_E))]

m_3 = df1["m_3"]
m_3 = m_3[np.logical_not(np.isnan(m_3))]
oc_3 = df1["oc_3"]
oc_3 = oc_3[np.logical_not(np.isnan(oc_3))]
m_3_E = df1["m_3_E"]
m_3_E = m_3_E[np.logical_not(np.isnan(m_3_E))]
oc_3_E = df1["oc_3_E"]
oc_3_E = oc_3_E[np.logical_not(np.isnan(oc_3_E))]

#GRAPH
#fig_5

plt.figure()
plt.plot(m_1,oc_1,'s',color='none',markeredgecolor='blue')
HI_g, = plt.plot(m_1_E,oc_1_E,'--b',label="HI")
HI_g_N, = plt.plot(m_space,oc_HI_NSW,'b',label="HI")

plt.plot(m_2,oc_2,'s',color='none',markeredgecolor='red')
HBr_g, = plt.plot(m_2_E,oc_2_E,'--r',label="HBr")
HBr_g_N, = plt.plot(m_space,oc_HBr_NSW,'r',label="HBr")

plt.plot(m_3,oc_3,'s',color='none',markeredgecolor='green')
KOH_g, = plt.plot(m_3_E,oc_3_E,'--g',label="KOH")
KOH_g_N, = plt.plot(m_space,oc_KOH_NSW,'g',label="KOH")

plt.xlabel(r'${\itm}$ / mol kg$^{-1}$',fontsize=20)
plt.ylabel(r'$\Phi$',fontsize=20)
plt.legend([HI_g_N,HBr_g_N,KOH_g_N],['HI','HBr','KOH'],loc="upper left",fontsize=20, frameon=False)
axes = plt.gca()
# axes.set_xlim([xmin,xmax])
axes.set_ylim([0.8,1.8])
size = 15
plt.xticks(fontsize=size)
plt.yticks(fontsize=size)
plt.tight_layout()
plt.show()   


#EXCELtoPYTHON
figure = 'fig_6'
df1 = pd.ExcelFile(file).parse(figure)

m_1 = df1["m_1"]
m_1 = m_1[np.logical_not(np.isnan(m_1))]
m_1_E = df1["m_1_E"]
m_1_E = m_1_E[np.logical_not(np.isnan(m_1_E))]
d_1 = df1["d_1"]
d_1 = d_1[np.logical_not(np.isnan(d_1))]
d_1_E = df1["d_1_E"]
d_1_E = d_1_E[np.logical_not(np.isnan(d_1_E))]
m_1_T = df1["m_1_T"]
m_1_T = m_1_T[np.logical_not(np.isnan(m_1_T))]
m_1_E_T = df1["m_1_E_T"]
m_1_E_T = m_1_E_T[np.logical_not(np.isnan(m_1_E_T))]
d_1_T = df1["d_1_T"]
d_1_T = d_1_T[np.logical_not(np.isnan(d_1_T))]
d_1_E_T= df1["d_1_E_T"]
d_1_E_T = d_1_E_T[np.logical_not(np.isnan(d_1_E_T))]

m_2 = df1["m_2"]
m_2 = m_2[np.logical_not(np.isnan(m_2))]
m_2_E = df1["m_2_E"]
m_2_E = m_2_E[np.logical_not(np.isnan(m_2_E))]
d_2 = df1["d_2"]
d_2 = d_2[np.logical_not(np.isnan(d_2))]
d_2_E = df1["d_2_E"]
d_2_E = d_2_E[np.logical_not(np.isnan(d_2_E))]

m_3 = df1["m_3"]
m_3 = m_3[np.logical_not(np.isnan(m_3))]
d_3 = df1["d_3"]
d_3 = d_3[np.logical_not(np.isnan(d_3))]
m_3_E = df1["m_3_E"]
m_3_E = m_3_E[np.logical_not(np.isnan(m_3_E))]
d_3_E = df1["d_3_E"]
d_3_E = d_3_E[np.logical_not(np.isnan(d_3_E))]


#GRAPH
#fig_6

plt.figure()
#LiI
#ExpData
plt.plot(m_1,d_1,'s',color='none',markeredgecolor='blue')
# plt.plot(m_1_T,d_1_T,'o',color='none',markeredgecolor='blue')
#Spherical 
# LiI_g, = plt.plot(m_1_E,d_1_E,'--b')
# LiI_g_2, = plt.plot(m_1_E,d_1_E,'sb')
# LiI_g_T, = plt.plot(m_1_E_T,d_1_E_T,':b')
# LiI_g_2_T, = plt.plot(m_1_E_T,d_1_E,'ob')
#nonSpherical
plt.plot(m_space,d_LiI,'--b',label="LiI")
LiI_g_N, = plt.plot(m_space,d_LiI_NSW,'b',label="LiI")
# LiI_g_N_2, = plt.plot(m_space,d_LiI_NSW,'bs',label="LiI")
# LiI_g_N_T, = plt.plot(m_space,d_LiI_NSW_T,'--b',label="LiI")
# LiI_g_N_2_T, = plt.plot(m_space,d_LiI_NSW_T,'bo',label="LiI")
#LiBr
#ExpData
plt.plot(m_2,d_2,'s',color='none',markeredgecolor='red')
#Spherical 
LiBr_g, = plt.plot(m_2_E,d_2_E,'--r',label="LiBr")
#nonSpherical
LiBr_g_N, = plt.plot(m_space,d_LiBr_NSW,'r',label="HBr")
#LiCl
#ExpData
plt.plot(m_3,d_3,'s',color='none',markeredgecolor='green')
LiCl_g, = plt.plot(m_3_E,d_3_E,'--g')
LiCl_g_N, = plt.plot(m_space,d_LiCl_NSW,'g')

plt.xlabel(r'${\itm}$ / mol kg$^{-1}$',fontsize=20)
plt.ylabel(r'${\rho}$ / g dm$^{-3}$ ',fontsize=20)
# plt.ylabel(r'$\Phi$')
# plt.ylabel(r'$\it{Rho}$/ $g dm^{-3}$')
plt.legend([LiI_g_N,LiBr_g_N,LiCl_g_N],['LiI','LiBr','LiCl'],loc="upper left",fontsize=20, frameon=False)
# plt.legend([LiI_g_N],['LiI'],loc="upper left")
axes = plt.gca()
# axes.set_xlim([xmin,xmax])
axes.set_ylim([0.9,1.4])
size = 15
plt.xticks(fontsize=size)
plt.yticks(fontsize=size)
plt.tight_layout()
plt.show()   


#EXCELtoPYTHON
figure = 'fig_7'
df1 = pd.ExcelFile(file).parse(figure)


m_1 = df1["m_1"]
m_1 = m_1[np.logical_not(np.isnan(m_1))]
m_1_E = df1["m_1_E"]
m_1_E = m_1_E[np.logical_not(np.isnan(m_1_E))]
p_1 = df1["p_1"]
p_1 = p_1[np.logical_not(np.isnan(p_1))]
p_1_E = df1["p_1_E"]
p_1_E = p_1_E[np.logical_not(np.isnan(p_1_E))]

m_2 = df1["m_2"]
m_2 = m_2[np.logical_not(np.isnan(m_2))]
m_2_E = df1["m_2_E"]
m_2_E = m_2_E[np.logical_not(np.isnan(m_2_E))]
p_2 = df1["p_2"]
p_2 = p_2[np.logical_not(np.isnan(p_2))]
p_2_E = df1["p_2_E"]
p_2_E = p_2_E[np.logical_not(np.isnan(p_2_E))]

m_3 = df1["m_3"]
m_3 = m_3[np.logical_not(np.isnan(m_3))]
p_3 = df1["p_3"]
p_3 = p_3[np.logical_not(np.isnan(p_3))]
m_3_E = df1["m_3_E"]
m_3_E = m_3_E[np.logical_not(np.isnan(m_3_E))]
p_3_E = df1["p_3_E"]
p_3_E = p_3_E[np.logical_not(np.isnan(p_3_E))]

m_4 = df1["m_4"]
m_4 = m_4[np.logical_not(np.isnan(m_4))]
p_4 = df1["p_4"]
p_4 = p_4[np.logical_not(np.isnan(p_4))]
m_4_E = df1["m_4_E"]
m_4_E = m_4_E[np.logical_not(np.isnan(m_4_E))]
p_4_E = df1["p_4_E"]
p_4_E = p_4_E[np.logical_not(np.isnan(p_4_E))]

#GRAPH
#fig_7

plt.figure()
plt.plot(m_1,p_1,'s',color='none',markeredgecolor='blue')
pT298_g, = plt.plot(m_1_E,p_1_E,'--b')
pT298_g_N, = plt.plot(m_space,pT298_g_N,'b')

plt.plot(m_2,p_2,'s',color='none',markeredgecolor='red')
pT323_g, = plt.plot(m_2_E,p_2_E,'--r')
pT323_g_N, = plt.plot(m_space,pT323_g_N,'r')

plt.plot(m_3,p_3,'s',color='none',markeredgecolor='green')
pT353_g, = plt.plot(m_3_E,p_3_E,'--g')
pT353_g_N, = plt.plot(m_space,pT353_g_N,'g')

plt.plot(m_4,p_4,'s',color='none',markeredgecolor='black')
pT373_g, = plt.plot(m_4_E,p_4_E,'--k')
pT373_g_N, = plt.plot(m_space,pT373_g_N,'k')


plt.xlabel(r'${\itm}$ / mol kg$^{-1}$',fontsize=20)
plt.ylabel(r'$\it{P}$ / kPa',fontsize=20)
# plt.legend([pT298_g_N,pT323_g_N,pT353_g_N,pT373_g_N],['298K','323K','353K','373K'],loc="upper right",fontsize=20, frameon=False)
plt.legend([pT298_g_N,pT323_g_N,pT353_g_N,pT373_g_N],['298K','323K','353K','373K'],loc="upper right",fontsize=20)
axes = plt.gca()
# axes.set_xlim([xmin,xmax])
axes.set_ylim([0,110])
size = 15
plt.xticks(fontsize=size)
plt.yticks(fontsize=size)
plt.tight_layout()
plt.show()   


#EXCELtoPYTHON
figure = 'fig_8'
df1 = pd.ExcelFile(file).parse(figure)


m_1 = df1["m_1"]
m_1 = m_1[np.logical_not(np.isnan(m_1))]
m_1_E = df1["m_1_E"]
m_1_E = m_1_E[np.logical_not(np.isnan(m_1_E))]
miac_1 = df1["miac_1"]
miac_1 = miac_1[np.logical_not(np.isnan(miac_1))]
miac_1_E = df1["miac_1_E"]
miac_1_E = miac_1_E[np.logical_not(np.isnan(miac_1_E))]

m_2 = df1["m_2"]
m_2 = m_2[np.logical_not(np.isnan(m_2))]
m_2_E = df1["m_2_E"]
m_2_E = m_2_E[np.logical_not(np.isnan(m_2_E))]
miac_2 = df1["miac_2"]
miac_2 = miac_2[np.logical_not(np.isnan(miac_2))]
miac_2_E = df1["miac_2_E"]
miac_2_E = miac_2_E[np.logical_not(np.isnan(miac_2_E))]

m_3 = df1["m_3"]
m_3 = m_3[np.logical_not(np.isnan(m_3))]
miac_3 = df1["miac_3"]
miac_3 = miac_3[np.logical_not(np.isnan(miac_3))]
m_3_E = df1["m_3_E"]
m_3_E = m_3_E[np.logical_not(np.isnan(m_3_E))]
miac_3_E = df1["miac_3_E"]
miac_3_E = miac_3_E[np.logical_not(np.isnan(miac_3_E))]

m_4 = df1["m_4"]
m_4 = m_4[np.logical_not(np.isnan(m_4))]
miac_4 = df1["miac_4"]
miac_4 = miac_4[np.logical_not(np.isnan(miac_4))]
m_4_E = df1["m_4_E"]
m_4_E = m_4_E[np.logical_not(np.isnan(m_4_E))]
miac_4_E = df1["miac_4_E"]
miac_4_E = miac_4_E[np.logical_not(np.isnan(miac_4_E))]

m_5 = df1["m_5"]
m_5 = m_5[np.logical_not(np.isnan(m_5))]
miac_5 = df1["miac_5"]
miac_5 = miac_5[np.logical_not(np.isnan(miac_5))]
m_5_E = df1["m_5_E"]
m_5_E = m_5_E[np.logical_not(np.isnan(m_5_E))]
miac_5_E = df1["miac_5_E"]
miac_5_E = miac_5_E[np.logical_not(np.isnan(miac_5_E))]

m_6 = df1["m_6"]
m_6 = m_6[np.logical_not(np.isnan(m_6))]
miac_6 = df1["miac_6"]
miac_6 = miac_6[np.logical_not(np.isnan(miac_6))]
m_6_E = df1["m_6_E"]
m_6_E = m_6_E[np.logical_not(np.isnan(m_6_E))]
miac_6_E = df1["miac_6_E"]
miac_6_E = miac_6_E[np.logical_not(np.isnan(miac_6_E))]

#GRAPH
#fig_8

plt.figure()
plt.plot(m_1,miac_1,'s',color='none',markeredgecolor='blue')
# LiBr_g, = plt.plot(m_1_E,miac_1_E,'--b')
LiBr_f8_N, = plt.plot(m_space,miac_LiBr_NSW)

plt.plot(m_2,miac_2,'s',color='none',markeredgecolor='red')
# LiCl_g, = plt.plot(m_2_E,miac_2_E,'--r')
LiCl_f8_N, = plt.plot(m_space,miac_LiCl_NSW,'r')

plt.plot(m_3,miac_3,'s',color='none',markeredgecolor='magenta')
# NaI_g, = plt.plot(m_3_E,miac_3_E,'--m')
NaI_f8_N, = plt.plot(m_space,miac_NaI_NSW,'m')

plt.plot(m_4,miac_4,'s',color='none',markeredgecolor='green')
# NaBr_g, = plt.plot(m_4_E,miac_4_E,'--g')
NaBr_f8_N, = plt.plot(m_space,miac_NaBr_NSW,'g')

plt.plot(m_5,miac_5,'s',color='none',markeredgecolor='black')
# NaCl_g, = plt.plot(m_5_E,miac_5_E,'--k',label="NaCl")
NaCl_f8_N, = plt.plot(m_space,miac_NaCl_NSW,'k',label="NaCl")

plt.plot(m_6,miac_6,'s',color='none',markeredgecolor='yellow')
# KF_g, = plt.plot(m_6_E,miac_6_E,'--y',label="KF")
KF_f8_N, = plt.plot(m_space,miac_KF_NSW,'y',label="KF")

plt.xlabel(r'${\itm}$ / mol kg$^{-1}$',fontsize=20)
# plt.ylabel(r'\it{$\gamma}$\pm,$ \it{m}')
plt.ylabel(r'$\it{\gamma}\pm , m$',fontsize=20)
plt.legend([LiBr_f8_N,LiCl_f8_N,NaI_f8_N,NaBr_f8_N,NaCl_f8_N,KF_f8_N],['LiBr','LiCl','NaI','NaBr','NaCl','KF'],loc="upper left",fontsize=20)
axes = plt.gca()
# axes.set_xlim([xmin,xmax])
axes.set_ylim([0,3])
size = 15
plt.xticks(fontsize=size)
plt.yticks(fontsize=size)
plt.tight_layout()
plt.show()   


#EXCELtoPYTHON
figure = 'fig_10'
df1 = pd.ExcelFile(file).parse(figure)

m_1 = df1["m_1"]
m_1 = m_1[np.logical_not(np.isnan(m_1))]
m_1_E = df1["m_1_E"]
m_1_E = m_1_E[np.logical_not(np.isnan(m_1_E))]
miac_1 = df1["miac_1"]
miac_1 = miac_1[np.logical_not(np.isnan(miac_1))]
miac_1_E = df1["miac_1_E"]
miac_1_E = miac_1_E[np.logical_not(np.isnan(miac_1_E))]

m_2 = df1["m_2"]
m_2 = m_2[np.logical_not(np.isnan(m_2))]
m_2_E = df1["m_2_E"]
m_2_E = m_2_E[np.logical_not(np.isnan(m_2_E))]
miac_2 = df1["miac_2"]
miac_2 = miac_2[np.logical_not(np.isnan(miac_2))]
miac_2_E = df1["miac_2_E"]
miac_2_E = miac_2_E[np.logical_not(np.isnan(miac_2_E))]

m_3 = df1["m_3"]
m_3 = m_3[np.logical_not(np.isnan(m_3))]
miac_3 = df1["miac_3"]
miac_3 = miac_3[np.logical_not(np.isnan(miac_3))]
m_3_E = df1["m_3_E"]
m_3_E = m_3_E[np.logical_not(np.isnan(m_3_E))]
miac_3_E = df1["miac_3_E"]
miac_3_E = miac_3_E[np.logical_not(np.isnan(miac_3_E))]

#GRAPH
#fig_10

plt.figure()
plt.plot(m_1,miac_1,'s',color='none',markeredgecolor='blue')
HI_g, = plt.plot(m_space,miac_HI,'--b')
HI_g_N, = plt.plot(m_space,miac_HI_NSW)

plt.plot(m_2,miac_2,'s',color='none',markeredgecolor='red')
HBr_g, = plt.plot(m_space,miac_HBr,'--r')
HBr_g_N, = plt.plot(m_space,miac_HBr_NSW,'r')

plt.plot(m_3,miac_3,'s',color='none',markeredgecolor='green')
KOH_g, = plt.plot(m_space,miac_KOH,'--g')
KOH_g_N, = plt.plot(m_space,miac_KOH_NSW,'g')

plt.xlabel(r'${\itm}$ / mol kg$^{-1}$',fontsize=20)
plt.ylabel(r'$\it{\gamma}\pm , m$',fontsize=20)
plt.legend([HI_g_N,HBr_g_N,KOH_g_N],['HI','HBr','KOH'],loc="upper left",fontsize=20,frameon=False)
axes = plt.gca()
# axes.set_xlim([xmin,xmax])
axes.set_ylim([0.5,4])
size = 15
plt.xticks(fontsize=size)
plt.yticks(fontsize=size)
plt.tight_layout()
plt.show()   

#EXCELtoPYTHON
figure = 'fig_11'
df1 = pd.ExcelFile(file).parse(figure)

m_1 = df1["m_1"]
m_1 = m_1[np.logical_not(np.isnan(m_1))]
m_1_E = df1["m_1_E"]
m_1_E = m_1_E[np.logical_not(np.isnan(m_1_E))]
p_1 = df1["p_1"]
p_1 = p_1[np.logical_not(np.isnan(p_1))]
p_1_E = df1["p_1_E"]
p_1_E = p_1_E[np.logical_not(np.isnan(p_1_E))]

m_2 = df1["m_2"]
m_2 = m_2[np.logical_not(np.isnan(m_2))]
m_2_E = df1["m_2_E"]
m_2_E = m_2_E[np.logical_not(np.isnan(m_2_E))]
p_2 = df1["p_2"]
p_2 = p_2[np.logical_not(np.isnan(p_2))]
p_2_E = df1["p_2_E"]
p_2_E = p_2_E[np.logical_not(np.isnan(p_2_E))]

m_3 = df1["m_3"]
m_3 = m_3[np.logical_not(np.isnan(m_3))]
p_3 = df1["p_3"]
p_3 = p_3[np.logical_not(np.isnan(p_3))]
m_3_E = df1["m_3_E"]
m_3_E = m_3_E[np.logical_not(np.isnan(m_3_E))]
p_3_E = df1["p_3_E"]
p_3_E = p_3_E[np.logical_not(np.isnan(p_3_E))]

m_4 = df1["m_4"]
m_4 = m_4[np.logical_not(np.isnan(m_4))]
p_4 = df1["p_4"]
p_4 = p_4[np.logical_not(np.isnan(p_4))]
m_4_E = df1["m_4_E"]
m_4_E = m_4_E[np.logical_not(np.isnan(m_4_E))]
p_4_E = df1["p_4_E"]
p_4_E = p_4_E[np.logical_not(np.isnan(p_4_E))]

#GRAPH
#fig_11_a&b

plt.figure()
plt.plot(m_1,p_1,'s',color='none',markeredgecolor='blue')
miac_T288_g, = plt.plot(m_space,miac_T288_g,'--b')
miac_T288_g_N, = plt.plot(m_space,miac_T288_g_N,'b')

plt.plot(m_2,p_2,'s',color='none',markeredgecolor='red')
# miac_T298_g_S, = plt.plot(m_2_E,p_2_E,'--r')
miac_T298_g, = plt.plot(m_space,miac_T298_g,'--r')
miac_T298_g_N, = plt.plot(m_space,miac_T298_g_N,'r')

plt.plot(m_3,p_3,'s',color='none',markeredgecolor='green')
# miac_T318_g_S, = plt.plot(m_3_E,p_3_E,'--g')
miac_T318_g, = plt.plot(m_space,miac_T318_g,'--g')
miac_T318_g_N, = plt.plot(m_space,miac_T318_g_N,'g')

plt.plot(m_4,p_4,'s',color='none',markeredgecolor='black')
# miac_T333_g_S, = plt.plot(m_4_E,p_4_E,'--k')
miac_T333_g_N, = plt.plot(m_space,miac_T333_g_N,'k')
miac_T333_g, = plt.plot(m_space,miac_T333_g,'--k')

plt.xlabel(r'${\itm}$ / mol kg$^{-1}$',fontsize=20)
plt.ylabel(r'$\it{\gamma}\pm , m$',fontsize=20)
plt.legend([miac_T288_g_N,miac_T298_g_N,miac_T318_g_N,miac_T333_g_N],['288K','298K','318K','333K'],loc="upper right",fontsize=20, frameon=False)
axes = plt.gca()
axes.set_xlim([0,0.1])
axes.set_ylim([0.75,1])
# axes.set_xlim([0,6])
# axes.set_ylim([0.5,1.1])
size = 15
plt.xticks(fontsize=size)
plt.yticks(fontsize=size)
plt.tight_layout()
plt.show() 

#EXCELtoPYTHON
figure = 'fig_12'
df1 = pd.ExcelFile(file).parse(figure)

m_1 = df1["m_1"]
m_1 = m_1[np.logical_not(np.isnan(m_1))]
m_1_E = df1["m_1_E"]
m_1_E = m_1_E[np.logical_not(np.isnan(m_1_E))]
p_1 = df1["p_1"]
p_1 = p_1[np.logical_not(np.isnan(p_1))]
p_1_E = df1["p_1_E"]
p_1_E = p_1_E[np.logical_not(np.isnan(p_1_E))]

m_2 = df1["m_2"]
m_2 = m_2[np.logical_not(np.isnan(m_2))]
m_2_E = df1["m_2_E"]
m_2_E = m_2_E[np.logical_not(np.isnan(m_2_E))]
p_2 = df1["p_2"]
p_2 = p_2[np.logical_not(np.isnan(p_2))]
p_2_E = df1["p_2_E"]
p_2_E = p_2_E[np.logical_not(np.isnan(p_2_E))]

m_3 = df1["m_3"]
m_3 = m_3[np.logical_not(np.isnan(m_3))]
p_3 = df1["p_3"]
p_3 = p_3[np.logical_not(np.isnan(p_3))]
m_3_E = df1["m_3_E"]
m_3_E = m_3_E[np.logical_not(np.isnan(m_3_E))]
p_3_E = df1["p_3_E"]
p_3_E = p_3_E[np.logical_not(np.isnan(p_3_E))]

m_4 = df1["m_4"]
m_4 = m_4[np.logical_not(np.isnan(m_4))]
p_4 = df1["p_4"]
p_4 = p_4[np.logical_not(np.isnan(p_4))]
m_4_E = df1["m_4_E"]
m_4_E = m_4_E[np.logical_not(np.isnan(m_4_E))]
p_4_E = df1["p_4_E"]
p_4_E = p_4_E[np.logical_not(np.isnan(p_4_E))]

m_5 = df1["m_5"]
m_5 = m_5[np.logical_not(np.isnan(m_5))]
p_5 = df1["p_5"]
p_5 = p_5[np.logical_not(np.isnan(p_5))]
m_5_E = df1["m_5_E"]
m_5_E = m_5_E[np.logical_not(np.isnan(m_5_E))]
p_5_E = df1["p_5_E"]
p_5_E = p_5_E[np.logical_not(np.isnan(p_5_E))]

m_6 = df1["m_6"]
m_6 = m_6[np.logical_not(np.isnan(m_6))]
p_6 = df1["p_6"]
p_6 = p_6[np.logical_not(np.isnan(p_6))]
m_6_E = df1["m_6_E"]
m_6_E = m_6_E[np.logical_not(np.isnan(m_6_E))]
p_6_E = df1["p_6_E"]
p_6_E = p_6_E[np.logical_not(np.isnan(p_6_E))]

#GRAPH
#fig_12

plt.figure()
plt.plot(m_1,p_1,'s',color='none',markeredgecolor='blue')
# oc_T288_g_S, = plt.plot(m_1_E,p_1_E,'--b')
oc_T288_g_N, = plt.plot(m_space,oc_T288_g_N,'b')
oc_T288_g, = plt.plot(m_space,oc_T288_g,'--b')

plt.plot(m_2,p_2,'s',color='none',markeredgecolor='red')
# oc_T298_g_S, = plt.plot(m_2_E,p_2_E,'--r')
oc_T298_g_N, = plt.plot(m_space,oc_T298_g_N,'r')
oc_T298_g, = plt.plot(m_space,oc_T298_g,'--r')

plt.plot(m_3,p_3,'s',color='none',markeredgecolor='green')
# oc_T310_g_S, = plt.plot(m_3_E,p_3_E,'--g')
oc_T310_g_N, = plt.plot(m_space,oc_T310_g_N,'g')
oc_T310_g, = plt.plot(m_space,oc_T310_g,'--g')

plt.plot(m_4,p_4,'s',color='none',markeredgecolor='black')
# oc_T323_g_S, = plt.plot(m_4_E,p_4_E,'--k')
oc_T323_g_N, = plt.plot(m_space,oc_T323_g_N,'k')
oc_T323_g, = plt.plot(m_space,oc_T323_g,'--k')

plt.plot(m_5,p_5,'s',color='none',markeredgecolor='magenta')
# oc_T348_g_S, = plt.plot(m_5_E,p_5_E,'--m')
oc_T348_g_N, = plt.plot(m_space,oc_T348_g_N,'m')
oc_T348_g, = plt.plot(m_space,oc_T348_g,'--m')

plt.plot(m_5,p_5,'s',color='none',markeredgecolor='yellow')
# oc_T373_g_S, = plt.plot(m_5_E,p_5_E,'--y')
oc_T373_g_N, = plt.plot(m_space,oc_T373_g_N,'y')
oc_T373_g, = plt.plot(m_space,oc_T373_g,'--y')

plt.xlabel(r'$\it{m}/ mol kg^{-1}$',fontsize=20)
plt.ylabel(r'$\Phi$',fontsize=20)
plt.legend([oc_T288_g_N,oc_T298_g_N,oc_T310_g_N,oc_T323_g_N,oc_T348_g_N,oc_T373_g_N],['288K','298K','310K','323K','348K','373K'],loc="upper left",fontsize=20)
axes = plt.gca()
axes.set_xlim([0,6])
axes.set_ylim([0.85,1.3])
size = 15
plt.xticks(fontsize=size)
plt.yticks(fontsize=size)
plt.tight_layout()
plt.show() 