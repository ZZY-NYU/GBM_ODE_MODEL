from scipy.integrate import solve_ivp,quad
import numpy as np
import matplotlib.pyplot as plt

from enum import IntEnum

class FitTreatment(IntEnum):
    Control = 0
    aPD1 = 1
    aCSF1R = 2
    dual = 3
N_FitTreatment = len(FitTreatment)

class ResultField(IntEnum):
    CD163 = 0
    CD154 = 1
    apop = 2
    TGFb = 3
    IL10 = 4
    TNFa = 5
    IFNg = 6
N_ResultField = len(ResultField)



class Cytokine_ProblemOutput:
    def __init__(self):
        self.block_result = np.full([N_ResultField, N_FitTreatment],np.nan)
        # self.monoclinic_CD154=True
        self.monoclinic_CD154=None
        # self.monoclinic_apop=True
        self.monoclinic_apop=None

class Detailed_Cytokine_ProblemOutput:
    def __init__(self):
        Cytokine_ProblemOutput.__init__(self)
        self.treatment_tracks = []
        for k1 in range(N_FitTreatment):
            self.treatment_tracks.append(Sol())
        self.fmu1t = None
    


# class ProblemOutput:
#     def __init__(self):
#         self.CD154_control = None
#         self.CD154_aPD1 = None
#         self.CD154_aCSF1 = None
#         self.CD154_dual = None

#         self.apop_control = None
#         self.apop_aPD1=  None
#         self.apop_aCSF1 = None
#         self.apop_dual = None

#         self.CD163_control = None
#         self.CD163_aCSF1 = None

# class ProblemDetailedOutput(ProblemOutput):
#     def __init__(self):
#         ProblemOutput.__init__(self)
#         self.sol_control = Sol()
#         self.sol_aPD1 = Sol()
#         self.sol_aCSF1 = Sol()
#         self.sol_dual = Sol()
#         self.get_csf_from_CD163=None
        



class Sol():
    def __init__(self):
        self.t = None
        self.y = None
        self.CD154 = None
        self.fCD163 = None
        self.CD163 = None
        self.apop = None
        self.TNFa = None
        self.TGFb = None
        self.IL10 = None
        self.IFNg = None
        self.ATFL = None
        self.influx = None
        self.CD8 = None # na+nn
        self.TME_CD154 = None # na/(na+nn)
        self.na_n_nn_1 = None # single number. PD1&control T cell influx amount (estimation)
        self.Us = None
        self.ImSup_excPD = None
        self.ImSup_incPD = None
        self.ImAct = None



def problem_input(
    r1,
    l1,
    # r2,
    r2_2_r1,
    l2,
    l3,
    e,
    kact_max,
    k12,
    k21_2_k12,
    k3,
    # k4,
    k4_base,
    k4_factor,
    k5,
    k6,
    # k7,
    k7_base,
    k7_factor,
    k8,
    k9,
    k10,
    k11,
    # DPDL_o_sPDL,
    DPDL_o_sPDL_base,
    DPDL_o_sPDL_factor,
    sPD_o_ref,
    cCSF_o_DCSF,
    cCSF_o_ref,
    # DCSF,
    Km2,
    Dag,
    KTGFb,
    KIL10,
    KTNFa,
    KIFNg,
    wTGFb,
    wIL10,
    wIFNg,
    wTNFa,
    wm2,
    rTGFb,
    rIL10,
    rIFNg,
    rTNFa,
    Tinh_CSF,
    Tinh_PD,
    nua_ratio,
    n0,
    na_n_nn_0,
    rCD154,
    tau,
    nM,
    wm1,
    Km1,
    KA,
    KU,
    Tt,
    ATFLin,
    DPDL1os_variation,
    sPD_variation,
    xinit,
    detailed: bool
):

    # FIXME Growth should be updated: not logistic
    # FIXME nuG/(nua+nuG) should be fixed

    DPDL_o_sPDL = DPDL_o_sPDL_base*DPDL_o_sPDL_factor
    k4 = k4_base * k4_factor
    k7=k7_base*k7_factor
    tend_influx = 3
    xUPD_null =[sPD_o_ref, DPDL_o_sPDL, lambda t: 1]

    # monoclinic_CD154=True
    # monoclinic_apop=True
    monoclinic_CD154=0
    monoclinic_apop=0
    monoclinic = [monoclinic_CD154,monoclinic_apop]

    r2 = r1 * r2_2_r1
    k21 = k21_2_k12*k12
    nG,m1_n_m2=[ n*n0 for n in xinit]
    na_n_nn = na_n_nn_0 * n0
    # nuM = n0*nM
    # m1 = m1_n_m2*km12*CD163_0


    # UCSF_0 = CSF_0/(CSF_0+KCSF)
    # UCSF_0 = cCSF_o_DCSF/(cCSF_o_DCSF+1)


    def get_inner_cytokines(nua,nun,mu1,mu2,nuG):
        # return
        cTGFb = k3*mu2 + k4*nuG + k5*mu2*nuG
        # cTGFb = k3*mu2 + k4*nuG
        cIL10 = k6*mu2 + k7*nuG + k8*mu2*nuG
        # cIL10 = k6*mu2 + k7*nuG

        cTNFa = k9*mu1 + k10* nua
        cIFNg = k11*nua

        return  cTGFb,cIL10,cTNFa,cIFNg
    
    def get_cytokines(nua,nun,mu1,mu2,nuG):
        # return 
        cTGFb,cIL10,cTNFa,cIFNg = get_inner_cytokines(nua,nun,mu1,mu2,nuG)
        return rTGFb*cTGFb, rIL10*cIL10, rTNFa*cTNFa, rIFNg*cIFNg


    def get_initial_Tcell_normalized(fm1,fm2,ng,na,T):
        def f_us(t):
            mu1=fm1(t)
            mu2=fm2(t)
            cTGFb,cIL10,cTNFa,cIFNg = get_inner_cytokines(na,None,mu1,mu2,ng)
            ## Production/ Bx
            # us_t = (1-wm2*mu2/(mu2+Km2))*(1+wm1*mu1/(mu1+Km1))* (1-wTGFb*cTGFb/(cTGFb+KTGFb)) * (1-wIL10* cIL10/(cIL10+KIL10)) * (1+ wTNFa* cTNFa/(cTNFa+KTNFa)) * (1+wIFNg* cIFNg/(cIFNg+KIFNg))
            ## Product/Bx, normalized
            # us_t = 1/8*(1-wm2*mu2/(mu2+Km2))*(1+wm1*mu1/(mu1+Km1))* (1-wTGFb*cTGFb/(cTGFb+KTGFb)) * (1-wIL10* cIL10/(cIL10+KIL10)) * (1+ wTNFa* cTNFa/(cTNFa+KTNFa)) * (1+wIFNg* cIFNg/(cIFNg+KIFNg))
            ## Product/Bx, truly normalized
            us_t = 1/((1+wm1)*(1+wTNFa)*(1+wIFNg))*(1-wm2*mu2/(mu2+Km2))*(1+wm1*mu1/(mu1+Km1))* (1-wTGFb*cTGFb/(cTGFb+KTGFb)) * (1-wIL10* cIL10/(cIL10+KIL10)) * (1+ wTNFa* cTNFa/(cTNFa+KTNFa)) * (1+wIFNg* cIFNg/(cIFNg+KIFNg))
            ## Production without baseline
            # us_t = 2*(1-wm2*mu2/(mu2+Km2))*(wm1*mu1/(mu1+Km1))* (1-wTGFb*cTGFb/(cTGFb+KTGFb)) * (1-wIL10* cIL10/(cIL10+KIL10)) * (wTNFa* cTNFa/(cTNFa+KTNFa)) * (wIFNg* cIFNg/(cIFNg+KIFNg))
            ## Summation
            # us_t = ( wm2*(1-mu2/(mu2+Km2))+wTGFb*(1-cTGFb/(cTGFb+KTGFb))+wIL10*(1-cIL10/(cIL10+KIL10)) )/(wm2+wTGFb+wIL10) * ( wm1*(1+mu1/(mu1+Km1))+wTNFa*(1+cTNFa/(cTNFa+KTNFa))+wIFNg*(1+cIFNg/(cIFNg+KIFNg)) )/(wm1+wTNFa+wIFNg)
            ## Complete Summation / summation exc. PD
            # us_t = (wm2*(1-mu2/(mu2+Km2))+wTGFb*(1-cTGFb/(cTGFb+KTGFb))+wIL10*(1-cIL10/(cIL10+KIL10)) + wm1*(1+mu1/(mu1+Km1))+wTNFa*(1+cTNFa/(cTNFa+KTNFa))+wIFNg*(1+cIFNg/(cIFNg+KIFNg)) )/(wm1+wTNFa+wIFNg+wm2+wTGFb+wIL10)
            ## Complete Summation w/o baseline
            # us_t = (wm2*(1-mu2/(mu2+Km2))+wTGFb*(1-cTGFb/(cTGFb+KTGFb))+wIL10*(1-cIL10/(cIL10+KIL10)) + wm1*(mu1/(mu1+Km1))+wTNFa*(cTNFa/(cTNFa+KTNFa))+wIFNg*(cIFNg/(cIFNg+KIFNg)) )/(wm1+wTNFa+wIFNg+wm2+wTGFb+wIL10)
            ## Cell only
            # us_t= (1-wm2*mu2/(mu2+Km2))*(1+wm1*mu1/(mu1+Km1))
            ## Cytokine only
            # us_t = (1-wTGFb*cTGFb/(cTGFb+KTGFb)) * (1-wIL10* cIL10/(cIL10+KIL10)) * (1+ wTNFa* cTNFa/(cTNFa+KTNFa)) * (1+wIFNg* cIFNg/(cIFNg+KIFNg))
            ## Root
            # us_t = us_t**(1/6)
            ## Michaelis
            # us_t = 2*us_t/(us_t+KA)
            # us_t = us_t/(us_t+KU)
            ## Tweaking
            # us_t = 1
            return us_t 
        kin_normalized = lambda t: f_us(t)*np.exp(-t/tau)
        nin_normalized = quad(kin_normalized,0,tend_influx)[0]
        return nin_normalized, kin_normalized

    

    def sysEq(t,x,fmu1,fmu2,Tinh_PD,monoclinic,fin,xUPD):
        nuG,nua,nun,nuGD=x
        inner_sPD_o_ref, inner_DPDL_o_sPDL, fUPD = xUPD
        mu1 = fmu1(t)
        mu2 = fmu2(t)
        dnin = fin(t)

        cTGFb,cIL10,cTNFa,cIFNg = get_cytokines(nua,nun,mu1,mu2,nuG)
        UPD = Tinh_PD*inner_sPD_o_ref* nuG / (nuG + inner_DPDL_o_sPDL)

        #FIXME This was the problem.
        ## Half
        # ATFL = (1-UPD)*(1-mu2/(m2+Km2))* (1-cTGFb/(cTGFb+KTGFb)) * (1-cIL10/(cIL10+KIL10))# * (2*cTNFa/(cTNFa+KTNFa)) * (2*cIFNg/(cIFNg+KIFNg))
        ## Product
        # ATFL = (1-UPD)*(1-wm2*mu2/(mu2+Km2))*(1+wm1*mu1/(mu1+Km1))* (1-wTGFb*cTGFb/(cTGFb+KTGFb)) * (1-wIL10* cIL10/(cIL10+KIL10)) * (1+ wTNFa* cTNFa/(cTNFa+KTNFa)) * (1+wIFNg* cIFNg/(cIFNg+KIFNg))
        ## Product, normalized
        # ATFL = 1/8*(1-UPD)*(1-wm2*mu2/(mu2+Km2))*(1+wm1*mu1/(mu1+Km1))* (1-wTGFb*cTGFb/(cTGFb+KTGFb)) * (1-wIL10* cIL10/(cIL10+KIL10)) * (1+ wTNFa* cTNFa/(cTNFa+KTNFa)) * (1+wIFNg* cIFNg/(cIFNg+KIFNg))
        ## Product, truly normalized
        ATFL = 1/((1+wm1)*(1+wTNFa)*(1+wIFNg))*(1-UPD)*(1-wm2*mu2/(mu2+Km2))*(1+wm1*mu1/(mu1+Km1))* (1-wTGFb*cTGFb/(cTGFb+KTGFb)) * (1-wIL10* cIL10/(cIL10+KIL10)) * (1+ wTNFa* cTNFa/(cTNFa+KTNFa)) * (1+wIFNg* cIFNg/(cIFNg+KIFNg))
        ## Product without baseline
        # ATFL = 2*(1-UPD)*(1-wm2*mu2/(mu2+Km2))*(wm1*mu1/(mu1+Km1))* (1-wTGFb*cTGFb/(cTGFb+KTGFb)) * (1-wIL10* cIL10/(cIL10+KIL10)) * ( wTNFa* cTNFa/(cTNFa+KTNFa)) * (wIFNg* cIFNg/(cIFNg+KIFNg))
        ## Summation
        # ATFL =( (1-UPD)+wm2*(1-mu2/(mu2+Km2))+wTGFb*(1-cTGFb/(cTGFb+KTGFb))+wIL10*(1-cIL10/(cIL10+KIL10)) )/(1+wm2+wTGFb+wIL10) * ( wm1*(1+mu1/(mu1+Km1))+wTNFa*(1+cTNFa/(cTNFa+KTNFa))+wIFNg*(1+cIFNg/(cIFNg+KIFNg)) )/(wm1+wTNFa+wIFNg)
        ## Complete Summation
        # ATFL =( (1-UPD)+wm2*(1-mu2/(mu2+Km2))+wTGFb*(1-cTGFb/(cTGFb+KTGFb))+wIL10*(1-cIL10/(cIL10+KIL10)) + wm1*(1+mu1/(mu1+Km1))+wTNFa*(1+cTNFa/(cTNFa+KTNFa))+wIFNg*(1+cIFNg/(cIFNg+KIFNg)) )/(wm1+wTNFa+wIFNg+1+wm2+wTGFb+wIL10)
        ## Summation exc. PD
        # ATFL =(1-UPD)*( wm2*(1-mu2/(mu2+Km2))+wTGFb*(1-cTGFb/(cTGFb+KTGFb))+wIL10*(1-cIL10/(cIL10+KIL10)) + wm1*(1+mu1/(mu1+Km1))+wTNFa*(1+cTNFa/(cTNFa+KTNFa))+wIFNg*(1+cIFNg/(cIFNg+KIFNg)) )/(wm1+wTNFa+wIFNg+wm2+wTGFb+wIL10)
        ## A + KA Bx Form
        # ATFL = ( (1-UPD) + KA*(1-wm2*mu2/(mu2+Km2))*(1+wm1*mu1/(mu1+Km1))* (1-wTGFb*cTGFb/(cTGFb+KTGFb)) * (1-wIL10* cIL10/(cIL10+KIL10)) * (1+ wTNFa* cTNFa/(cTNFa+KTNFa)) * (1+wIFNg* cIFNg/(cIFNg+KIFNg)) )/(1+KA)
        ## Complete Summation w/o baseline
        # ATFL =( (1-UPD)+wm2*(1-mu2/(mu2+Km2))+wTGFb*(1-cTGFb/(cTGFb+KTGFb))+wIL10*(1-cIL10/(cIL10+KIL10)) + wm1*(mu1/(mu1+Km1))+wTNFa*(cTNFa/(cTNFa+KTNFa))+wIFNg*(cIFNg/(cIFNg+KIFNg)) )/(wm1+wTNFa+wIFNg+1+wm2+wTGFb+wIL10)
        ## Cell only
        # ATFL = (1-UPD)*(1-wm2*mu2/(mu2+Km2))*(1+wm1*mu1/(mu1+Km1))
        ## Cytokine only
        # ATFL = (1-UPD)*(1-wTGFb*cTGFb/(cTGFb+KTGFb)) * (1-wIL10* cIL10/(cIL10+KIL10)) * (1+ wTNFa* cTNFa/(cTNFa+KTNFa)) * (1+wIFNg* cIFNg/(cIFNg+KIFNg))
        ## Root
        # ATFL = ATFL**(1/7)
        ## M-M
        # ATFL = 2*ATFL/(ATFL+KA)
        # ATFL = ATFL/(ATFL+KA)
        # print(t,ATFL)
        if not ATFLin == None: ATFL = ATFLin
        ## Killing regulated
        # dnuG = r1*nuG*(1-nuG/nuM)-ATFL*e*nua*nuG/(nua+nuG)-l1*nuG
        ## Killing not regulated
        # dnuG = r1*nuG*(1-nuG/nuM)-e*nua*nuG/(nua+nuG)-l1*nuG
        # dnuG = r1*nuG-e*nua*nuG/(nua+nuG)-l1*nuG
        dnuG = r1*nuG-e*nua*nuG/(nua+nun+nuG+mu1+mu2)-l1*nuG
        ## Killing regulated by PD-1 only
        # dnuG = r1*nuG*(1-nuG/nuM)-(1-UPD)*e*nua*nuG/(nua+nuG)-l1*nuG
        ## Proliferation not even regulated
        # dnua = r2*nua+ ATFL* nuG/(nuG+Dag) * kact_max * nun - l2*nua
        ## Proliferation not influenced by antigen
        dnua = nua_ratio*dnin+r2*nua+ fUPD(t)*ATFL* nuG/(nuG+Dag) * kact_max * nun - l2*nua
        ## Proliferation influenced by antigen
        # dnua = ATFL*r2*nua*nuG/(nuG+Dag) + ATFL* nuG/(nuG+Dag) * kact_max * nun - l2*nua
        # dnua = ATFL
        dnun = (1-nua_ratio)*dnin - fUPD(t)*ATFL* nuG/(nuG+Dag) * kact_max * nun-l3*nun
        #dnn = AIAR*(kpro-kdif)*nn-l3*nn
        dnuGD = -dnuG + r1*nuG
        # dnuGD = -dnuG + r1*nuG*(1-nuG/nuM)

        # Judge if monoclinic
        # monoclinic[0]*= (dnua>=0)
        # monoclinic[1]*= (dnuGD*nuG - nuGD*dnuG>=0)
        monoclinic[0] = np.min([0,monoclinic[0],dnua])
        monoclinic[1] = np.min([0,monoclinic[0],(dnuGD*nuG-nuGD*dnuG)/(nuG+nuGD)**2])

        return [dnuG,dnua,dnun,dnuGD]


    rtol, atol = (1e-8,1e-8)
    ## Day 0--3

    Tinh_CSF_0 = 1
    Tinh_PD_0 = 1

    UCSF_0 = Tinh_CSF_0 * cCSF_o_ref* cCSF_o_DCSF / (cCSF_o_DCSF + 1)
    m1 = k21*(1-UCSF_0)/(k21*(1-UCSF_0)+k12)*m1_n_m2
    m2 = m1_n_m2-m1
    nuG0,nua0,nun0 = nG, 0, 0
    fmu10 = lambda t: m1
    fmu20 = lambda t: m2
    fin0 = lambda t: 0
    mu10 = m1
    mu20 = m2
    nuGD0 = 0
    x0 = [nuG0,nua0,nun0,nuGD0]
    duration = [0,3]
    sol1 = solve_ivp(sysEq,duration,x0,rtol=rtol,atol=atol,args=(fmu10,fmu20,Tinh_PD_0,monoclinic,fin0,xUPD_null))
    xt1 = [y[-1] for y in sol1.y]

    ## Day 3--6
    
    nuGt1,nuat1,nunt1,nuGDt1 = xt1
    # xt2 = [nuGt1,na,nn,0] 
    durationt2 = [0,3]
    Tinh_CSF_2 = Tinh_CSF
    Tinh_PD_2 = Tinh_PD
    UCSF_2 = Tinh_CSF_2 * cCSF_o_ref * cCSF_o_DCSF / (cCSF_o_DCSF + 1)
    k21new = (1-UCSF_2)*k21
    fmu1t = lambda t: (m1+m2)*k21new/(k21new+k12)-(m2*k21new-m1*k12)/(k21new+k12)*np.exp(-(k21new+k12)*t)
    fmu2t = lambda t: m1_n_m2-fmu1t(t)


    nin_control,kin_control = get_initial_Tcell_normalized(fmu10,fmu20,nuGt1,0,tau)
    nin_treat,kin_treat = get_initial_Tcell_normalized(fmu1t,fmu2t,nuGt1,0,tau)

    ## fixing control influx
    # nua_n_nun_control = na_n_nn
    # # nua_n_nun_treat = na_n_nn*nin_treat/nin_control 
    # ## To suppress error when nin_control=0
    # nua_n_nun_treat = na_n_nn*nin_treat/nin_control if nin_control and nin_treat else 1e-4

    ## fixing treat influx may be better
    # nua_n_nun_treat = na_n_nn
    # nua_n_nun_control = na_n_nn/nin_treat*nin_control if nin_control and nin_treat else 1e-4

    if nin_treat: 
        fin_treat = lambda t: kin_treat(t)/nin_treat*na_n_nn
        fin_control = lambda t: kin_control(t)/nin_treat*na_n_nn
    else:   
        fin_treat = lambda t: 0
        fin_control = lambda t: 0
    


    xt2_con = [nuGt1,0,0,0]
    xt2_tre = [nuGt1,0,0,0]


    treat_sol2 = [None,None,None,None]
    treat_sol2[FitTreatment.Control] = solve_ivp(sysEq,duration,xt2_con,rtol=rtol,atol=atol,args=(fmu10,fmu20,Tinh_PD_0,monoclinic,fin_control,xUPD_null))
    treat_sol2[FitTreatment.aPD1] = solve_ivp(sysEq,duration,xt2_con,rtol=rtol,atol=atol,args=(fmu10,fmu20,Tinh_PD_2,monoclinic,fin_control,xUPD_null))
    treat_sol2[FitTreatment.aCSF1R] = solve_ivp(sysEq,duration,xt2_tre,rtol=rtol,atol=atol,args=(fmu1t,fmu2t,Tinh_PD_0,monoclinic,fin_treat,xUPD_null))
    treat_sol2[FitTreatment.dual] = solve_ivp(sysEq,duration,xt2_tre,rtol=rtol,atol=atol,args=(fmu1t,fmu2t,Tinh_PD_2,monoclinic,fin_treat,xUPD_null))

    ## what is needed for output
    # CD154 level
    if not detailed:
        ret = Cytokine_ProblemOutput()
        ret.monoclinic_CD154, ret.monoclinic_apop= monoclinic
    else:
        ret = Detailed_Cytokine_ProblemOutput()
        # treat_sol2[FitTreatment.dual] = solve_ivp(sysEq,duration,xt2_tre,rtol=rtol,atol=atol,args=(fmu1t,fmu2t,Tinh_PD_2,monoclinic))

    blk = ret.block_result
    mu1fin = fmu1t(durationt2[-1])
    mu2fin = fmu2t(durationt2[-1])
    # for treat in FitTreatment:
    for treat in [FitTreatment.Control,FitTreatment.aPD1]:
        sol = treat_sol2[treat]
        if sol:
            # blk[ResultField.CD154,treat] = sol.y[1][-1]/(sol.y[1][-1]+sol.y[2][-1])*100
            blk[ResultField.CD154,treat] = sol.y[1][-1]/(sol.y[1][-1]+sol.y[2][-1]+Tt-nin_control)*100
            blk[ResultField.apop,treat]  = sol.y[3][-1]/(sol.y[3][-1]+sol.y[0][-1])*100
            sol_TGFb,sol_IL10,sol_TNFa,sol_IFNg  = get_cytokines(sol.y[1][-1],sol.y[2][-1],m1,m2,sol.y[0][-1])
            blk[[ResultField.TGFb, ResultField.IL10, ResultField.TNFa, ResultField.IFNg],treat]=(sol_TGFb,sol_IL10,sol_TNFa,sol_IFNg)

    for treat in [FitTreatment.aCSF1R,FitTreatment.dual]:
        sol = treat_sol2[treat]
        if sol:
            # blk[ResultField.CD154,treat] = sol.y[1][-1]/(sol.y[1][-1]+sol.y[2][-1])*100
            blk[ResultField.CD154,treat] = sol.y[1][-1]/(sol.y[1][-1]+sol.y[2][-1]+Tt-nin_treat)*100
            blk[ResultField.apop,treat]  = sol.y[3][-1]/(sol.y[3][-1]+sol.y[0][-1])*100
            sol_TGFb,sol_IL10,sol_TNFa,sol_IFNg  = get_cytokines(sol.y[1][-1],sol.y[2][-1],mu1fin,mu2fin,sol.y[0][-1])
            blk[[ResultField.TGFb, ResultField.IL10, ResultField.TNFa, ResultField.IFNg],treat]=(sol_TGFb,sol_IL10,sol_TNFa,sol_IFNg)
    blk[ResultField.CD163, FitTreatment.Control] = m2/m1_n_m2 * 100
    
    blk[ResultField.CD163, FitTreatment.aCSF1R] = mu2fin/m1_n_m2 * 100
    blk[ResultField.CD163, FitTreatment.dual] = mu2fin/m1_n_m2 * 100

    # blk[ResultField.CD154,:]*=rCD154

    # assert(np.all(np.isfinite(blk)))

    # ret.CD154_control = sol2_control.y[1][-1]/(sol2_control.y[1][-1]+sol2_control.y[2][-1])*100
    

    # ret.CD163_control = m2/m1_n_m2*100
    # ret.CD163_aCSF1 = fmu2t(durationt2[1])/m1_n_m2*100


    if not detailed:
        return ret
    else:
        import scipy.interpolate as inp
        time_step_exp = 0.2
        duration_exp2 = [0,4]
        time_eval2 = [ *np.arange(*duration_exp2,time_step_exp), duration_exp2[1] ]
        time_eval =  [ *np.arange(*duration,time_step_exp), duration[1] ]
        sol1 = solve_ivp(sysEq,duration,x0,t_eval=time_eval,rtol=rtol,atol=atol,args=(fmu10,fmu20,Tinh_PD_0,monoclinic,fin0,xUPD_null))
        treat_sol2[FitTreatment.Control] = solve_ivp(sysEq,duration_exp2,xt2_con,t_eval=time_eval2,rtol=rtol,atol=atol,args=(fmu10,fmu20,Tinh_PD_0,monoclinic,fin_control,xUPD_null))
        if DPDL1os_variation is not None or sPD_variation is not None:
            if DPDL1os_variation is not None:
                new_DPDL1os = DPDL1os_variation
            else:
                new_DPDL1os = DPDL_o_sPDL
            if sPD_variation is not None:
                new_sPD = sPD_variation
            else:
                new_sPD = sPD_o_ref

            t_2_arr = treat_sol2[FitTreatment.Control].t
            nuG_control_arr = treat_sol2[FitTreatment.Control].y[0]
            UPD_old_arr = Tinh_PD_0*sPD_o_ref* nuG_control_arr / (nuG_control_arr + DPDL_o_sPDL)
            UPD_variation_arr = Tinh_PD_0*new_sPD* nuG_control_arr / (nuG_control_arr + new_DPDL1os)
            xUPD_var = [new_sPD,new_DPDL1os, inp.interp1d(t_2_arr,(1-UPD_old_arr)/(1-UPD_variation_arr),kind='cubic')]
        else:
            xUPD_var=xUPD_null
        treat_sol2[FitTreatment.aPD1] = solve_ivp(sysEq,duration_exp2,xt2_con,t_eval=time_eval2,rtol=rtol,atol=atol,args=(fmu10,fmu20,Tinh_PD_2,monoclinic,fin_control,xUPD_var))
        treat_sol2[FitTreatment.aCSF1R] = solve_ivp(sysEq,duration_exp2,xt2_tre,t_eval=time_eval2,rtol=rtol,atol=atol,args=(fmu1t,fmu2t,Tinh_PD_0,monoclinic,fin_treat,xUPD_var))
        treat_sol2[FitTreatment.dual] = solve_ivp(sysEq,duration_exp2,xt2_tre,t_eval=time_eval2,rtol=rtol,atol=atol,args=(fmu1t,fmu2t,Tinh_PD_2,monoclinic,fin_treat,xUPD_var))

        ret.fmu1t = fmu1t

        for sol in (sol1,*treat_sol2):
            for ndx in [0,1,2,3]:
                sol.y[ndx] /= n0    
        # for sol in (sol1,sol2_aCSF1,sol2_aPD1,sol2_control,sol2_dual):
        #     for ndx in [0,1,2,6]: # index of number densities.
        #         sol.y[ndx] /= n0

        
        sol2s = treat_sol2
        tys = ret.treatment_tracks
        csf_ins = [None,None,None,None]
        csf_ins[FitTreatment.Control]=False; csf_ins[FitTreatment.aPD1]=False; csf_ins[FitTreatment.aCSF1R]=True; csf_ins[FitTreatment.dual]=True
        pd1_ins = [None,None,None,None]
        pd1_ins[FitTreatment.Control]=False; pd1_ins[FitTreatment.aPD1]=True; pd1_ins[FitTreatment.aCSF1R]=False; pd1_ins[FitTreatment.dual]=True
        for sol2,ty,csf_in,pd1_in in zip(sol2s,tys,csf_ins,pd1_ins):
            ty.t = np.concatenate([sol1.t,sol2.t+duration[-1]])
            ty.y = []
            for k_y in range(len(sol1.y)):
                ty.y.append(np.concatenate([sol1.y[k_y],sol2.y[k_y]]))
            # ty.CD154 =rCD154*np.concatenate([ np.full_like(sol1.y[1],np.nan), 100*sol2.y[1]/(sol2.y[1]+sol2.y[2]) ])
            
            if csf_in:
                ty.fCD163 = lambda t: 100*fmu2t(t-duration[1])/m1_n_m2 if t>duration[1] else 100*m2/m1_n_m2
                ty.CD154 =np.concatenate([ np.full_like(sol1.y[1],np.nan), 100*sol2.y[1]/(sol2.y[1]+sol2.y[2]+Tt-nin_treat) ])
                ty.na_n_nn_1 = na_n_nn
            else:
                ty.fCD163 = lambda t: 100*fmu20(t-duration[1])/m1_n_m2 if t>duration[1] else 100*m2/m1_n_m2  
                ty.CD154 =np.concatenate([ np.full_like(sol1.y[1],np.nan), 100*sol2.y[1]/(sol2.y[1]+sol2.y[2]+Tt-nin_control) ])
                ty.na_n_nn_1 = nin_control/nin_treat*na_n_nn
        
            complete_nuGD = np.concatenate([sol1.y[3],sol2.y[3][1:]+sol1.y[3][-1]])
            fnuGD = inp.interp1d(np.unique(ty.t),complete_nuGD,kind='cubic')
            nuGD_t = ty.t.copy()-3
            nuGD_t[nuGD_t<0]=0
            nuGD = np.array([fnuGD(tf)-fnuGD(ti) for tf,ti in zip(ty.t,nuGD_t)])
            ty.apop = nuGD/(nuGD+ty.y[0])*100

            ty.CD163 =np.full(len(ty.t),np.nan) 
            ty.IL10 = np.full(len(ty.t),np.nan)
            ty.TGFb = np.full(len(ty.t),np.nan)
            ty.TNFa = np.full(len(ty.t),np.nan)
            ty.IFNg = np.full(len(ty.t),np.nan)
            ty.ATFL = np.full(len(ty.t),np.nan)
            ty.influx = np.full(len(ty.t),np.nan)
            ty.CD8 = np.full(len(ty.t),np.nan)
            ty.TME_CD154 = np.full(len(ty.t),np.nan)
            ty.Us = np.full(len(ty.t),np.nan)
            ty.ImSup_excPD = np.full(len(ty.t),np.nan)
            ty.ImSup_incPD = np.full(len(ty.t),np.nan)
            ty.ImAct = np.full(len(ty.t),np.nan)

            boundary_k = np.nonzero(ty.t==duration[1])[0]
            for k1,t in enumerate(ty.t):
                is_in_duration2= (k1>=boundary_k[1])

                ty.CD163[k1] = ty.fCD163(t)
                nuG,nua,nun,nuGD = [ty.y[k2][k1] for k2 in range(4)]

                ty.CD8[k1] = nua+nun
                ty.TME_CD154[k1] = nua/(nua+nun) * 100

                if csf_in:
                    mu2 = fmu2t(t-duration[1]) if is_in_duration2 else m2
                    if is_in_duration2:ty.influx[k1] = kin_treat(t-duration[1]) 
                else:
                    mu2 = m2
                    if is_in_duration2:ty.influx[k1] = kin_control(t-duration[1]) 
                mu1 = m1_n_m2-mu2
                cTGFb,cIL10,cTNFa,cIFNg = get_cytokines(nua,nun,mu1,mu2,nuG)
                ty.IL10[k1]= cIL10; ty.TGFb[k1]=cTGFb; ty.TNFa[k1]=cTNFa; ty.IFNg[k1]=cIFNg
                Tinh_PD_tmp = None
                if is_in_duration2:
                    if pd1_in:
                        Tinh_PD_tmp = Tinh_PD_2
                    else:
                        Tinh_PD_tmp = Tinh_PD_0
                else:
                    Tinh_PD_tmp = Tinh_PD_0
                # print(Tinh_PD_tmp)
                UPD = Tinh_PD_tmp*sPD_o_ref* nuG / (nuG + DPDL_o_sPDL)
                # ATFL = (1-UPD)*(1-wm2*mu2/(mu2+Km2))*(1+wm1*mu1/(mu1+Km1))* (1-wTGFb*cTGFb/(cTGFb+KTGFb)) * (1-wIL10* cIL10/(cIL10+KIL10)) * (1+ wTNFa* cTNFa/(cTNFa+KTNFa)) * (1+wIFNg* cIFNg/(cIFNg+KIFNg))
                ## Product, normalized
                # ATFL = 1/8*(1-UPD)*(1-wm2*mu2/(mu2+Km2))*(1+wm1*mu1/(mu1+Km1))* (1-wTGFb*cTGFb/(cTGFb+KTGFb)) * (1-wIL10* cIL10/(cIL10+KIL10)) * (1+ wTNFa* cTNFa/(cTNFa+KTNFa)) * (1+wIFNg* cIFNg/(cIFNg+KIFNg))
                ## Product, truly normalized
                ATFL = 1/((1+wm1)*(1+wTNFa)*(1+wIFNg))*(1-UPD)*(1-wm2*mu2/(mu2+Km2))*(1+wm1*mu1/(mu1+Km1))* (1-wTGFb*cTGFb/(cTGFb+KTGFb)) * (1-wIL10* cIL10/(cIL10+KIL10)) * (1+ wTNFa* cTNFa/(cTNFa+KTNFa)) * (1+wIFNg* cIFNg/(cIFNg+KIFNg))
                # ATFL =( (1-UPD)+wm2*(1-mu2/(mu2+Km2))+wTGFb*(1-cTGFb/(cTGFb+KTGFb))+wIL10*(1-cIL10/(cIL10+KIL10)) )/(1+wm2+wTGFb+wIL10) * ( wm1*(1+mu1/(mu1+Km1))+wTNFa*(1+cTNFa/(cTNFa+KTNFa))+wIFNg*(1+cIFNg/(cIFNg+KIFNg)) )/(wm1+wTNFa+wIFNg)
                ## Complete Summation
                # ATFL =( (1-UPD)+wm2*(1-mu2/(mu2+Km2))+wTGFb*(1-cTGFb/(cTGFb+KTGFb))+wIL10*(1-cIL10/(cIL10+KIL10)) + wm1*(1+mu1/(mu1+Km1))+wTNFa*(1+cTNFa/(cTNFa+KTNFa))+wIFNg*(1+cIFNg/(cIFNg+KIFNg)) )/(wm1+wTNFa+wIFNg+1+wm2+wTGFb+wIL10)
                ## Summation exc. PD
                # ATFL =(1-UPD)*( wm2*(1-mu2/(mu2+Km2))+wTGFb*(1-cTGFb/(cTGFb+KTGFb))+wIL10*(1-cIL10/(cIL10+KIL10)) + wm1*(1+mu1/(mu1+Km1))+wTNFa*(1+cTNFa/(cTNFa+KTNFa))+wIFNg*(1+cIFNg/(cIFNg+KIFNg)) )/(wm1+wTNFa+wIFNg+wm2+wTGFb+wIL10)
                ## A + KA Bx Form
                # ATFL = ( (1-UPD) + KA*(1-wm2*mu2/(mu2+Km2))*(1+wm1*mu1/(mu1+Km1))* (1-wTGFb*cTGFb/(cTGFb+KTGFb)) * (1-wIL10* cIL10/(cIL10+KIL10)) * (1+ wTNFa* cTNFa/(cTNFa+KTNFa)) * (1+wIFNg* cIFNg/(cIFNg+KIFNg)) )/(1+KA)
                ## Cell only
                # ATFL = (1-UPD)*(1-wm2*mu2/(mu2+Km2))*(1+wm1*mu1/(mu1+Km1))
                ## Cytokine only
                # ATFL = (1-UPD)* (1-wTGFb*cTGFb/(cTGFb+KTGFb)) * (1-wIL10* cIL10/(cIL10+KIL10)) * (1+ wTNFa* cTNFa/(cTNFa+KTNFa)) * (1+wIFNg* cIFNg/(cIFNg+KIFNg))
                ## M-M
                # ATFL = 2*ATFL/(ATFL+KA)
                # ATFL = ATFL/(ATFL+KA)
                ## Root
                # ATFL = ATFL**(1/7)
                if not ATFLin == None: ATFL = ATFLin
                ty.ATFL[k1]=ATFL
                ### Us
                us_t = 1/((1+wm1)*(1+wTNFa)*(1+wIFNg))*(1-wm2*mu2/(mu2+Km2))*(1+wm1*mu1/(mu1+Km1))* (1-wTGFb*cTGFb/(cTGFb+KTGFb)) * (1-wIL10* cIL10/(cIL10+KIL10)) * (1+ wTNFa* cTNFa/(cTNFa+KTNFa)) * (1+wIFNg* cIFNg/(cIFNg+KIFNg))
                ty.Us[k1] = us_t
                ### Immnunosuppressive factor excluding PD
                ImSup = (1-wm2*mu2/(mu2+Km2))* (1-wTGFb*cTGFb/(cTGFb+KTGFb)) * (1-wIL10* cIL10/(cIL10+KIL10))
                ty.ImSup_excPD[k1] = ImSup

                ### Immnunosuppressive factor including PD
                ty.ImSup_incPD[k1] = ImSup * (1-UPD)

                ### Immunoactive factor
                ImAct = 1/((1+wm1)*(1+wTNFa)*(1+wIFNg))*(1+wm1*mu1/(mu1+Km1))* (1+ wTNFa* cTNFa/(cTNFa+KTNFa)) * (1+wIFNg* cIFNg/(cIFNg+KIFNg))
                ty.ImAct[k1] = ImAct

        return ret
        





    









