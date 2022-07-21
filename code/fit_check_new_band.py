import numpy as np
import matplotlib.pyplot as plt
import pickle

import new_fit33 as new_fit
# import Inputs_ATFL_influx_slope_us_problem as problem
# import PD_Inputs_ATFL_influx_slope_us_problem as problem
# import PDL1_Integral_Inputs_ATFL_influx_slope_us_problem as problem
import PD_PDL1_Integral_Inputs_ATFL_influx_slope_us_problem as problem
# problem = new_fit.problem

import datetime
import shutil
import inspect
cur_timestamp = datetime.datetime.strftime(datetime.datetime.now(),'%Y%m%d_%H%M%S')
export_dir = 'fig_exp.d/fit_check/'+'exp_'+cur_timestamp+'_'
log_dir = 'log.d/fit_check/'+'log_'+cur_timestamp+'_'
def dubf():
    pass 
shutil.copyfile(inspect.getsourcefile(dubf),log_dir+'fit_check.py.txt')
shutil.copyfile(inspect.getsourcefile(new_fit.fit),log_dir+'fit_check_new_fit.py.txt')
shutil.copyfile(inspect.getsourcefile(problem.problem_input),log_dir+'problem.py.txt')


draw_track=True
draw_pdl1=False
draw_pd1 = False
draw_sens=False
draw_datamask=False

def get_params():
    # fname = 'fig_exp.d/new_fit/exp_20210123_202600_res.bin'
    # fname = 'fig_exp.d/new_fit/exp_20210126_202647_res.bin'
    fname = 'fig_exp.d/new_fit/exp_20210202_210953_res.bin'
    with open(fname,'rb') as f:
        params = pickle.load(f)[1]


    # PRE_Modification
    for p in params:
        # p['e']=p['e']*5
        # p['wTGFb']=1
        # p['wIL10']=1
        # p['wTNFa']=0
        # p['wIFNg']=0
        # p['wm1']=1
        # p['wm2']=1
        # p['r2']=1
        # p['DPDL_o_sPDL']=1
        # p['k3']=1250
        p['ATFLin']=None
        p['sPD_variation']=None
        p['DPDL1os_variation']=None
        pass
    ### Only works for certain subtype division
    for p in params[2:4]:
        p['DPDL_o_sPDL_factor'] = 0.408/p['DPDL_o_sPDL_base']

    return params

def get_cCSF_o_DCSF(CD163_0, p):
    k12 = p['k12']
    k21 = p['k21_2_k12'] * k12
    cCSF_o_ref = p['cCSF_o_ref']
    Tinh_CSF_0 = 1
    x = CD163_0/100
    UCSF_0 = 1-(1-x)/x * k12/k21
    cCSF_o_DCSF = UCSF_0/(Tinh_CSF_0*cCSF_o_ref-UCSF_0)
    return cCSF_o_DCSF

def draw_pats():

    pat2fit = [0,2,4,5,6,8]
    pat2check =[1,3,7]

    exp_data = new_fit.exp_data
    exp_errorbar = new_fit.original_exp_error
    params = get_params()
    # p0 = params[0]
    paramD = {}
    paramD_up = {}
    paramD_dn = {}

    for pat in pat2check:
        
        sbtp =  new_fit.subtype_of_total_pat[pat]
        model_pat = new_fit.subtype_fit_list[sbtp][0] 
        p0 = params[ model_pat ]
        print('Patient {} is initialized with parameter of fitting pat {}, in subtype {}'.format(str(pat),str(model_pat),str(sbtp)))
        data_CD163 = exp_data[pat,problem.FitTreatment.Control,problem.ResultField.CD163]
        errorbar_CD163 = exp_errorbar[pat,problem.FitTreatment.Control,problem.ResultField.CD163]
        paramD[pat] = p0.copy()
        paramD[pat]['cCSF_o_DCSF'] = get_cCSF_o_DCSF(data_CD163,p0)
        paramD_up[pat] = p0.copy()
        paramD_up[pat]['cCSF_o_DCSF'] = get_cCSF_o_DCSF(data_CD163+errorbar_CD163,p0)
        paramD_dn[pat] = p0.copy()
        paramD_dn[pat]['cCSF_o_DCSF'] = get_cCSF_o_DCSF(max(data_CD163-errorbar_CD163,0),p0)

    for k1,pat in enumerate(pat2fit):
        paramD[pat] = params[k1]
    
    for kp in [*pat2check,*pat2fit]:
        plt.close('all')
        param = paramD[kp]
        param['detailed']=True
        r = problem.problem_input(**param)
        with open(export_dir+'Parameters_pat{}.txt'.format(str(kp)),'w') as f:
            f.write('# fixed parameters\n')
            for key in list(new_fit.param_fixed):
                f.write(str(key)+'\t\t\t\t'+str(param[key])+'\n')
            f.write('# global variable parameters\n')
            for key in list(new_fit.param_shared):
                f.write(str(key)+'\t\t\t\t'+str(param[key])+'\n')
            f.write('# subtype variable parameters\n')
            for key in list(new_fit.param_subtype):
                f.write(str(key)+'\t\t\t\t'+str(param[key])+'\n')    
            f.write('# individual variable parameters\n')
            for key in list(new_fit.param_individual):
                f.write(str(key)+'\t\t\t\t'+str(param[key])+'\n')

        if kp in pat2check:
            param_up = paramD_up[kp]
            param_up['detailed']=True
            r_up = problem.problem_input(**param_up)
            param_dn = paramD_dn[kp]
            param_dn['detailed']=True
            r_dn = problem.problem_input(**param_dn)

        with open(export_dir+'TcellInflux.txt','a') as ftin:
            for k2,treat in enumerate(problem.FitTreatment):
                tin_str = 'Pat_{}\'s estimated T cell influx amount for {} is {}'.format(str(kp),treat.name,r.treatment_tracks[treat].na_n_nn_1)
                print(tin_str)
                ftin.write(tin_str+'\n')


        # titles = ['apop','CD154','CD163','TNF-alpha','IFN-gamma','IL-10','TGF-beta']
        if draw_track:
            legends = [ x.name for x in problem.FitTreatment ]
            f_datas = lambda sol: [sol.CD163, sol.CD154, sol.apop, sol.TGFb, sol.IL10, sol.TNFa, sol.IFNg]
            points = r.block_result
            tracks = r.treatment_tracks
            if kp in pat2check:
                tracks_up = r_up.treatment_tracks
                tracks_dn = r_dn.treatment_tracks
            for k1, fld in enumerate(problem.ResultField):
                fig = plt.figure(figsize=(5,5),dpi=72)
                sols = tracks
                if kp in pat2check:
                    sols_up = tracks_up
                    sols_dn = tracks_dn
                ys = [f_datas(sol)[k1] for sol in sols] 
                if kp in pat2check:
                    ys_up = [f_datas(sol)[k1] for sol in sols_up] 
                    ys_dn = [f_datas(sol)[k1] for sol in sols_dn] 
                ts = [sol.t for sol in sols]
                if kp in pat2check:
                    ts_up = [sol.t for sol in sols_up]
                    ts_dn = [sol.t for sol in sols_dn]
                # predict_point_treatments = predict_points[k1]
                for k2,treat in enumerate(problem.FitTreatment):
                    plt.plot(ts[treat],ys[treat])
                if kp in pat2check:
                    for k2,treat in enumerate(problem.FitTreatment):
                    ## Only works for homogeneous ts_up and ts_dn
                        plt.fill_between(ts_up[treat],ys_dn[treat],ys_up[treat],color='C'+str(k2),alpha=0.3)
                for k2,treat in enumerate(problem.FitTreatment):
                    plt.scatter(6,new_fit.exp_data[kp,treat,fld])
                    delta_errorbar_x = 0.05 # to avoid errorbar overlapping for a clearer view
                    plt.errorbar(6+k2*delta_errorbar_x,new_fit.exp_data[kp,treat,fld],yerr=new_fit.original_exp_error[kp,treat,fld],capsize=3,ecolor='C'+str(k2))
                    # print(predict_point_treatments)
                for k2,treat in enumerate(problem.FitTreatment):
                    plt.scatter(6,points[fld,treat],marker='*')
                plt.legend(legends)
                title = fld.name
                plt.title(title)
                plt.tight_layout()
                plt.savefig(export_dir+'pat_{}_{}.png'.format(str(kp),title))
                
                # headr = [x.name for x in problem.FitTreatment]
                # with open(export_dir+'pat_{}_{}.png'.format(str(kp),title), 'w') as f:
                #     for k2,treat in enumerate(problem.FitTreatment):

                ## Only works for homogeneous output
                assert(all([all(ts[0]==ts[x]) for x in range(problem.N_FitTreatment)]))
                headr = ['time', *[x.name for x in problem.FitTreatment]]
                out_arr = np.full([len(ts[0]),problem.N_FitTreatment+1],np.nan)
                out_arr[:,0] = ts[0]
                for k2,treat in enumerate(problem.FitTreatment):
                    out_arr[:,k2+1] = ys[treat][:]
                np.savetxt(export_dir+'pat_{}_{}.csv'.format(str(kp),title), out_arr,delimiter=', ',header=str(headr))
                        
                
            
            # Tcell Act plot
            fig = plt.figure(figsize=(5,5),dpi=72)
            title='na'
            sols = tracks
            ys = [sol.y[1] for sol in sols] 
            ts = [sol.t for sol in sols]
            for k2,treat in enumerate(problem.FitTreatment):
                plt.plot(ts[treat],ys[treat])
                # plt.plot(ts[treat][1:],(ys[treat][1:]-ys[treat][:-1])/(ts[treat][1:]-ts[treat][:-1]))
            plt.legend(legends)
            plt.title(title)
            plt.tight_layout()
            plt.savefig(export_dir+'pat_{}_{}.png'.format(str(kp),title))
                    ## Only works for homogeneous output
            assert(all([all(ts[0]==ts[x]) for x in range(problem.N_FitTreatment)]))
            headr = ['time', *[x.name for x in problem.FitTreatment]]
            out_arr = np.full([len(ts[0]),problem.N_FitTreatment+1],np.nan)
            out_arr[:,0] = ts[0]
            for k2,treat in enumerate(problem.FitTreatment):
                out_arr[:,k2+1] = ys[treat][:]
            np.savetxt(export_dir+'pat_{}_{}.csv'.format(str(kp),title), out_arr,delimiter=', ',header=str(headr))

            # Tcell Nonact plot
            fig = plt.figure(figsize=(5,5),dpi=72)
            title='nn'
            sols = tracks
            ys = [sol.y[2] for sol in sols] 
            ts = [sol.t for sol in sols]
            for k2,treat in enumerate(problem.FitTreatment):
                plt.plot(ts[treat],ys[treat])
            plt.legend(legends)
            plt.title(title)
            plt.tight_layout()
            plt.savefig(export_dir+'pat_{}_{}.png'.format(str(kp),title))
                    ## Only works for homogeneous output
            assert(all([all(ts[0]==ts[x]) for x in range(problem.N_FitTreatment)]))
            headr = ['time', *[x.name for x in problem.FitTreatment]]
            out_arr = np.full([len(ts[0]),problem.N_FitTreatment+1],np.nan)
            out_arr[:,0] = ts[0]
            for k2,treat in enumerate(problem.FitTreatment):
                out_arr[:,k2+1] = ys[treat][:]
            np.savetxt(export_dir+'pat_{}_{}.csv'.format(str(kp),title), out_arr,delimiter=', ',header=str(headr))



            # Tcell Nonact plot
            fig = plt.figure(figsize=(5,5),dpi=72)
            title='ATFL'
            sols = tracks
            ys = [sol.ATFL for sol in sols] 
            ts = [sol.t for sol in sols]
            for k2,treat in enumerate(problem.FitTreatment):
                plt.plot(ts[treat],ys[treat])
            plt.legend(legends)
            plt.title(title)
            plt.tight_layout()
            # plt.show()
            plt.savefig(export_dir+'pat_{}_{}.png'.format(str(kp),title))
                    ## Only works for homogeneous output
            assert(all([all(ts[0]==ts[x]) for x in range(problem.N_FitTreatment)]))
            headr = ['time', *[x.name for x in problem.FitTreatment]]
            out_arr = np.full([len(ts[0]),problem.N_FitTreatment+1],np.nan)
            out_arr[:,0] = ts[0]
            for k2,treat in enumerate(problem.FitTreatment):
                out_arr[:,k2+1] = ys[treat][:]
            np.savetxt(export_dir+'pat_{}_{}.csv'.format(str(kp),title), out_arr,delimiter=', ',header=str(headr))

            # Tcell Nonact plot
            fig = plt.figure(figsize=(5,5),dpi=72)
            title='ImSup_excPD'
            sols = tracks
            ys = [1-sol.ImSup_excPD for sol in sols] 
            ts = [sol.t for sol in sols]
            for k2,treat in enumerate(problem.FitTreatment):
                plt.plot(ts[treat],ys[treat])
            plt.legend(legends)
            plt.title(title)
            plt.tight_layout()
            # plt.show()
            plt.savefig(export_dir+'pat_{}_{}.png'.format(str(kp),title))
                    ## Only works for homogeneous output
            assert(all([all(ts[0]==ts[x]) for x in range(problem.N_FitTreatment)]))
            headr = ['time', *[x.name for x in problem.FitTreatment]]
            out_arr = np.full([len(ts[0]),problem.N_FitTreatment+1],np.nan)
            out_arr[:,0] = ts[0]
            for k2,treat in enumerate(problem.FitTreatment):
                out_arr[:,k2+1] = ys[treat][:]
            np.savetxt(export_dir+'pat_{}_{}.csv'.format(str(kp),title), out_arr,delimiter=', ',header=str(headr))

            # Tcell Nonact plot
            fig = plt.figure(figsize=(5,5),dpi=72)
            title='ImSup_incPD'
            sols = tracks
            ys = [1-sol.ImSup_incPD for sol in sols] 
            ts = [sol.t for sol in sols]
            for k2,treat in enumerate(problem.FitTreatment):
                plt.plot(ts[treat],ys[treat])
            plt.legend(legends)
            plt.title(title)
            plt.tight_layout()
            # plt.show()
            plt.savefig(export_dir+'pat_{}_{}.png'.format(str(kp),title))
                    ## Only works for homogeneous output
            assert(all([all(ts[0]==ts[x]) for x in range(problem.N_FitTreatment)]))
            headr = ['time', *[x.name for x in problem.FitTreatment]]
            out_arr = np.full([len(ts[0]),problem.N_FitTreatment+1],np.nan)
            out_arr[:,0] = ts[0]
            for k2,treat in enumerate(problem.FitTreatment):
                out_arr[:,k2+1] = ys[treat][:]
            np.savetxt(export_dir+'pat_{}_{}.csv'.format(str(kp),title), out_arr,delimiter=', ',header=str(headr))

            # Tcell Nonact plot
            fig = plt.figure(figsize=(5,5),dpi=72)
            title='ImAct'
            sols = tracks
            ys = [sol.ImAct for sol in sols] 
            ts = [sol.t for sol in sols]
            for k2,treat in enumerate(problem.FitTreatment):
                plt.plot(ts[treat],ys[treat])
            plt.legend(legends)
            plt.title(title)
            plt.tight_layout()
            # plt.show()
            plt.savefig(export_dir+'pat_{}_{}.png'.format(str(kp),title))
                    ## Only works for homogeneous output
            assert(all([all(ts[0]==ts[x]) for x in range(problem.N_FitTreatment)]))
            headr = ['time', *[x.name for x in problem.FitTreatment]]
            out_arr = np.full([len(ts[0]),problem.N_FitTreatment+1],np.nan)
            out_arr[:,0] = ts[0]
            for k2,treat in enumerate(problem.FitTreatment):
                out_arr[:,k2+1] = ys[treat][:]
            np.savetxt(export_dir+'pat_{}_{}.csv'.format(str(kp),title), out_arr,delimiter=', ',header=str(headr))

            # Tcell Nonact plot
            fig = plt.figure(figsize=(5,5),dpi=72)
            title='VTAM'
            sols = tracks
            ys = [1-sol.Us for sol in sols] 
            ts = [sol.t for sol in sols]
            for k2,treat in enumerate(problem.FitTreatment):
                plt.plot(ts[treat],ys[treat])
            plt.legend(legends)
            plt.title(title)
            plt.tight_layout()
            # plt.show()
            plt.savefig(export_dir+'pat_{}_{}.png'.format(str(kp),title))
                    ## Only works for homogeneous output
            assert(all([all(ts[0]==ts[x]) for x in range(problem.N_FitTreatment)]))
            headr = ['time', *[x.name for x in problem.FitTreatment]]
            out_arr = np.full([len(ts[0]),problem.N_FitTreatment+1],np.nan)
            out_arr[:,0] = ts[0]
            for k2,treat in enumerate(problem.FitTreatment):
                out_arr[:,k2+1] = ys[treat][:]
            np.savetxt(export_dir+'pat_{}_{}.csv'.format(str(kp),title), out_arr,delimiter=', ',header=str(headr))

            # Tcell Nonact plot
            fig = plt.figure(figsize=(5,5),dpi=72)
            title='CD8 (in TME)'
            sols = tracks
            ys = [sol.CD8 for sol in sols] 
            ts = [sol.t for sol in sols]
            for k2,treat in enumerate(problem.FitTreatment):
                plt.plot(ts[treat],ys[treat])
            plt.legend(legends)
            plt.title(title)
            plt.tight_layout()
            # plt.show()
            plt.savefig(export_dir+'pat_{}_{}.png'.format(str(kp),title))
                    ## Only works for homogeneous output
            assert(all([all(ts[0]==ts[x]) for x in range(problem.N_FitTreatment)]))
            headr = ['time', *[x.name for x in problem.FitTreatment]]
            out_arr = np.full([len(ts[0]),problem.N_FitTreatment+1],np.nan)
            out_arr[:,0] = ts[0]
            for k2,treat in enumerate(problem.FitTreatment):
                out_arr[:,k2+1] = ys[treat][:]
            np.savetxt(export_dir+'pat_{}_{}.csv'.format(str(kp),title), out_arr,delimiter=', ',header=str(headr))


            # Tcell Nonact plot
            fig = plt.figure(figsize=(5,5),dpi=72)
            title='CD154 (in TME)'
            sols = tracks
            ys = [sol.TME_CD154 for sol in sols] 
            ts = [sol.t for sol in sols]
            for k2,treat in enumerate(problem.FitTreatment):
                plt.plot(ts[treat],ys[treat])
            plt.legend(legends)
            plt.title(title)
            plt.tight_layout()
            # plt.show()
            plt.savefig(export_dir+'pat_{}_{}.png'.format(str(kp),title))
                    ## Only works for homogeneous output
            assert(all([all(ts[0]==ts[x]) for x in range(problem.N_FitTreatment)]))
            headr = ['time', *[x.name for x in problem.FitTreatment]]
            out_arr = np.full([len(ts[0]),problem.N_FitTreatment+1],np.nan)
            out_arr[:,0] = ts[0]
            for k2,treat in enumerate(problem.FitTreatment):
                out_arr[:,k2+1] = ys[treat][:]
            np.savetxt(export_dir+'pat_{}_{}.csv'.format(str(kp),title), out_arr,delimiter=', ',header=str(headr))


            # Tcell Nonact plot
            fig = plt.figure(figsize=(5,5),dpi=72)
            title='T cell influx'
            sols = tracks
            ys = [sol.influx for sol in sols] 
            ts = [sol.t for sol in sols]
            for k2,treat in enumerate(problem.FitTreatment):
                plt.plot(ts[treat],ys[treat])
            plt.legend(legends)
            plt.title(title)
            plt.tight_layout()
            # plt.show()
            plt.savefig(export_dir+'pat_{}_{}.png'.format(str(kp),title))
                    ## Only works for homogeneous output
            assert(all([all(ts[0]==ts[x]) for x in range(problem.N_FitTreatment)]))
            headr = ['time', *[x.name for x in problem.FitTreatment]]
            out_arr = np.full([len(ts[0]),problem.N_FitTreatment+1],np.nan)
            out_arr[:,0] = ts[0]
            for k2,treat in enumerate(problem.FitTreatment):
                out_arr[:,k2+1] = ys[treat][:]
            np.savetxt(export_dir+'pat_{}_{}.csv'.format(str(kp),title), out_arr,delimiter=', ',header=str(headr))


            # Tcell Nonact plot
            fig = plt.figure(figsize=(5,5),dpi=72)
            title='Tumor'
            sols = tracks
            ys = [sol.y[0] for sol in sols] 
            ts = [sol.t for sol in sols]
            for k2,treat in enumerate(problem.FitTreatment):
                plt.plot(ts[treat],ys[treat])
            plt.legend(legends)
            plt.title(title)
            plt.tight_layout()
            # plt.show()
            plt.savefig(export_dir+'pat_{}_{}.png'.format(str(kp),title))
                    ## Only works for homogeneous output
            assert(all([all(ts[0]==ts[x]) for x in range(problem.N_FitTreatment)]))
            headr = ['time', *[x.name for x in problem.FitTreatment]]
            out_arr = np.full([len(ts[0]),problem.N_FitTreatment+1],np.nan)
            out_arr[:,0] = ts[0]
            for k2,treat in enumerate(problem.FitTreatment):
                out_arr[:,k2+1] = ys[treat][:]
            np.savetxt(export_dir+'pat_{}_{}.csv'.format(str(kp),title), out_arr,delimiter=', ',header=str(headr))
    
# def draw_var_PD():
    # Variation of sPD-1 over sPD-1-ref
    if draw_pd1:
        pat2draw = 2
        p = paramD[pat2draw]
        p['detailed']=True
        timepoint = 7
        sPD_arr = np.linspace(0.5,1,21)
        sol_list =  []
        for sPD_var in sPD_arr:
            p['sPD_variation']=sPD_var
            s = problem.problem_input(**p)
            sol_list.append(s)
        ### only works for homogeneous time steps  
        kk = np.searchsorted(sol_list[0].treatment_tracks[problem.FitTreatment.Control].t,timepoint)

        fig = plt.figure(figsize=(5,5),dpi=72)
        title = 'Apoptosis ratio varing sPD1/sPD1ref,\n and fixing VTME for control'
        plt.plot(sPD_arr,[sol_list[kk1].treatment_tracks[problem.FitTreatment.Control].apop[kk] for kk1 in range(len(sPD_arr))],'o-')
        plt.plot(sPD_arr,[sol_list[kk1].treatment_tracks[problem.FitTreatment.aPD1].apop[kk] for kk1 in range(len(sPD_arr))],'o-')
        legends=['Control','aPD']
        plt.legend(legends)
        plt.xlabel('sPD1/sPD1ref')
        plt.ylabel('Apoptosis %')
        plt.title(title)
        plt.tight_layout()
        plt.savefig(export_dir+'var_sPD_{}.png'.format(str(pat2draw)))
        out_arr = np.full([3,len(sPD_arr)],np.nan)
        out_arr[0,:]=sPD_arr[:]
        out_arr[1,:]=np.array([sol_list[kk1].treatment_tracks[problem.FitTreatment.Control].apop[kk] for kk1 in range(len(sPD_arr))])
        out_arr[2,:]=np.array([sol_list[kk1].treatment_tracks[problem.FitTreatment.aPD1].apop[kk] for kk1 in range(len(sPD_arr))])
        hdr = ['sPD1/sPD1ref','control', 'aPD1']
        np.savetxt(export_dir+'var_sPD_{}.txt'.format(str(pat2draw)),out_arr.T, delimiter=', ',header=str(hdr))

        fig = plt.figure(figsize=(5,5),dpi=72)
        title = 'Apoptosis ratio (time series),\n varing sPD1/sPD1ref, and fixing VTME for control'
        legends=[]
        for kk2,sPD_var in enumerate(sPD_arr):
            sol = sol_list[kk2].treatment_tracks[problem.FitTreatment.aPD1]
            plt.plot(sol.t,sol.apop,color='b',alpha=(kk2+1)/len(sPD_arr))
            legends.append('{}%'.format(100*sPD_var))
        # legends=['Control','aPD']
        plt.legend(legends,title='sPD1/sPD1ref')
        plt.xlabel('time(day)')
        plt.ylabel('Apoptosis %')
        plt.title(title)
        plt.tight_layout()
        plt.savefig(export_dir+'var_sPD_{}_time.png'.format(str(pat2draw)))
        ### only works for homogeneous timestep
        t_arr = sol_list[0].treatment_tracks[problem.FitTreatment.aPD1].t
        out_arr = np.full([len(sPD_arr)+1, len(t_arr)],np.nan)
        out_arr[0,:]=t_arr
        for kk2,sPD_var in enumerate(sPD_arr):
            sol = sol_list[kk2].treatment_tracks[problem.FitTreatment.aPD1]
            out_arr[kk2+1,:]=sol.apop
        hdr = ['time',*legends]
        np.savetxt(export_dir+'var_sPD_{}_time.txt'.format(str(pat2draw)),out_arr.T, delimiter=', ',header=str(hdr))

# def draw_var_PDL1():
    if draw_pdl1:
    # Variation of sPDL1/DPDL1
        pat2draw = 2
        p = paramD[pat2draw]
        p['detailed']=True
        timepoint = 7
        sPDLoD_ratio_arr = np.linspace(1,5,21)
        DPDLos_ratio_arr = 1/sPDLoD_ratio_arr
        sol_list =  []
        old_DPDLos = p['DPDL_o_sPDL_base']*p['DPDL_o_sPDL_factor']
        old_sPDLoD = 1/old_DPDLos
        old_s = problem.problem_input(**p)
        for DPDLos_var in DPDLos_ratio_arr:
            # p['DPDL1os_variation']=DPDLos_var * old_DPDLos
            p['DPDL1os_variation']=DPDLos_var
            s = problem.problem_input(**p)
            sol_list.append(s)
        ### only works for homogeneous time steps  
        kk = np.searchsorted(sol_list[0].treatment_tracks[problem.FitTreatment.Control].t,timepoint)

        fig = plt.figure(figsize=(5,5),dpi=72)
        title = 'Apoptosis ratio varing sPDL1/KD-PD(L)1,\n and fixing VTME for control'
        plt.plot(sPDLoD_ratio_arr,[sol_list[kk1].treatment_tracks[problem.FitTreatment.Control].apop[kk] for kk1 in range(len(sPDLoD_ratio_arr))],'o-')
        plt.plot(sPDLoD_ratio_arr,[sol_list[kk1].treatment_tracks[problem.FitTreatment.aPD1].apop[kk] for kk1 in range(len(sPDLoD_ratio_arr))],'o-')
        plt.scatter(old_sPDLoD,old_s.treatment_tracks[problem.FitTreatment.aPD1].apop[kk],color='g')
        legends=['Control','aPD','This patient']
        plt.legend(legends)
        plt.xlabel('sPDL1/KD-PD(L)1')
        # plt.xlabel('sPDL1 fold change')
        plt.ylabel('Apoptosis %')
        plt.title(title)
        plt.tight_layout()
        plt.savefig(export_dir+'var_sPDL_{}.png'.format(str(pat2draw)))
        out_arr = np.full([3,len(sPDLoD_ratio_arr)],np.nan)
        out_arr[0,:]=sPDLoD_ratio_arr[:]
        out_arr[1,:]=np.array([sol_list[kk1].treatment_tracks[problem.FitTreatment.Control].apop[kk] for kk1 in range(len(sPDLoD_ratio_arr))])
        out_arr[2,:]=np.array([sol_list[kk1].treatment_tracks[problem.FitTreatment.aPD1].apop[kk] for kk1 in range(len(sPDLoD_ratio_arr))])
        hdr = ['sPDL1/KD-PD(L)1','control', 'aPD1']
        np.savetxt(export_dir+'var_sPDL_{}.txt'.format(str(pat2draw)),out_arr.T, delimiter=', ',header=str(hdr))

        fig = plt.figure(figsize=(5,5),dpi=72)
        title = 'Apoptosis ratio (time series),\n varing sPDL1/KD-PD(L)1, and fixing VTME for control'
        legends=[]
        for kk2,sPDLoD_var in enumerate(sPDLoD_ratio_arr):
            sol = sol_list[kk2].treatment_tracks[problem.FitTreatment.aPD1]
            plt.plot(sol.t,sol.apop,color='b',alpha=(kk2+1)/len(sPDLoD_ratio_arr))
            legends.append('{:.1f}'.format(sPDLoD_var))

        old_sol = old_s.treatment_tracks[problem.FitTreatment.aPD1]
        plt.plot(old_sol.t,old_sol.apop,color='r')
        legends.append('This patient')

        plt.legend(legends,title='sPDL1/KD-PD(L)1')
        # plt.legend(legends,title='sPDL1 fold change')
        plt.xlabel('time(day)')
        plt.ylabel('Apoptosis %')
        plt.title(title)
        plt.tight_layout()
        plt.savefig(export_dir+'var_sPDL_{}_time.png'.format(str(pat2draw)))
        ### only works for homogeneous timestep
        t_arr = sol_list[0].treatment_tracks[problem.FitTreatment.aPD1].t
        out_arr = np.full([len(sPDLoD_ratio_arr)+1, len(t_arr)],np.nan)
        out_arr[0,:]=t_arr
        for kk2,sPD_var in enumerate(sPDLoD_ratio_arr):
            sol = sol_list[kk2].treatment_tracks[problem.FitTreatment.aPD1]
            out_arr[kk2+1,:]=sol.apop
        hdr = ['time',*legends]
        np.savetxt(export_dir+'var_sPDL_{}_time.txt'.format(str(pat2draw)),out_arr.T, delimiter=', ',header=str(hdr))


# def sensitivity():
    if draw_sens:
        pat2draw=2
        pat2draw = 2
        p = paramD[pat2draw]
        p['detailed']=True
        timepoint = 7
        p_dev_factor = 0.01
        
        p_var_names = [*list(new_fit.param_shared),*list(new_fit.param_subtype),*list(new_fit.param_individual)]
        ### trim any parameter not needed for sensitivity analysis
        p_var_names.remove('DPDL_o_sPDL_base')
        p_var_names.remove('k4_base')
        p_var_names.remove('k7_base')

        # p_var_names_natural = get_natural_names(p_var_names)

        p_var_plus_list = np.full(len(p_var_names),np.nan)
        p_var_minus_list = np.full(len(p_var_names),np.nan)


        sol_pivot = problem.problem_input(**p).treatment_tracks[problem.FitTreatment.aPD1]
        kk = np.searchsorted(sol_pivot.t,timepoint)
        ### only works for homogeneous timestep
        apop_pivot = sol_pivot.apop[kk]

        for kk1,key in enumerate(p_var_names):
            p_plus = p.copy()
            p_minus = p.copy()

            p_value_pivot = p[key]
            p_value_plus = p_value_pivot*(1+p_dev_factor)
            p_value_minus = p_value_pivot*(1-p_dev_factor)

            p_plus[key] = p_value_plus
            p_minus[key] = p_value_minus

            apop_plus = problem.problem_input(**p_plus).treatment_tracks[problem.FitTreatment.aPD1].apop[kk]
            apop_minus = problem.problem_input(**p_minus).treatment_tracks[problem.FitTreatment.aPD1].apop[kk]

            rel_dev_apop_plus = (apop_plus-apop_pivot)/apop_pivot*100
            rel_dev_apop_minus = (apop_minus-apop_pivot)/apop_pivot*100
            p_var_plus_list[kk1]=rel_dev_apop_plus
            p_var_minus_list[kk1] = rel_dev_apop_minus
        
        ndxsort = np.argsort(p_var_plus_list)
        concise_ndxsort = np.concatenate((ndxsort[:5],ndxsort[-5:]))
        ndx_meaningful = []

        ndxs = [ndxsort,concise_ndxsort,ndx_meaningful]
        filename = ['sensitivity','sensitivity_concise','sensitivity_meaningful']
        
        for kk2,ndx in enumerate(ndxs):
            plt.figure(figsize=(7,7),dpi=72)
            y_pos = np.arange(len(ndx))
            plt.barh(y_pos,p_var_plus_list[ndx],color='r',alpha=0.5)
            plt.barh(y_pos,p_var_minus_list[ndx],color='b',alpha=0.5)
            ax = plt.gca()
            ax.set_yticks(y_pos)
            ax.set_yticklabels([p_var_names[x] for x in ndx])
            plt.legend(['+1%','-1%'])
            plt.xlabel('Percent variance of Apoptosis rate')
            plt.tight_layout()
            plt.savefig(export_dir+'{}_pat{}_apop.svg'.format(filename[kk2],pat2draw))
            plt.savefig(export_dir+'{}_pat{}_apop.png'.format(filename[kk2],pat2draw))
        
        def get_natural_names(p_var_names):
            ret = {}
            for name in p_var_names:
                ret[name]=name
            # ret['']
            return None



def data_scrutiny():
    pat2check =[1,3,7]
    exp_data_mask = new_fit.AllPSR_mask
    pat_check_mask = np.full_like(exp_data_mask,True)
    for pat in pat2check:
        pat_check_mask[pat,:,:]=False
    titles =  {'apop':'Apoptosis %', 'CD154':'CD8+CD154+ %', 'CD163':'CD163+ %'}
    ylabels = {'aPD1':'Nivolumab','aCSF1R':'BLZ945','Control':'Control','dual':'Nivolumab+BLZ945'}
    for fld in [problem.ResultField.CD154,problem.ResultField.CD163,problem.ResultField.apop]:
        plt.figure(figsize=(7,4),dpi=72)
        plt.imshow(exp_data_mask[:,:,fld].T==False,cmap='bwr', norm=plt.Normalize(vmin=-2,vmax=2))
        plt.imshow(pat_check_mask[:,:,fld].T==False,cmap='Oranges', norm=plt.Normalize(vmin=0,vmax=2),alpha=(pat_check_mask[:,:,fld].T==False).astype(float))
        ax = plt.gca()
        # plt.imshow(exp_data_mask[:,:,fld])
        # plt.show()
        ax.set_xticks(np.arange(new_fit.N_total_pat))
        ax.set_yticks(np.arange(len(problem.FitTreatment)))
        ax.set_xticklabels([str(x) for x in range(1,new_fit.N_total_pat+1)])
        ax.set_yticklabels([ylabels[x.name] for x in problem.FitTreatment])
        ax.xaxis.set_ticks_position('none')
        ax.yaxis.set_ticks_position('none')
        thetitle = titles[fld.name]
        plt.title(thetitle)
        plt.tight_layout()
        plt.savefig(export_dir+'data_scrutiny_{}.svg'.format(fld.name))
        plt.savefig(export_dir+'data_scrutiny_{}.png'.format(fld.name))

# def check_newpat():
#     patient2fit = [0,2]
#     patient2check = [1]
#     exp_data = new_fit.exp_data
#     exp_error = new_fit.exp_error
#     params = get_params()
#     p = params[0]
    


    





if __name__ == '__main__':
    draw_pats()
    if draw_datamask:
        data_scrutiny()
    # check_newpat()