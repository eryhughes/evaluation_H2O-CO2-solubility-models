# import python packages
import numpy as np
import pandas as pd
import tqdm
import VESIcal as vc

########################################
### measured parameters within error ### 
########################################

def comp_within_error_one_sample(run,setup,iterations):

    for iter in range(iterations): 
        
        for x in ['T_C', 'P_bar', 'H2O_wt%', 'CO2_ppm', 'SiO2_wt%','TiO2_wt%','Al2O3_wt%','MnO_wt%','MgO_wt%','FeOT_wt%','CaO_wt%','Na2O_wt%','K2O_wt%','P2O5_wt%','Fe3_FeT']:
            
            if setup.loc[run,x+' sd type'] == 'NR':
                if x == 'T_C':
                    setup.loc[run,x+' sd'] = 10. # 'C
                    setup.loc[run,x+' sd type'] = 'A'
                elif x == 'P_bar':
                    if setup.loc[run,'Experimental apparatus'] == 'IHPV':
                        setup.loc[run,x+' sd'] = 50 # bar
                    elif setup.loc[run,'Experimental apparatus'] == 'PC':
                        setup.loc[run,x+' sd'] = 100 # bar
                    elif setup.loc[run,'Experimental apparatus'] == 'TZM':
                        setup.loc[run,x+' sd'] = 50 # bar
                    setup.loc[run,x+' sd type'] = 'A'
                elif x == 'H2O_wt%' or x == 'CO2_ppm' :
                    if setup.loc[run,x] == 0.:
                        setup.loc[run,x+' sd'] = 0.
                        setup.loc[run,x+' sd type'] = 'A'
                    else:
                        setup.loc[run,x+' sd'] = 0.1
                        setup.loc[run,x+' sd type'] = 'R'
                elif x == 'SiO2_wt%' or x == 'TiO2_wt%' or x == 'Al2O3_wt%' or x == 'MnO_wt%' or x == 'MgO_wt%' or x == 'FeOT_wt%' or x == 'CaO_wt%' or x == 'Na2O_wt%' or x == 'K2O_wt%' or x == 'P2O5_wt%':
                    setup.loc[run,x+' sd'] = 0.01
                    setup.loc[run,x+' sd type'] = 'R'
                elif x == 'Fe3_FeT':
                    setup.loc[run,x+' sd'] = 0.03
                    setup.loc[run,x+' sd type'] = 'A'                   
            if setup.loc[run,x+' sd type'] == "A": # absolute
                sd = setup.loc[run,x+" sd"]
            else: # relative
                sd = setup.loc[run,x+" sd"]*setup.loc[run,x]

            if setup.loc[run,x+' sd type'] == 'ND':
                x_value = 0.
            else:
                x_value = float(np.random.normal(setup.loc[run,x],sd,1))

            if x == 'T_C':
                if x_value < -273.15:
                    x_value = -273.15
            else:
                if x_value < 0:
                    x_value = 0
            if x == 'T_C':
                T_C = x_value
            elif x == 'P_bar':
                P_bar = x_value
            elif x == 'H2O_wt%':
                H2O = x_value
            elif x == 'CO2_ppm':
                CO2 = x_value
            elif x == 'SiO2_wt%':
                SiO2 = x_value
            elif x == 'TiO2_wt%':
                TiO2 = x_value
            elif x == 'Al2O3_wt%':
                Al2O3 = x_value
            elif x == 'MnO_wt%':
                MnO = x_value
            elif x == 'MgO_wt%':
                MgO = x_value
            elif x == 'FeOT_wt%':
                FeOT = x_value
            elif x == 'CaO_wt%':
                CaO = x_value
            elif x == 'Na2O_wt%':
                Na2O = x_value
            elif x == 'K2O_wt%':
                K2O = x_value
            elif x == 'P2O5_wt%':
                P2O5 = x_value
            elif x == 'Fe3_FeT':
                Fe3FeT = x_value

            
        results1 = pd.DataFrame([[P_bar, T_C, H2O, CO2, SiO2, TiO2, Al2O3, MnO, MgO, FeOT, CaO, Na2O, K2O, P2O5,Fe3FeT]])
            
        if iter == 0:
            results = results1
        else:
            results = pd.concat([results, results1])
    
    results.columns = ["P_bar","T_C","H2O_wt%","CO2_ppm","SiO2_wt%","TiO2_wt%","Al2O3_wt%","MnO_wt%","MgO_wt%","FeOT_wt%","CaO_wt%","Na2O_wt%","K2O_wt%","P2O5_wt%",'Fe3_FeT']
    results = results.reset_index(drop=True)
    return results

#######################
### CALCULATE PVSAT ###
#######################

def VESIcal_MC(setup,iterations,model,normalisation='fixedvolatiles',redox='N'):
    
    with tqdm.tqdm(total=len(setup)) as tqdmsteps:
        for run in range(0,len(setup),1):
            random_comps = comp_within_error_one_sample(run,setup,iterations)
            random_comps_vc = random_comps.rename(columns = {"SiO2_wt%":"SiO2","TiO2_wt%":"TiO2","Al2O3_wt%":"Al2O3","FeOT_wt%":"FeOT","MnO_wt%":"MnO","MgO_wt%":"MgO","CaO_wt%":"CaO",
                                                "Na2O_wt%":"Na2O","K2O_wt%":"K2O","P2O5_wt%":"P2O5","H2O_wt%":"H2O"})
            if redox=='Y':
                random_comps_vc['Fe2O3'] = (random_comps_vc['FeOT']*random_comps_vc['Fe3_FeT'])*(159.69/71.844)
                random_comps_vc['FeO'] = random_comps_vc['FeOT']*(1. - random_comps_vc['Fe3_FeT'])
                random_comps_vc = random_comps_vc.drop(columns = 'FeOT')
            random_comps_vc['CO2'] = random_comps['CO2_ppm']/10000.
            random_comps_vc.to_csv('random_comps.csv', index=True, header=True)
            myfile = vc.BatchFile('random_comps.csv')
            norm_data = myfile.get_data(normalization=normalisation)
            fixed_norm_data = vc.BatchFile(filename=None, dataframe=norm_data)
            fluid_comp = fixed_norm_data.calculate_equilibrium_fluid_comp(temperature='T_C',model=model)
            pvsat_av = fluid_comp['SaturationP_bars_VESIcal'].mean()
            pvsat_sd = fluid_comp['SaturationP_bars_VESIcal'].std()
            XH2O_av = fluid_comp['XH2O_fl_VESIcal'].mean()
            XH2O_sd = fluid_comp['XH2O_fl_VESIcal'].std()
            XCO2_av = fluid_comp['XCO2_fl_VESIcal'].mean()
            XCO2_sd = fluid_comp['XCO2_fl_VESIcal'].std()
            results1 = pd.DataFrame([[run,model,pvsat_av,pvsat_sd,XH2O_av,XH2O_sd,XCO2_av,XCO2_sd]])
            if run == 0:
                results = results1
            else:
                results = pd.concat([results, results1])
            #if run % 50 == 0:
                results.to_csv('results_pvsat_MC_VESIcal_interim'+model+normalisation+redox+'.csv', index=False, header=True)
            tqdmsteps.update(1)
    
    results.columns = ['run','model option','P_bar model','P_bar sd model',"xgH2O model",'xgH2O sd model',"xgCO2 model",'xgCO2 sd model']
    results = results.set_index('run', drop=True)
    results = results.reset_index(drop=True)
    all_results = pd.concat([setup,results], axis=1)
    all_results.to_csv('results_pvsat_MC_VESIcal'+model+normalisation+redox+'.csv', index=False, header=True)
    return all_results