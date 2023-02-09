#! /usr/bin/env python

''' Acab PYthon PArser by ESS-Bilbao Dic 2021
A python tool to get info from ACAB output fort.6 by 
 Mr. Miguel Magan and Dr. Octavio Gonzalez
'''

# import linecache
import io
import numpy as np
# import array
# import sys
import pandas as pd

# =============================================================================

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

intervals = (
    ('years', 31536000), # 60 * 60 * 24 * 365
    ('months', 2592000), # 60 * 60 * 24 * 30
    ('weeks', 604800),  # 60 * 60 * 24 * 7
    ('days', 86400),    # 60 * 60 * 24
    ('hours', 3600),    # 60 * 60
    ('minutes', 60),
    ('seconds', 1),
    )

def display_time(seconds, granularity=2):
    ''' Change any time expresion in seconds to human readable units \n
    granularity: level of detail of the transformation\n'''
    if not is_number(seconds):
        print('Non-numerical time input')
        return None
    result = []
    for name, count in intervals:
        value = seconds // count
        if value:
            seconds-=value * count
            if value == 1:
                name = name.rstrip('s')
            result.append("{}_{}".format(int(value), name))
    if not result:
        result.append('Shutdown')
    return '+'.join(result[:granularity])

def human_time(dataframe,axis,detail=2):
    '''Changes any time expresions list in a pandas dataframe to human values\n
    axis: must be 0 for index or 1 for columns names\n
    detail: level of detail of the transformation\n'''
    times = []
    if axis == 0:
        times = dataframe.index.copy()
        # for index in dataframe.index:
            # times.append(index)
    elif axis == 1:
        times = dataframe.columns.copy()
        # for index in dataframe.columns:
            # times.append(index)
    else:
        print('Axis must be 0 for index or 1 for columns names')
        return dataframe
    print('Before:',times)
    htimes = [display_time(float(time), detail) for time in times]
    print('Now:',htimes)
    newdataframe = dataframe.copy()
    if axis == 0:
        newdataframe.index = htimes
    elif axis == 1:
        newdataframe.columns = htimes
    return newdataframe

def GetDataTable(table):
    """ Get the heat/activity/dose info from data, an ACAB output table.
    Returns 't' as the times, 'iso' as the isotopes, and 'data' as actual data
    """
    l0 = table[0].split()
    l0 = [x if x not in ['RESTART', 'SHUTDOWN'] else '1' for x in l0]
    l0 = [x.replace('S', '') for x in l0]
    t = [float(x) for x in l0[1:]]
    nt = len(t)
    for i, line in enumerate(table):
        if 'TOTAL' in line[0:7]:
            niso = i
            break
    else:
        print('Fatal error reading ACAB output')
        return None
    iso = np.zeros((niso), dtype='6str')
    data = np.zeros((niso, len(t)), dtype='float64')
    for i, line in enumerate(table[1:niso+1]):
        iso[i] = line[1:7].replace(' ','').capitalize()
        palabras = line.split()
        data[i] = list(palabras[-nt:])
    return {'iso': iso, 'time': t, 'data': data}


def __GetStartDatalin(datalines, key):
    """ returns the line position from list datalines where the data relevant
    for key starts  possible keys are 'DISINTEGRATIONS/SEC' for Bq,
    PHOTONS/CCM/SEC for gamma  emission, etc. """
    timesets = []
    for i, lin in enumerate(datalines):
        if lin[0] == '1':
            palabras = lin.split()
            if palabras[-1] == 'SET':
                if key in datalines[i+4]:
                    timeset = palabras[-3].replace('.', '')
                    timesets.append(int(timeset))

    if len(timesets) == 0:
        print('No relevant data found, nothing to do here...')
        return None
    datalin = []
    srchstr = [str(ts)+'. TIME SET' for ts in timesets]
    for i, lin in enumerate(datalines):
        if not (lin[0] == '1' and any(s in lin for s in srchstr)):
            continue
        if key not in datalines[i+4]:
            continue
        datalin.append(i+7)
    return datalin


def __GetVolume(datalines):
    """ Get the volume from datalines. Assumes single value. Will probably
    fail for multiple zones but who uses that"""
    data = iter(datalines)
    for line in data:
        if 'BLOCK #2,  CARD #1' in line:
            line = next(data)
            volume = float(line.split()[0])
            break
    return volume

def Volume_cell(entrada):
    """ Get the volume from file"""
    with open(entrada, 'r', encoding="UTF-8") as datafile:
        line = datafile.readlines()
    volume = __GetVolume(line)
    return volume


def __getNOGG(lines):
    """Internal function to get the number of gamma groups from text object"""
    for line in lines:
        if 'NOGG   NUMBER OF ACTIVATION GAMMA GROUPS' in line:
            NOGG = line.split()[-1]
            break
    else:
        return None  #If somehow, we don't find the string (prolly it is not a ACAB file)
    return int(NOGG)


def __getEGRP(lines):
    """Internal function to get the gamma energies from text object"""
    EGRP = []
    for i, line in enumerate(lines):
        if 'BLOCK #2,  CARD #6' in line:
            nline = i
    data = iter(lines[nline+1:])
    line = next(data)
    while line.strip():
        values = line.split()
        for v in values:
            EGRP.append(float(v))
        line = next(data)
    return EGRP


def getEGRP(entrada):
    """ Get the Activation gamma groups from file"""
    with open(entrada, 'r', encoding="UTF-8") as datafile:
        line = datafile.readlines()
    EGRP = __getEGRP(line)
    return EGRP


def GetTimeSets(entrada):
    ''' Conocer el identificador de los time set de un fort.6 de ACAB '''
    bl7 = np.zeros((0), dtype='8i')
    bl8 = []
    timesets = []
    with open(entrada, 'r', encoding="UTF-8") as datafile:
        lines = datafile.readlines()
    lines = iter(lines)
    for line in lines:
        if 'BLOCK #7,  CARD #1' in line:
            bl7_1 = np.zeros((8), dtype='i')
            line = next(lines)
            bl7_1[0:3] = [int(l) for l in line.split()]
            line = next(lines)
            bl7_1[3:6] = [int(l) for l in line.split()]
            line = next(lines)
            bl7_1[6:8] = [int(l) for l in line.split()]
            bl7 = np.append(bl7, [bl7_1], axis=0)

        if 'BLOCK #8,  CARD #1' in line:
            bl8_1 = np.zeros((bl7[-1, 1]), dtype='f')
            i = 0
            while any(str.strip(line)):
                line = next(lines)
                t = [float(l) for l in line.split()]
                bl8_1[i:i+len(t)] = t
                i += len(t)
            bl8.append(bl8_1)
        if 'BLOCK #13,  CARD #2' in line:
            bl13_2 = np.zeros(len(bl8), dtype=bool)
            i = 0
            while any(str.strip(line)):
                line = next(lines)
                t = [bool(int(l)) for l in line.split()]
                bl13_2[i:i+len(t)] = t
                i += len(t)

    for i, timesteps in enumerate(bl8):
        if bl13_2[i]:
            timesets.append(timesteps)

    return timesets

def purge_dataframe(indata, isotope_list, Threshold):
    '''Keep just the isotopes of interest.
       If isotopes = 'All' just eliminates isotopes that are 0 at every time '''
# One axis is times other is strigs, it should work over strings.
    if all(isinstance(x, str) for x in indata.columns) == True:
        data = indata
    else:
        data = indata.T
#    print(data.columns)
    if 0 >= Threshold or Threshold > 1:
        return None
    if 'Subtot' in data.columns:
        data = data.drop('Subtot',axis= 1)
# These loops list the main isotopes of any time to generata a list of item common fot every index
    main_cols = []
    for i, col in enumerate(data.T.columns):
        part_sum = 0
        j = 1
        stop_count = data.T.sort_values(col, ascending=False)[col][0]*Threshold
        while part_sum < stop_count:
#            print(col, i, j, len(data.columns))
            part_sum += data.T.sort_values(col, ascending=False)[col][j]
            if data.T.sort_values(col, ascending=False).index[j] not in main_cols:
                main_cols.append(data.T.sort_values(col, ascending=False).index[j])
            j += 1
            if j == len(data.columns):
                break
#            print(i,data.T.sort_values(col,ascending=False).index[i],part_sum,stop_count)
#    print(main_cols,len(main_cols))
# Let's start the purge...
    for colum in  data.columns:
        if colum not in isotope_list and colum != 'Total' and colum not in main_cols:
            data = data.drop(columns = colum)
    return data

# =========================================================================
# =========================================================================

def RadActIsotopes_full_pd(entrada, isotope_list='All', Threshold=0.9):
    ''' Gets the Radionulcide decay of a cell from ACAB\n
    isotope_list: List of isotopes of interest (All: check all of them)\n
    Threshold: Pay attentions to isotopes up to represent Total*Threshold at any time'''
    print('RadActIsotopes_full_pd',entrada, isotope_list, Threshold)
    with open(entrada, 'r', encoding="UTF-8") as datafile:
        line = datafile.readlines()
    datalin = __GetStartDatalin(line, 'DISINTEGRATIONS/SEC')
#    print(datalin)
    if not datalin:
        print('NO radioactive material, nothing to do here...\n')
        return 0, 0
    raw_data = [GetDataTable(line[dl:]) for dl in datalin]
# the same using pandas! sort this a bit, put isotopes together & eliminate redundant RESTARTS
    data = pd.DataFrame(raw_data[0]['data'],raw_data[0]['iso'],raw_data[0]['time'])
    additionalframes = []
    for raw_i in raw_data[1:]:
        additionalframes.append(pd.DataFrame(raw_i['data'],raw_i['iso'],raw_i['time']))
    for frame in additionalframes:
        for colum in frame.columns:
            if colum not in data.columns:
                datacolum=frame[colum]
                data = data.join(datacolum)
    dataT = data.T
    dataT.set_index = 'decay_times_s'
# Fun starts here: Keep just the isotopes of interest if isotopes = 'All' just
## eliminates isotopes that are 0 at everytime
    dataT = purge_dataframe(dataT,isotope_list,Threshold)
    if not isinstance(dataT,pd.DataFrame):
        print('Wrong threshold limit, it should be between 0 and 1')
        return None
#    Tool to generate molar files for interest nuclides
    print("Easy print: \n dataT.to_csv('Summary_decay_Nuclides_Bq_s.csv',sep=',', index_label = 'decay_times_s',float_format='%.4E')")

    print("Easy plot :\n data.plot(y=['isotope'],xlabel='decay_times_s',ylabel='Nuclide radioactivity Bq/s',logx=True,logy=True,colormap='jet').get_figure().savefig('Summary_decay_Nuclides_Bq_s.png',bbox_inches='tight',dpi=200)")

    return dataT

#===================================================================================================

def Isotopes_molarity(entrada,isotope_list='All',Threshold=0.9):
    ''' Gets the mol composition of a cell from ACAB\n
    isotope_list: List of isotopes of interest (All: check all of them)\n
    Threshold: Pay attentions to isotopes up to represent Total*Threshold at any time'''
    print('Isotopes_molarity ',entrada, isotope_list, Threshold)
    with open(entrada, 'r', encoding="UTF-8") as datafile:
        line = datafile.readlines()
    datalin = __GetStartDatalin(line, 'NUCLIDE CONCENTRATIONS,  AT.GR.')
#    print(datalin)

    if not datalin:
        print('No concentration results, nothing to do here...')
        return 0, 0

    raw_data = [GetDataTable(line[dl:]) for dl in datalin]

# the same using pandas! sort this a bit, put isotopes together & eliminate redundant RESTARTS
    data = pd.DataFrame(raw_data[0]['data'],raw_data[0]['iso'],raw_data[0]['time'])
    additionalframes = []
    for raw_i in raw_data[1:]:
        additionalframes.append(pd.DataFrame(raw_i['data'],raw_i['iso'],raw_i['time']))
    for frame in additionalframes:
        for colum in frame.columns:
            if colum not in data.columns:
                datacolum=frame[colum]
                data = data.join(datacolum)
    dataT = data.T
    dataT.set_index = 'decay_times_s'
# Fun starts here: Keep just the isotopes of interest if isotopes = 'All' just
## eliminates isotopes that are 0 at everytime
    dataT = purge_dataframe(dataT,isotope_list,Threshold)
    if not isinstance(dataT,pd.DataFrame):
        print('Wrong threshold limit, it should be between 0 and 1')
        return None

    print("Easy print: \n dataframe.to_csv('Nuclide_quantities_in_Mol.csv',sep=',', index_label = 'decay_times_s',float_format='%.4E')")

    print("Easy plot :\n data.plot(y=['isotope'],xlabel='decay_times_s',ylabel='mol',logx=True,logy=True,colormap='jet').get_figure().savefig('prueba.png',bbox_inches='tight',dpi=200)")

    return dataT

# =============================================================================
# =============================================================================

def gammas_full_pd(entrada):
    ''' Permite conocer la emision gamma a lo largo de todo el tiempo de decaimiento '''
    print('gammas_full_pd', entrada)
    with open(entrada, 'r', encoding="UTF-8") as datafile:
        line = datafile.readlines()
    NOGG = __getNOGG(line)
    datalin = __GetStartDatalin(line, 'PHOTONS/CCM/SEC')

    if not datalin:
        print("NO Gamma data found, nothing to do here...\n")
        return None
    datalin = [dl + 1 for dl in datalin]  # due to fort.6 format
    raw_data = [line[dl:dl+NOGG+1] for dl in datalin]

    data = pd.read_csv(io.StringIO('\n'.join(raw_data[0]).replace('RESTART', '1.0').
                                   replace('S', '')), delim_whitespace=True)
    add_data = [pd.read_csv(io.StringIO('\n'.join(raw).replace('RESTART', '1.0').replace('S', '')),
                          delim_whitespace=True) for raw in raw_data[1:]]

    for tab in add_data:
        for col in tab.columns:
            if col not in data.columns:
                datacolum = tab[col]
                data = data.join(datacolum)
    data = data.set_index('(MEV)')

    total = [data[time].sum() for time in data.columns]
    total_df = pd.DataFrame([total], columns=data.columns, index=['Total'])
    data = data.append(total_df)

    print("Easy print: \n data.to_csv('Gamma_release_rates_Gammas_cm2_s.csv',sep=',',index_label = 'Energy_MeV',float_format='%.4E')")
    return data

def gammas_by_Egroups(data, EGROUPS, plot=False):
    ''' Distribute gamma data by energy groups list '''
    dataT = data.T
    total = list(dataT['total'])
    dataT.drop(['total'],axis = 1, inplace = True)
    data = dataT.T

    data_by_Egroups = pd.DataFrame(np.zeros((len(data.columns),len(EGROUPS))),
                                   columns=EGROUPS,index=data.columns)
    for i,j in np.ndindex(len(data.columns),len(data.index)):
        if data.index[j] >= EGROUPS[1]:
            data_by_Egroups[EGROUPS[2]][data.columns[i]] = data_by_Egroups[EGROUPS[2]][data.columns[i]] + data[data.columns[i]][data.index[j]] / total[i]
        if data.index[j] < EGROUPS[0]:
            data_by_Egroups[EGROUPS[0]][data.columns[i]] = data_by_Egroups[EGROUPS[0]][data.columns[i]] + data[data.columns[i]][data.index[j]] / total[i]
        if EGROUPS[0] <= data.index[j] < EGROUPS[1]:
            data_by_Egroups[EGROUPS[1]][data.columns[i]] = data_by_Egroups[EGROUPS[1]][data.columns[i]] + data[data.columns[i]][data.index[j]] / total[i]

    Index_Egroups = []
    for i in EGROUPS:
        Index_Egroups.append('<=_'+str(i)+'MeV')
    data_by_Egroups.columns = Index_Egroups

    data_by_EgroupsT = data_by_Egroups.T
    print("Easy print: \n data_by_EgroupsT.to_csv('Gamma_release_rates_Gammas_cm3_s_byE_groups.csv',sep=',', index_label = 'Energy_MeV',float_format='%.4E')")

#dataE.T.plot.bar(xlabel='decay_times_s',stacked=True,colormap='jet').get_figure().savefig('prueba.png',bbox_inches='tight',dpi=200)

    if plot is True:
        dataE = pd.DataFrame(np.zeros((len(data.columns),len(EGROUPS))),
                                   columns=EGROUPS,index=data.columns)
        for i,j in np.ndindex(len(data.columns),len(data.index)):
            if data.index[j] >= EGROUPS[1]:
                dataE[EGROUPS[2]][data.columns[i]] = dataE[EGROUPS[2]][data.columns[i]] + data[data.columns[i]][data.index[j]]
            if data.index[j] < EGROUPS[0]:
                dataE[EGROUPS[0]][data.columns[i]] = dataE[EGROUPS[0]][data.columns[i]] + data[data.columns[i]][data.index[j]]
            if EGROUPS[0] <= data.index[j] < EGROUPS[1]:
                dataE[EGROUPS[1]][data.columns[i]] = dataE[EGROUPS[1]][data.columns[i]] + data[data.columns[i]][data.index[j]]
        if all(is_number(x) for x in dataE.index):
            dataE = human_time(dataE,0,2)
        dataE.columns = Index_Egroups
        plotname = "Gamma_release_rates_Gammas_cm2_s_byE_groups.png"
        dataE.plot.bar(xlabel='decay_times_s', ylabel='Gamma_release_rates_Gammas/cm3/s',
                       stacked=True, colormap='jet').get_figure().savefig(plotname, bbox_inches='tight', dpi=200, figsize=(16, 9))
    print("Easy plot :use kwarg plot=True")
    return data_by_EgroupsT

#=============================================================================
def HeatIsotopes_full_pd(entrada, isotope_list='All', Threshold=0.9):
    ''' Gets the Radionuclide afterheat generation of a cell from ACAB\n
    isotope_list: List of isotopes of interest (All: check all of them)\n
    Threshold: Pay attention to isotopes up to represent Total*Threshold at any time'''
    print('HeatIsotopes_full',entrada, isotope_list,Threshold)
    with open(entrada, 'r', encoding="UTF-8") as datafile:
        line = datafile.readlines()
    datalin = __GetStartDatalin(line, 'NUCLIDE THERMAL POWER, WATTS')
    volume = __GetVolume(line)
    if not datalin:
        print('No afterheat results, nothing to do here...')
        return None
    raw_data = [GetDataTable(line[dl:]) for dl in datalin]
# the same using pandas! sort this a bit, put isotopes together & eliminate redundant RESTARTS 
    data = pd.DataFrame(raw_data[0]['data'],raw_data[0]['iso'],raw_data[0]['time'])
    additionalframes = [pd.DataFrame(r['data'], r['iso'], r['time']) for r in raw_data]
    for frame in additionalframes:
        for colum in frame.columns:
            if colum not in data.columns:
                datacolum = frame[colum]
                data = data.join(datacolum)
    dataT = data.T/volume # Change of units to w/cm3!!!
# Fun starts here: Keep just the isotopes of interest if isotopes = 'All' just 
## eliminates isotopes that are 0 at everytime
    dataT = purge_dataframe(dataT, isotope_list, Threshold)
    if not isinstance(dataT, pd.DataFrame):
        print('Wrong threshold limit, it should be between 0 and 1')
        return None

#    Tool to generate molar files for interest nuclides
    print("Easy print: \n dataT.to_csv('Summary_afterheat_Watts_cm3.csv',sep=',', index_label = 'decay_times_s',float_format='%.4E')")
    print("Easy plot :\n data.plot(y=['isotope'],xlabel='decay_times_s',ylabel='Nuclide afterheat W/cm3',logx=True,logy=True,colormap='jet').get_figure().savefig('prueba.png',bbox_inches='tight',dpi=200)")

    return dataT
# =========================================================================

def GammaDoseIsotopes_full_pd(entrada, isotope_list='All', Threshold=0.9):
    ''' Gets the Radionuclide gamma generation of a cell from ACAB\n
    isotope_list: List of isotopes of interest (All: check all of them)\n
    Threshold: Pay attentions to isotopes up to represent Total*Threshold at any time'''
    print('GammaDoseIsotopes_full ',entrada, isotope_list, Threshold)
    with open(entrada, 'r', encoding="UTF-8") as datafile:
        line = datafile.readlines()
    datalin = __GetStartDatalin(line, 'SURFACE GAMMA DOSE RATES DUE TO')
#    print(datalin)
    if not datalin:
        print('No Dose results, nothing to do here...')
        return None
    raw_data = [GetDataTable(line[dl:]) for dl in datalin]
# the same using pandas! sort this a bit, put isotopes together & eliminate redundant RESTARTS
    data = pd.DataFrame(raw_data[0]['data'],raw_data[0]['iso'],raw_data[0]['time'])
    additionalframes = []
    for raw_i in raw_data[1:]:
        additionalframes.append(pd.DataFrame(raw_i['data'],raw_i['iso'],raw_i['time']))
    for frame in additionalframes:
        for colum in frame.columns:
            if colum not in data.columns:
                datacolum = frame[colum]
                data = data.join(datacolum)
    dataT = data.T*1000 # Change of units to mSv/h!!!!
# Fun starts here: Keep just the isotopes of interest if isotopes = 'All' just
## eliminates isotopes that are 0 at everytime
    dataT = purge_dataframe(dataT, isotope_list, Threshold)
    if not isinstance(dataT, pd.DataFrame):
        print('Wrong threshold limit, it should be between 0 and 1')
        return None

    print("Easy print: \n dataT.to_csv('Summary_gamma_Eq_dose_mSv_h.csv',sep=',', index_label = 'decay_times_s',float_format='%.4E')")
    print("Easy plot :\n data.plot(y=['isotope'],xlabel='decay_times_s',ylabel='Surface equivalent gamma dose rates mSv/h',logx=True,logy=True,colormap='jet').get_figure().savefig('Summary_gamma_Eq_dose_mSv_h.png',bbox_inches='tight',dpi=200)")

    return dataT

# =============================================================================
# ============================================================================
