#!/suzakuproc/data1/conda/anaconda3/bin/python3

import os
#import lines
from sklearn import cluster
from subprocess import check_output
#from string import find
import numpy as np
#from bokeh.charts import Histogram,Scatter
#from bokeh.plotting import figure,output_file,show
from pandas import DataFrame
from astropy.coordinates import SkyCoord
import astropy.units as u
import pandas as pd
from collections import Counter
from astropy.table import Table
from astropy.io import ascii
from astroquery.simbad import Simbad
import json


# read dictionary with SNRcat info for dictance and age

f = open("snr_dict.txt","r")
snr_dict = json.loads(f.read())
f.close()

with open("make_snr_catalog.config","r") as f:
    config = json.loads(f.read())
#f.close()


# reading pulsar dictionary

f = open("snr_psr_dict.txt","r")
snr_psr_dict = json.loads(f.read())
f.close()


#data = lines.linedata()
#mysb = Simbad()
#mysb.add_votable_fields('distance')


# get LMC SMC coordinates

#smc_tab = mysb.query_object("smc")
#lmc_tab = mysb.query_object("lmc")

#smc_tab.write("smc_tab.txt",format='ascii')
#lmc_tab.write("lmc_tab.txt",format='ascii')
smc_tab = ascii.read("smc_tab.txt")
lmc_tab = ascii.read("lmc_tab.txt")

#quit()

smc = SkyCoord(ra=smc_tab[0]["RA"],dec=smc_tab[0]["DEC"],unit=(u.hourangle,u.deg),frame='fk5')
lmc = SkyCoord(ra=lmc_tab[0]["RA"],dec=lmc_tab[0]["DEC"],unit=(u.hourangle,u.deg),frame='fk5')

#print(smc,lmc)

master = ascii.read("master.dat",delimiter="|")
#print(master.head())

#snr_sources = pd.read_pickle("snrs_sources.pkl")
snr_sources = pd.read_excel("snrs_m012417.xlsx")
snr_sources["TYPE"].fillna("NA")

#snr_obsids = pd.read_pickle("snrs_obsids.pkl")
snr_obsids = pd.read_excel("snrs_obsids.xlsx")

snr_coord = SkyCoord(ra=snr_sources["RA"].values,dec=snr_sources["DEC"].values,unit=(u.deg,u.deg),frame='fk5')


#print snr_sources

#print snr_sources
#quit()
#data_ids = lines.id_list()
#print data


xarray = []

ew = []
ew_max = []
ew_min = []
linee = []
linee_max = []
linee_min = []
sigma = []
sigma_max = []
sigma_min = []
norm = []
norm_max = []
norm_min = []
sarr = []
snrarr = []
spec_arr = []
stype = []
#iddict = {"si_he":0,"s_he":1,"mg_he":2,"o_he":3,"o_h":4,"ne_he":5,"ne_h":6,"mg_he":7}
ew_arr = []
obsids = []
pl_model = []
pl_index = []
pl_minind = []
pl_maxind = []
nh_model = []
nh_arr = []
nh_minind = []
nh_maxind = []
inst_arr = []
#str_ids = []
ra_obs = []
dec_obs = []
ra_snr = []
dec_snr = []
sep_arr = []
d_min = []
d_max = []
a_min = []
a_max = []
psr_arr = []
psr_age_arr = []
dtd_arr = []
region_arr = []
image_arr = []

# arrays to hold fake line records
ew_0 = []
ew_max_0 = []
ew_min_0 = []
linee_0 = []
linee_max_0 = []
linee_min_0 = []
sigma_0 = []
sigma_max_0 = []
sigma_min_0 = []
norm_0 = []
norm_max_0 = []
norm_min_0 = []
sarr_0 = []
snrarr_0 = []
spec_arr_0 = []
stype_0 = []
#iddict = {"si_he":0,"s_he":1,"mg_he":2,"o_he":3,"o_h":4,"ne_he":5,"ne_h":6,"mg_he":7}
ew_arr_0 = []
obsids_0 = []
pl_model_0 = []
pl_index_0 = []
pl_minind_0 = []
pl_maxind_0 = []
nh_model_0 = []
nh_arr_0 = []
nh_minind_0 = []
nh_maxind_0 = []
inst_arr_0 = []
ra_obs_0 = []
dec_obs_0 = []
ra_snr_0 = []
dec_snr_0 = []
sep_arr_0 = []
d_min_0 = []
d_max_0 = []
a_min_0 = []
a_max_0 = []
psr_arr_0 = []
psr_age_arr_0 = []
dtd_arr_0 = []
region_arr_0 = []
image_arr_0 = []

obsids_unique = []
data_ind = 0
spec_number = 0
line_count = 0
line_triaged = 0

json_file_list = os.listdir(config['json_dir'])

#for fit in data:
for json_file_name in json_file_list:

    try:
        with open(os.path.join(config['json_dir'],json_file_name),'r') as fjson:
            fit = json.loads(fjson.read())
    except Exception as e:
        print("Can't read file {0}".format(json_file_name))
        continue
        
#    source = fit["source"]
#    obsdata = data[fit]
    oid = fit["obsid"]
#    print(fit)
#    id_data = data_ids[data_ind]
    data_ind = data_ind + 1
    obsids_unique.append(oid)
    try:
        image_file = fit["imagefile"]
    except:
        image_file = "none"
    print(image_file)
    out = check_output(["grep",'|'+oid+'|',"master.txt"])
    spl_out = out.decode("utf-8").split("|")
#    print spl_out
#    i = [ i for i in range(len(mas_ids)) if mas_ids[i] == id]
#    print i,mas_ids[i]
    source = spl_out[1]
    if "source" in fit.keys(): source = fit["source"]
    snr_source = source


#    st = lines.stype_num(source)

    ra = spl_out[2]
    dec = spl_out[3]
    time = spl_out[7]
    coord_obs =  SkyCoord(ra=ra,dec=dec,unit=(u.deg,u.deg),frame='fk5')

    dtd = "L"
    sep = 0.0
    st = "NA"
    try:

        source_names = snr_obsids[snr_obsids["OBSID"] == int(oid)]["SOURCE"].values
        sz =  source_names.shape[0]
#        print sz,oid

        if "snr" in fit.keys():

            snr_source = fit["snr"]
#            print snr_source
            sr = snr_sources[snr_sources["NAME"] == snr_source]
            st = sr["TYPE"].values[0]
            ra_src = sr["RA"].values[0]
            dec_src = sr["DEC"].values[0]
#            print(snr_source,ra_src,dec_src)

        elif sz == 0:

            isnr,d2d,d3d =  coord_obs.match_to_catalog_sky(snr_coord)
#           print isnr,snr_sources["NAME"].values[isnr]
            if d2d.deg < 10.0:
                snr_source = snr_sources["NAME"].values[isnr]
#           source_index = isnr
#           print source_ind.astype(int),d2d.deg,d3d,source
#           source_index = source_ind.astype(int)
#           source_index = snr_sources.iloc[source_index]
 #          print source_index
                ra_src = snr_sources["RA"].values[isnr]
                st = sr["TYPE"].values[0]
                dec_src = snr_sources["DEC"].values[isnr]
            else:
                snr_source = source
                ra_src = ra
                dec_src = dec

#           print snr_source,isnr,sz

        else:

            source_index = source_names[0]
            snr_source = snr_sources["NAME"].iloc[source_index]

#            print snr_source,source_index,sz
#       st = lines.stype_num(snr_source)

            ra_src = snr_sources["RA"].iloc[source_index]
            dec_src = snr_sources["DEC"].iloc[source_index]
            st = snr_sources["TYPE"].iloc[source_index]

#        dist_response = mysb.query_object(snr_source)
#       print dist_response
        dist_min = np.nan
        dist_max = np.nan
        age_min = np.nan
        age_max = np.nan
        psr_name = ""
        psr_sep = np.nan
        psr_age = np.nan


        if snr_source[4] != "G":
#

            c = SkyCoord(ra=ra_src,dec=dec_src,unit=(u.hourangle,u.deg),frame='fk5')
            ra_src = c.ra.degree
            dec_src = c.dec.degree

#           print snr_source,snr_source[0], ra_src, snr_sources["RA"].iloc[source_index], dec_src, snr_sources["DEC"].iloc[source_index]

            lmc_sep = c.separation(lmc).degree
            smc_sep = c.separation(smc).degree
            if lmc_sep > smc_sep:
#               print "SMC",smc_sep,lmc_sep
                dist_min = 60.0
                dist_max = 62.0
            else:

#               print "LMC",smc_sep,lmc_sep
                dist_min = 50.0
                dist_max = 50.0



        else:

            try:
#               print snr_dict[snr_source[4:]],snr_source[4:]
                psr_name = snr_psr_dict[snr_source]["pulsar"]
                psr_sep = snr_psr_dict[snr_source]["sep"]
                psr_age = snr_psr_dict[snr_source]["age"]
                sd = snr_dict[snr_source[4:]]
                dist_min = sd['dmin']
                dist_max = sd['dmax']
                age_min = sd['amin']
                age_max = sd['amax']

            except Exception as e:
                pass
                #print source,oid,e.message

            c =  SkyCoord(ra=ra_src,dec=dec_src,unit=(u.deg,u.deg),frame='fk5')

        sep = coord_obs.separation(c).degree

    except Exception as e:
        print(e)
        pass

#    st = lines.stype_num(snr_source)
#    if "sntype" in fit.keys(): st = fit["sntype"]
#    st = 
    try:
        if np.isnan(st): st = "NA"
    except:
        pass
#    print source,st


    xarr = np.zeros(8,dtype=float)

#    try:

    power_laws = [ c for c in fit["fit"] if fit["fit"][c]["model"] in ["powerlaw","pegpwrlw"]]
    nhs = [ c for c in fit["fit"] if fit["fit"][c]["model"] in ["phabs","wabs"]]

    print(len(power_laws))

    if len(power_laws):
        pl = fit["fit"][power_laws[0]]

#    for comp in fit["fit"].keys():

#        pl = fit['powerlaw']
        pl_m = pl["model"]
        pl_i = pl["params"]["PhoIndex"]
        pl_v = pl_i["value"]
        pl_mn = pl_i["min"]
        pl_mx = pl_i["max"]

    else:

        print(oid," - No PL component")
        pl_m = 0
        pl_v = None
        pl_mn = None
        pl_mx = None

#    try:
    if len(nhs):

        nh = fit["fit"][nhs[0]]
#       pl_m = fit["pl_model"]
        nh_par = nh["params"]["nH"]
        nh_v = nh_par["value"]
        nh_mn = nh_par["min"]
        nh_mx = nh_par["max"]


    else:

        print(oid," - No nH component")
        nh_v = None
        nh_mn = None
        nh_mx = None


    xiss = ["xi"+f["spectrum"].split("_")[0][-1] for f in fit["data"]]
    regcode = fit["data"][0]["spectrum"].split("_")[1]
    inst = " ".join(xiss)
    ga_count = 0

    for comp in fit["fit"].keys():

        if comp in power_laws or comp in nhs: continue
#       print comp

        try:
            c = fit["fit"][comp]
            pars = c["params"]
            if c["model"] != "gaussian": continue

#   rejecting bad lines 

            if pars["LineE"]["value"] > 10 or pars["LineE"]["value"] < 0.35 or\
                abs(c["eqwidth"]["value"]) > 10.0 or \
                pars["LineE"]["value"]/pars["Sigma"]["value"] < 7.5 or \
                pars["norm"]["value"] < 1.0e-8:
                line_triaged += 1
                continue

            xarray.append((pars["LineE"]["value"],pars["Sigma"]["value"]))
            stype.append(st)
            spec_arr.append(spec_number)
            inst_arr.append(inst)
            ew_arr.append(c["eqwidth"]["value"])
            obsids.append(oid)
            sarr.append(source)
            snrarr.append(snr_source)
            pl_model.append(pl_m)
            pl_index.append(pl_v)
            pl_minind.append(pl_mn)
            pl_maxind.append(pl_mx)
            nh_arr.append(nh_v)
            nh_minind.append(nh_mn)
            nh_maxind.append(nh_mx)
            d_min.append(dist_min)
            d_max.append(dist_max)
            a_min.append(age_min)
            a_max.append(age_max)
            psr_arr.append(psr_name)
            psr_age_arr.append(psr_age)
            region_arr.append(regcode)
            image_arr.append(image_file)
            ra_obs.append(ra)
            dec_obs.append(dec)
            ra_snr.append(ra_src)
            dec_snr.append(dec_src)
            sep_arr.append(sep)

            linee.append(pars["LineE"]["value"])
            linee_max.append(pars["LineE"]["max"])
            linee_min.append(pars["LineE"]["min"])
            ew.append(c["eqwidth"]["value"])
            ew_max.append(c["eqwidth"]["max"])
            ew_min.append(c["eqwidth"]["min"])
            sigma.append(pars["Sigma"]["value"])
            sigma_max.append(pars["Sigma"]["max"])
            sigma_min.append(pars["Sigma"]["min"])
            norm.append(pars["norm"]["value"])
            norm_max.append(pars["norm"]["max"])
            norm_min.append(pars["norm"]["min"])
            dtd_arr.append("L")
            ga_count = ga_count + 1


        except Exception as e:
            print(e)
            pass

#    spec_number = spec_number + 1

    if ga_count == 0:


        stype_0.append(st)
        spec_arr_0.append(spec_number)
        inst_arr_0.append(inst)
        ew_arr_0.append(0.0)
        obsids_0.append(oid)
        sarr_0.append(source)
        snrarr_0.append(snr_source)
        pl_model_0.append(pl_m)
        pl_index_0.append(pl_v)
        pl_minind_0.append(pl_mn)
        pl_maxind_0.append(pl_mx)
        nh_arr_0.append(nh_v)
        nh_minind_0.append(nh_mn)
        nh_maxind_0.append(nh_mx)

        ra_obs_0.append(ra)
        dec_obs_0.append(dec)
        ra_snr_0.append(ra_src)
        dec_snr_0.append(dec_src)
        sep_arr_0.append(sep)
        psr_arr_0.append(psr_name)
        psr_age_arr_0.append(psr_age)

        linee_0.append(0.0)
        linee_max_0.append(0.0)
        linee_min_0.append(0.0)
        ew_0.append(0.0)
        ew_max_0.append(0.0)
        ew_min_0.append(0.0)
        sigma_0.append(0.0)
        sigma_max_0.append(0.0)
        sigma_min_0.append(0.0)
        norm_0.append(0.0)
        norm_max_0.append(0.0)
        norm_min_0.append(0.0)
        d_min_0.append(dist_min)
        d_max_0.append(dist_max)
        a_min_0.append(age_min)
        a_max_0.append(age_max)
        image_arr_0.append(image_file)

        region_arr_0.append(regcode)

        if pl_v > 0:
            dtd_arr_0.append("DT")
        else:
            dtd_arr_0.append("N")

        print(oid,dtd_arr_0[-1])

    spec_number = spec_number + 1



print("Triaged lines:",line_triaged)


ions = ['','fe_20','fe_17','fe_21','o_7','o_8','fe_20','fe_17','fe_19','ar_16',
        "ni_20","fe_21","ne_9","fe_23",'fe_22','ne_10','fe_21','fe_22','fe_23','fe_24',
        'ne_10','mg_11','mg_12','fe_22','mg_11',"al_12","mg_11",'al_13','fe_23/24','si_13',
        'si_13','si_14','si_13','si_13','si_14','s_15','s_16','s_15','ar_17','ar_18','ar_16','ca_19','ca_20',
       'ca_19','fe_24','cr_22','mn_23','fe_18','fe_24','fe_26','co_25',
       'ni_26','ni_26','fe_26','fe_22/25','NA']
transition = ['','19-287','3-143','22-455','1-7','1-4','14-42','1-2','1-11','3-20',
              "3-29","7-28","1-7","9-15",'6-23','1-4','1-40','1-20','5-20','3-8',
              '1-7','1-7','1-4','1-270','1-17','1-7',"1-31",'1-4','1-161/18','1-2',
              '1-7','1-4','1-17','1-49','1-7','1-7','1-4',
             '1-17','1-7','1-4','1-85','1-7','1-4','1-17','7-67','1-69',
             '1-30','1-337','1-69','1-4','2-25','16-151','1-69','1-7','2-blend','']
energies = [0.0,0.357,0.421,0.494,0.574,0.653,0.707,0.725,0.821,0.778,
            0.856,0.884,0.921,0.979,0.971,1.021,1.009,1.053,1.056,1.109,
            1.21,1.351,1.472,1.531,1.579,1.597,1.65,1.728,1.728,1.838,
            1.864,2.005,2.181,2.344,2.375,2.459,2.621,
           2.882,3.137,3.321,3.619,3.899,4.104,4.579,5.451,5.656,
           6.151,6.43,6.657,6.968,7.096,7.517,7.758,8.246,1.15,0.0]

init_arr = np.array([[0.0,0.01],[0.4,0.01],[0.42,0.02],[0.49,0.01],[0.56,0.001],[0.66,0.001],[0.68,0.02],[0.72,0.001],
                     [0.82,0.001],[0.84,0.01],[0.86,0.01],[0.88,0.01],[0.91,0.001],[0.96,0.01],[0.98,0.02],[1.01,0.001],[1.01,0.02],[1.03,0.02],[1.05,0.001],[1.07,0.001],[1.2,0.01],
                     [1.33,0.005],[1.47,0.01],[1.5,0.05],[1.55,0.01],[1.6,0.01],[1.65,0.01],[1.71,0.01],[1.73,0.01],[1.85,0.001],
                     [1.86,0.02],[2.0,0.01],[2.18,0.01],[2.33,0.01],[2.32,0.1],[2.44,0.01],
                     [2.6,0.01],[2.88,0.01],[3.11,0.01],[3.4,0.01],[3.6,0.01],[3.9,0.01],
                     [4.12,0.01],[4.5,0.01],[5.46,0.01],[5.56,0.01],[6.1,0.01],[6.43,0.03],[6.7,0.2],
                     [6.9,0.01],[7.1,0.01],[7.5,0.1],[7.8,0.1],[8.2,0.01],[1.09,0.02]])

#init_arr = np.array([[0.42,0.01],[0.49,0.01],[0.56,0.001],[0.66,0.001],[0.72,0.001],[0.8,0.001],[0.92,0.001],[1.01,0.001],[1.05,0.001],[1.2,0.01],
#        [1.33,0.005],[1.47,0.01],[1.5,0.05],[1.6,0.01],[1.73,0.01],[1.85,0.001],[1.86,0.02],[2.0,0.01],[2.18,0.01],[2.33,0.01],[2.32,0.1],[2.44,0.01],
#        [2.6,0.01],[2.88,0.01],[3.11,0.01],[3.4,0.01],[3.6,0.01],[3.9,0.01],[4.12,0.01],[4.5,0.01],[5.46,0.01],[6.1,0.01],
#        [6.43,0.03],[6.7,0.2],[6.9,0.01],[7.1,0.01],[7.5,0.1],[7.8,0.1],[8.2,0.01],[1.09,0.02]])

ncluster = len(init_arr)
#print "Numberof clusteors:"+ncluster

k_means = cluster.KMeans(ncluster,init=init_arr)
k_means.fit(xarray)

#print(k_means.cluster_centers_)
#labels_tmp = [ 15 if i == 16 else i for i in k_means.labels_ ]
#labels = [ 8 if i == 39 else i for i in labels_tmp ]

labels = k_means.labels_
count_dict = Counter(labels)

counts = np.zeros(ncluster,dtype=int)
centers = k_means.cluster_centers_[:,0]
centers_sigmas = k_means.cluster_centers_[:,1]

for tup in count_dict.most_common(ncluster):
    i = tup[0]
    cnt = count_dict[i]
    print("{:2d} {:.3f}  {:.3f}  {:3d}".format(i,centers[i],centers_sigmas[i],cnt))


str_ids = ["{} ({})".format(ions[lab],transition[lab]) for lab in labels]


#print len(ra_snr),len(dec_snr)
#print sep_arr_0,sep_arr
line_dict = {
    "obsid":obsids,
    "spectrum":spec_arr,
    "linee":linee,
    "linee_max":linee_max,
    "linee_min":linee_min,
    "sigma":sigma,
    "sigma_max":sigma_max,
    "sigma_min":sigma_min,
    "norm":norm,
    "norm_max":norm_max,
    "norm_min":norm_min,
    "norm":norm,
    "norm_max":norm_max,
    "norm_min":norm_min,
    "ew":ew,
    "ew_max":ew_max,
    "ew_min":ew_min,
    "pl_model":pl_model,
    "pl_index":pl_index,
    "pl_indmin":pl_minind,
    "pl_indmax":pl_maxind,
    "nh":nh_arr,
    "nh_min":nh_minind,
    "nh_max":nh_maxind,
    "label":labels,
    "id":str_ids,
    "type":stype,
    "source":sarr,
    "snr":snrarr,
    "instrument":inst_arr,
    "ra_obs":ra_obs,
    "dec_obs":dec_obs,
    "ra_snr":ra_snr,
    "dec_snr":dec_snr,
    "sep":sep_arr,
    "dist_min":d_min,
    "dist_max":d_max,
    "age_min":a_min,
    "age_max":a_max,
    "pulsar":psr_arr,
    "psr_age":psr_age_arr,
    "detect":dtd_arr,
    "regcode":region_arr,
    "imagefile":image_arr
}

str_ids_0 = ["" for label in range(len(obsids_0))]

line_dict_0 = {
    "obsid":obsids_0,
    "spectrum":spec_arr_0,
    "linee":linee_0,
    "linee_max":linee_max_0,
    "linee_min":linee_min_0,
    "sigma":sigma_0,
    "sigma_max":sigma_max_0,
    "sigma_min":sigma_min_0,
    "norm":norm_0,
    "norm_max":norm_max_0,
    "norm_min":norm_min_0,
    "norm":norm_0,
    "norm_max":norm_max_0,
    "norm_min":norm_min_0,
    "ew":ew_0,
    "ew_max":ew_max_0,
    "ew_min":ew_min_0,
    "pl_model":pl_model_0,
    "pl_index":pl_index_0,
    "pl_indmin":pl_minind_0,
    "pl_indmax":pl_maxind_0,
    "nh":nh_arr_0,
    "nh_min":nh_minind_0,
    "nh_max":nh_maxind_0,
    "label":-1*np.ones(len(obsids_0),dtype=int),
    "id":str_ids_0,
    "type":stype_0,
    "source":sarr_0,
    "snr":snrarr_0,
    "instrument":inst_arr_0,
    "ra_obs":ra_obs_0,
    "dec_obs":dec_obs_0,
    "ra_snr":ra_snr_0,
    "dec_snr":dec_snr_0,
    "sep":sep_arr_0,
    "dist_min":d_min_0,
    "dist_max":d_max_0,
    "age_min":a_min_0,
    "age_max":a_max_0,
    "pulsar":psr_arr_0,
    "psr_age":psr_age_arr_0,
    "detect":dtd_arr_0,
    "regcode":region_arr_0,
    "imagefile":image_arr_0
}

#print line_dict
line_data_1 = DataFrame(data=line_dict)
line_data_0 = DataFrame(data=line_dict_0)
line_data = pd.concat([line_data_1,line_data_0])
#line_data.loc[np.isnan(line_data['type']),'type'] = "NA"
line_data.to_pickle("learned_lines.pkl")
line_data.to_excel("learned_lines.xlsx")

fits_dat = Table.from_pandas(line_data)
fits_dat.write("lines_list_id.fits",format='fits')

print(stype_0,obsids_0)
#print line_data
