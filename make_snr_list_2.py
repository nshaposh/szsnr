#!/suzakuproc/data1/conda/anaconda3/bin/python3

#import aplpy
import os
from numpy import isnan
import subprocess as sp
#import pyregion
import sys
import astropy.io.fits as pyfits
from math import isnan
#from astroquery.simbad import Simbad
from astropy import units as u
from astropy.coordinates import SkyCoord,FK5
from astropy.table import Table
from astropy.io import ascii
import pandas as pd
#import lines
import copy
import json
from collections import Counter
#obsid_list = lines.get_obsid_list()
#print obsid_list


with open('make_snr_catalog.config','r') as f:
    config = json.loads(f.read())

spec_list = open('spec_list.txt','w')
#obsids  = pd.read_pickle("snrs_obsids.pkl")

#sources = pd.read_pickle("snrs_sources.pkl")

sources = pd.read_excel("snrs.xlsx")


df_lines = pd.read_pickle("learned_lines.pkl")
#print(df_lines.columns.values)
#quit()


bkgs = pd.read_csv("backgrounds.csv")
#print(bkgs.head())

obsids_ex  = pd.read_excel("snrs_obsids.xlsx")
obsids = obsids_ex

#print(len(obsids_ex["OBSID"].unique()))

snrgreen = ascii.read("snrgreen.dat",delimiter="|")
regions = ascii.read("regions.txt",delimiter=",").to_pandas()

main_catalog_table = os.path.join(config['root_url'],"catalog_snr_table.html")


ftp_url  = config['ftp_url']
img_url  = config['img_url']
html_url = config['html_url']

#print obsids.head()

#sources.to_excel("snrs.xlsx")

#df = df_lines[df_lines["obsid"] == "100003010" ]

#print df

#quit()

master = ascii.read("master.dat").to_pandas()
htmlf = ascii.read("html_files.dat").to_pandas()
#jsf = ascii.read("js_files.dat").to_pandas()

anobs = df_lines['obsid'].unique()

all_obs = obsids_ex['OBSID'].astype(str).unique()

total_obs = len(all_obs)

print(("Number of spectra:",len(df_lines['spectrum'].unique())))
print(("Total Number of spectra:",total_obs))
print(("Number of observations:",len(anobs)))
irrel_obs = obsids_ex.loc[obsids_ex['RELEVANCE'] == 0,'OBSID'].astype(str).unique()
irrel = len(obsids_ex.loc[obsids_ex['RELEVANCE'] == 0,'OBSID'].astype(str).unique())
#print("Number of irrelevant observations:",irrel)
galcen_obs = obsids_ex.loc[obsids_ex['RELEVANCE'] == 2,'OBSID'].astype(str).unique()
gc_obs = len(obsids_ex.loc[obsids_ex['RELEVANCE'] == 2,'OBSID'].unique())
contam_obs = obsids_ex.loc[obsids_ex['RELEVANCE'] == 4,'OBSID'].astype(str).unique()
cont_obs = len(obsids_ex.loc[obsids_ex['RELEVANCE'] == 4,'OBSID'].unique())



crab_obs = obsids_ex.loc[obsids_ex["SNR"] == "SNR G184.6-05.8",'OBSID'].astype(str).unique()
num_crab_obs = len(crab_obs)

#print("Number of Crab observations:",num_crab_obs)


#print("Number of Galactic Center observations:",gc_obs)


nodata = obsids_ex.loc[obsids_ex['RELEVANCE'] == 3,'OBSID'].astype(str).unique()
nodata_obs = len(obsids_ex.loc[obsids_ex['RELEVANCE'] == 3,'OBSID'].unique())

#print("Number of obids without real data:",nodata_obs)

analyse_list_1 = [ob for ob in anobs if ob not in all_obs]

#print(analyse_list_1)


obs = []
status = []
name = []
snrs = []
for o in all_obs:

    obs.append(o)
    name.append(",".join(obsids_ex.loc[obsids_ex['OBSID'] == int(o),"Name"].astype(str).unique()))
    snrs.append(",".join(obsids_ex.loc[obsids_ex['OBSID'] == int(o),"SNR"].astype(str).unique()))

    if o in anobs:
        status.append("Analysed")
        continue
    if o in irrel_obs:
        status.append("Irrelevant")
        continue
    if o in crab_obs:
        status.append("Crab")
        continue
    if o in galcen_obs:
        status.append("Galactic center")
        continue
    if o in nodata:
        status.append("No data")
        continue
    if o in contam_obs:
        status.append("Contaminated/Bright Source in the FOV")
        continue
    status.append("Need to analyse")

s = Counter(status)
print(s)
pdf = pd.DataFrame({'obsid':obs,"status":status,'name':name,"snr":snrs})
pdf.to_excel("obsid_to_analyse.xlsx")

#print(len(analyse_list_2))
#quit()

anobsint = [int(o) for o in anobs]

#print df_lines.column
#grouped_lines = df_lines.loc[:,['obsid','linee']].groupby('obsid').count()
grouped_spec = df_lines.loc[:,['obsid','spectrum']].groupby('obsid').count()
#grouped_snr = df_lines.loc[:,['obsid','linee']].groupby('obsid').value_counts()

#print(grouped_lines.index)
#obsids_unique = obsids_ex["OBSID"].unique()


obsids_ex['OBSID']= obsids['OBSID'].astype(str)

#grouped2 = obsids_ex.loc[:,['OBSID','RELEVANCE']].groupby('OBSID').sum()
#print(grouped2.index)
#grouped2['obsid'] = grouped2['OBSID']
#joined1 = grouped2.join(grouped_spec,how='left')
#joined2 = joined1.join(grouped_spec,how='left')

#joined1.to_excel('obsids_joined.xlsx')

#print(obsids.describe())
#print("Unique Observations: {}".format(len(obsids["OBSID"].unique())))
#print("Unique Observations with Analysis: {}".format(len(obsids[obsids["OBSID"].isin(anobsint)]["OBSID"].unique())))
#print(obsids.head())

#print(sources.head())
#for s in sources.itertuples():

#    dfn = df_lines[df_lines["snr"] == s[1]]
#    if len(dfn["spectrum"].unique()) > len(dfn["obsid"].unique()):
#       print s[1],len(dfn["spectrum"].unique()),len(dfn["obsid"].unique())


#print(obsids.head())

#for s in obsids["OBSID"].unique():
#    print(s[1])
#    dfn = df_lines[df_lines["obsid"] == str(s)]
#    if len(dfn["snr"].unique()) > 1:
#       print s," ".join(dfn["snr"].unique()),len(dfn["snr"].unique())-1

#quit()

#casa = df_lines[df_lines["snr"] == "SNR G111.7-02.1"]
#tycho = df_lines[df_lines["snr"] == "SNR G120.1+01.4"]
#tycho = df_lines[df_lines["snr"] == "SNR G120.1+01.4"]

#n807008010 = df_lines[df_lines["obsid"] == "807008010"]
#n500006010 = df_lines[df_lines["obsid"] == "500006010"]
#n707020010 = df_lines[df_lines["obsid"] == "707020010"]



#print("Unique OBSIDs:",len(df_lines['obsid'].unique()))
#quit()

#print("Unique CASA Observations: {}".format(len(obsids["OBSID"].unique())))
#print("Unique CasA Spectra with Analysis: {}".format(len(casa["spectrum"].unique())))
#print("Unique Tycho Spectra with Analysis: {}".format(len(tycho["spectrum"].unique())))
#print("Unique Tycho Spectra with Analysis: {}".format(len(tycho["spectrum"].unique())))
#print("Unique 807008010 Spectra with Analysis: {}".format(len(n807008010["spectrum"].unique())))
#print("Unique 500006010 Spectra with Analysis: {}".format(len(n500006010["spectrum"].unique())))
#print("Unique 707020010 Spectra with Analysis: {}".format(len(n707020010["spectrum"].unique())))

#print(casa["obsid"].unique())

#quit()




#print df_lines.head()


reg_desc_dict = {
    "std":{0:"Entire XIS",1:"Circle 100\" centered on the XIS center",2:"XIS - central 100 \"circle excluded",3:"Circle 260\" centered on the xis center",2:"XIS - central 260 \"circle excluded"},
    "stdprod":{0:"Entire XIS",1:"Circle 100\" centered on the XIS center",2:"XIS - central 100 \"circle excluded",3:"Circle 260\" centered on the xis center",2:"XIS - central 260 \"circle excluded"},
    "snr":{0:"Background (Entire SIMBAD sources excluded)",1:"SIMBAD source 1",2:"SIMBAD source 2",3:"SIMBAD source 3"},
    "pie":{0:"central region",1:"North",2:"North East",3:"East",4:"South East",5:"North East",6:"South",7:"South West",8:"East"},
    "ann":{0:"central region",1:"Annulus 1",2:"Annulus 2"}
}

img_css = """
<style>
div.img {
    margin: 5px;
    border: 1px solid #ccc;
    left;
    width: 400px;
}

div.img:hover {
    border: 1px solid #777;
}

div.img img {
    width: 100%;
    height: auto;
}

div.desc {
    padding: 20px 20px 20px 20px;
    text-align: center;
}

div.imsp {
    padding: 5px;
    text-align: center;
    float:left;
}

</style>
"""


table_css = """
<style>
#snrtable {
    border-collapse:collapse;
    padding: 5px;
}

table {
    width:90%;
    border: 1px solid #ddd;
    border-collapse:collapse;
    padding: 5px;
}

th, tr, td {
    font-size:12px;
    border: 1px solid #ddd;
    padding: 5px;
}

tr.badobs {
color:red;
}

tr:hover {background-color: #f5f5f5}
tr:nth-child(even) {background-color: #f2f2f2}

td.heading {
    font-family:"Helvetica";

}
td.link:hover {background-color: #ff9999;}

.bottom {
    font-style: italic;
}

a.link {
    padding:5px 2px;
    display:block;
    color:blue;
    text-decoration:none;
}

a.tree_link {
    padding:5px 2px;
}

a.tree_link:hover { 
    text-color: red;
}

.tree_container {

    font-family: "Arial", Arial, sans-serif;
    font-size: 20px;
    
}


</style>
"""

#print obsids
#print sources

#print htmlf

fh = open(main_catalog_table,"w")
#fjs = open("snr_obs_list.js","w")

fh.write(table_css)
fh.write("<center><h2>SNRs observed by Suzaku</h2></center>\n")
fh.write("<table><tr>\n")
fh.write("<td><a href=\"guide/suzaku_snr_catalog_guide.html\">Documentation</a></td>")
#fh.write("<td><a href=\"suzakuproc4.html\">SNR Table 1</a></td>")
#fh.write("<td><a href=\"suzaku_snr_catalog.html\">SNR Table 2</a><br><br><hr></td>")
fh.write("</tr></table>\n")
fh.write("""<span class="tree_container">/CATALOG</span>""")

fh.write('<table border="1" id="snrtable">\n')
fh.write("""    
    <tr>
        <th width="14%">SNR</th>
        <th width="14%">ALT. NAME</th>
        <th width="10%">RA (h:m:s)</th>
        <th width="10%">DEC (d:m:s)</th>
        <th width="5%">SIZE, arcmin</th>
        <th width="10%">AGE, yr</th>
        <th width="10%">DISTANCE, kpc</th>
        <th width="5%">SN TYPE</th>
        <th width="8%">MORPHOLOGY</th>
        <th width="3%">N. OBS.</th>
        <th width="3%">Contam.</th>
    </tr>""")

#fjs.write('var point_src = "P";\n')
#fjs.write('var vext_src = "D";\n')
#fjs.write('var ext_src = "E";\n')
#fjs.write('\nvar sources = [\n')


#src_list = []
#obs_list = []
src_ind = 0

#def analysis_link(f):
img_dict = {"mos":"M","std":"S"}
spec_dict = {"ann":"A","snr":"R","stdprod":"S","std":"S","pie":"P","total":"T","snr1":"R","snr1amin":"R"}

#    print code,f
#    if f.find("img") > 0:
#    else:

#    return '<a href="{}">{}</a>'.format(f,letter)

for s in sources.itertuples():

    name = s[1]
#    print s[]

    greensource = snrgreen[snrgreen["NAME"] == name]
#    print greensource["CLASS"][0]

    if s[2].__class__.__name__ == "float64" and \
        s[3].__class__.__name__ == "float64":
        coord = SkyCoord(ra=s[2],dec=s[3],unit=(u.deg,u.deg),frame="fk5").to_string('hmsdms').split()
        ra = coord[0].replace('h',' ').replace('m','\' ').replace('s','"')
        dec = coord[1].replace('d',' ').replace('m','\' ').replace('s','"')
        print(ra,dec)

    else:
        ra = s[2]
        dec = s[3]
        coord = SkyCoord(ra=s[2],dec=s[3],unit=(u.deg,u.deg),frame="fk5").to_string('dms').split()
        ra = coord[0].replace('d',' ').replace('m','\' ').replace('s','"')
        dec = coord[1].replace('d',' ').replace('m','\' ').replace('s','"')
        print(ra,dec)


    try:
        rad = "{:2.2f}".format(s[7])
    except:
        rad = s[7]
    try:
        if isnan(rad): rad = ""
    except:
        pass

    if rad == "nan": rad = ""

    common_name = str(s[9])

    if common_name == "nan":
        try:
            common_name = s[5].replace("NAME ","").split(",")[0].replace("/","_").replace("(","").replace(")","")
        except:
            common_name = s[1]


    obs_sel = obsids[obsids["SOURCE"] == src_ind]

    obs_ind = 0
    num_obs = 0

#    print common_name
    sname = "_".join(common_name.lower().split())
#    if common_name == "":
#       sname = "_".join(name.lower().split())

    link = os.path.join(config['root_url'],"snr_html/"+sname+".html")
    fhobs = open(link,"w")
    fhobs.write(table_css)
    fhobs.write("<center><h2>Suzaku Observations of {0} ({1})</h2></center>\n".format(name,common_name))
#    fhobs.write("<center><a href=""><a href=\"../snr_ref_table.html#{0}ref\">References</a></center><br>".format(s[0]))
#    fhobs.write("")
    fhobs.write("""<span class="tree_container"><a href=\"{0}\" class="tree_link">/CATALOG</a>/{1}</span>""".format(main_catalog_table,common_name))
    fhobs.write('<table border="1" id="snrtable">')
    fhobs.write("""
    <tr>
        <th>OBSID</th>
        <th>SUZAKU TARGET</th>
        <th>RA POINT.</th>
        <th>DEC POINT.</th>
        <th>EXPOSURE</th>
        <th>XIS</th>
        <th>N. LINES</th>
    </tr>""")

    age = ""
    dist = ""
    sn_type = ""
    morph = ""

    df_snr = df_lines[(df_lines["snr"] == name)]

    try:
        morph = " ".join(greensource["CLASS"][0].split()[1:])
    except:
        pass
#       morph = snrgreen[snrgreen["NAME"] == name]["CLASS"].values[0]
    try:
        sn_type = lines.source_types[name]
    except Exception as e:
#       print e,e.message
        pass
    try:

        inst = df_snr["instrument"].values[0].replace("xi","")
        age_min = df_snr["age_min"].values[0]
        age_max = df_snr["age_max"].values[0]
        age = "{:.0f} - {:.0f}".format(age_min,age_max)
        if abs(age_min - age_max) < 0.0001: age = "{:.0f}".format(age_min)
#        print(abs(age_min - age_max))
        if isnan(age_min) and isnan(age_max): age = ""


        dist_min = df_snr["dist_min"].values[0]
        dist_max = df_snr["dist_max"].values[0]
        dist = "{} - {}".format(dist_min,dist_max)
        if abs(dist_min - dist_max) < 0.0001: dist = "{}".format(dist_min)
        if isnan(dist_min) and isnan(dist_max): dist = ""
#       sn_type = df_snr["type"].values[0]

    except:
        inst = ""
        lnum = ""



    for obs in obs_sel.itertuples():

        num_obs += 1
        bad = -1
        try:
            bad = int(obsids_ex.ix[obs[0],'RELEVANCE'])
        except:
            pass
#       print bad
#       print obs
        obsid = int(obs[1])
        mas = master[master["OBSID"] == obsid ]
#       obs_snrs = obsids[obsids["OBSID"] == obsid]["NAME"].unique()
#       obs_snrs_num = len(obs_snrs)


        if inst == "":

            try:
                inst = regions.loc[ (regions["obsid"] == obsid) & (regions["snr"] == name) ,"xis"].values[0]
            except:
                pass


        obs_snrs = sources.iloc[obsids[obsids["OBSID"] == obsid]["SOURCE"].unique()]['NAME'].values
        obs_snrs_num = len(obs_snrs)
#       print obs_snrs,obs_snrs_num


        szsrc = mas["NAME"].values[0]
        date_mjd = mas["TIME"].values[0]

#       ra_obs = mas["RA"].values[0]
#       dec_obs = mas["DEC"].values[0]

        coord_obs = SkyCoord(ra=mas["RA"].values[0],dec=mas["DEC"].values[0],unit=(u.deg,u.deg),frame="fk5").to_string('hmsdms').split()
        ra_obs = coord_obs[0].replace('h',' ').replace('m','\' ').replace('s','"')
        dec_obs = coord_obs[1].replace('d',' ').replace('m','\' ').replace('s','"')

        exp = mas["EXPOSURE"].values[0]

#       print obsid.__class__.__name__

        df_l = df_snr[ (df_snr["obsid"] == str(obsid))]

        regcode = ""
        regnum = 1

        lnum = df_l[df_l["linee"] > 0.35 ]["linee"].count()
        lnum_link = ""

        obs_regions = df_l["regcode"].unique()

        print(sname,obsid,name,inst,obs_regions)

        products_url = os.path.join(config['html_url'],"{0}/html/ae{0}_{1}.html".format(obsid,sname))


#  ******** EMMISSION MODEL (LINE LISTS) HTML FILES

        for r in obs_regions:

            df_reg = df_l[df_l["regcode"] == r]
            df_snr_name = df_reg["snr"].values[0]
            df_image_file = df_reg["imagefile"].values[0]

            if df_image_file == "none":
                fit_image_link = ""
                xcm_script = ""
            else:
                spec_list.write("mkdir {0}/spec\n".format(obsid))
                new_image_file = "spec/ae{0}_spec_{1}_img.gif".format(obsid,r)
                fit_image_link = """<a href="{2}/{1}/{0}"><img src="{2}/{1}/{0}" width=300 height=300></a>""".format(new_image_file,obsid,img_url)
                spec_list.write("cp {0}/fit/{1} {0}/{2}\n".format(obsid,df_image_file,new_image_file))
                xcm_script_file = "fit/{0}.xcm".format(df_image_file.split(".")[1])
                new_script_file = "ae{0}_spec_{1}.xcm".format(obsid,r)
                xcm_script = "<tr><td>XSPEC SCRIPT</td><td><a href=\"{2}/{1}/spec/{0}\">{0}</a></td></tr>".format(new_script_file,obsid,ftp_url)
                spec_list.write("cp {0}/{1} {0}/spec/{2}\n".format(obsid,xcm_script_file,new_script_file))

            df_for_name = sources[sources["NAME"] == df_snr_name]

            if len(df_for_name) == 0:
#               print "***{0}*****".format(df_snr_name)
                df_common_name = "NETY"
            else:
                df_common_name = str(df_for_name["commonname"].values[0])

                if df_common_name == "nan":
                    try:
                        df_common_name = df_for_name["ALT_NAMES"].values[0].replace("NAME ","").split(",")[0].replace("/","_").replace("(","").replace(")","")
                    except:
                        df_common_name = df_for_name["NAME"].values[0]


            lines_link_url = os.path.join(config['html_url'],"{0}/html/lines_{1}.html".format(obsid,r))
            fhl = open(lines_link_url,"w")
            fhl.write(table_css)

            contam_html = ""
            cnum = 0
            for s in obs_snrs:
                if s != df_snr_name:
                    contam_html +=  "<tr><td colspan=\"2\"><center>{}</center></td></tr>".format(s)
                    cnum += 1
            if cnum > 0: contam_html = "<tr><td colspan=\"2\">Contaminating sources</td></tr>".format(cnum) + contam_html


#           if obs_snrs_num > 1:


#               snr_arr = copy.copy(obs_snrs)
#               print snr_arr
#               while df_snr_name in snr_arr: snr_arr.remove(df_snr_name)
#               contam_html = "<td>{}</td>".format(" ".join(snr_arr))

#           fhl.write("<center><h2>Spectral lines observed for {0} ({1}) in </h2></center>\n".format(df_snr_name,df_common_name))
#           fhl.write("<center><h2>Observation sequence <a href=\"ae{0}_{2}.html\">{0}</a>, Region {1}</h2></center>\n".format(obsid,r,"_".join(df_common_name.lower().split())))
#           fhl.write("<br><a href=\"ae{0}_{1}.html\"><- Back to observations list for {0}</a>".format(obsid))

            fhl.write("""<colspan  class="tree_container"><a href=\"{0}\" class="tree_link">/CATALOG</a>/<a href=\"{1}\" class="tree_link">{4}</a>/<a href=\"{3}\" class="tree_link">{2}</a>/SPECTRAL FIT </span>""".format(main_catalog_table,link,obsid,products_url))
            fhl.write("<table>")
#            fhl.write("""<tr ><td colspan=2  class="tree_container"><a href=\"../../../{0}\" class="tree_link">CATALOG</a> * <a href=\"../../../{1}\" class="tree_link">{4}</a> * <a href=\"ae{2}_{3}.html\" class="tree_link">{2}</a> * SPECTRAL FIT </td></tr>""".format(main_catalog_table,link,obsid,sname,common_name))
            fhl.write("<tr ><td colspan=2  class=\"heading\"> <center> <h2>Spectral lines information</h2> </center></td></tr>")
            fhl.write("""<tr>
    <td><table>
            <tr><td>SNR</td><td class="link"><a href="{8}" class="link">{0}</a></td></tr>
            <tr><td>ALT. NAME</td><td>{1}</td></tr>
            <tr><td>Observation sequence</td><td class="link"><a href={9}" class="link">{2}</a></td></tr>
            <tr><td>Region</td><td>{3}</td></tr>
            {6}
            {7}
        </table>
    </td>
    <td>{5}</td>
    </tr>
""".format(df_snr_name,
           df_common_name,
           obsid,
           r,
           "_".join(df_common_name.lower().split()),
           fit_image_link,
           xcm_script,
           contam_html,link,products_url))

            fhl.write("</table><br>")

# formatters
            def fmt1(x): 
                if x == None: return ""
                return ("%.3e"%x).strip(".e+00")



            df_html = df_reg.loc[:,["id","linee","linee_min","linee_max","sigma","sigma_min","sigma_max","norm","norm_min","norm_max","ew","ew_min","ew_max"]].sort("linee")
#            df_html.rename({"id":"Line ID","linee":"Energy, keV"}, inplace=True)
            df_html.columns = ["Line ID","Energy, keV","Emin, keV","Emax, keV","Sigma, keV","Sigma min","Sigma max","norm","norm min","norm max","ew","EW min","EW max"]
            fhl.write(df_html.to_html(index=False,
                                  formatters=[
                                  lambda x:x,
                                  fmt1,fmt1,fmt1,
                                  fmt1,fmt1,fmt1,
                                  fmt1,fmt1,fmt1,
                                  fmt1,fmt1,fmt1]
            
            ))
            fhl.write("<center><h2>Continuum charactristics</h2></center>\n")
            fhl.write(df_reg.loc[:,["pl_model","pl_index","pl_indmin","pl_indmax","nh","nh_min","nh_max"]][0:1].
                      to_html(index=False,
                              formatters=[
                                  lambda x:x,
                                  fmt1,fmt1,fmt1,
                                  fmt1,fmt1,fmt1,
                                  fmt1,fmt1,fmt1,
                                  fmt1,fmt1,fmt1]
                     )
            )
#       fhl.write("""    <tr>
#        <th>Line ID</th>
#        <th>LineE</th>
#        <th>E Min.</th>
#        <th>E Max.</th>
#        <th>Sigma</th>
#        <th>Sigma Min.</th>
#        <th>Sigma Max.</th>
#        <th>Norm</th>
#        <th>Norm Min.</th>
#        <th>Norm Max.</th>
#    </tr>""")

#           fhl.write("<a href=\"../../../{}\"><- Back to the main SNR table</a><br>".format(main_catalog_table))
#            fhl.write("<br><a href=\"ae{0}_{1}.html\"><- {1}/{0} data products page</a>".format(obsid,sname))
#            fhl.write("<br><hr><span class=\"bottom\">Suzaku Supernova Remnants Catalog</span>".format(obsid))




#        if lnum>0:

#           lnum_link = "<a href=\"../obsids/{0}/html/lines_{2}.html\">{1}</a>".format(str(obsid),lnum,sname)

#           regcode = df_l["regcode"].values[0]
#           try:
#               regnum = int(regcode[1])
#           except Exception as e:
#               print e.message





        fhspec = open(products_url,"w")
        fhspec.write(img_css)
        fhspec.write(table_css)
        fhspec.write("<center><h2>Data products for {0} from observation {1}</h2></center>\n".format(common_name,str(obsid)))
        fhspec.write("""<span class="tree_container"><a href=\"{0}\" class="tree_link">/CATALOG</a>/<a href=\"{2}\" class="tree_link">{3}</a>/{1} </span>""".format(main_catalog_table,obsid,link,common_name))
        fhspec.write("""<table style="width:90%;">
        <tr><td>
        <table>
                <tr><td>Observation sequence</td><td></td><td>{0}</td></tr>
                <tr><td>RA </td><td>deg.</td><td>{1}</td></tr>
                <tr><td>DEC </td><td>deg.</td><td>{2}</td></tr>
                <tr><td>DATE </td><td>MJD</td><td>{3}</td></tr>
                <tr><td>EXPOSURE </td><td>sec.</td><td>{4}</td></tr>
                <tr><td>TARGET</td><td></td><td>{7}</td></tr>
        </table>
        </td>
        <td>
        """.format(obsid,ra_obs,dec_obs,date_mjd,exp,link,common_name,szsrc))

        html_files = htmlf[htmlf["file"].str.contains("ae"+str(obs[1]))]
        js_files = jsf[jsf["file"].str.contains("ae"+str(obs[1]))]
        js_files_snr = js_files[js_files["file"].str.contains("ae"+str(obs[1])+"_snr")]
        

# insert mosaic image if available

        for f in html_files[html_files["file"].str.contains("mos")]["file"].values:
            fhspec.write("""

        <div class="img">
                <a target="_blank" href="ae{0}_mos.png">
                        <img src="ae{0}_mos.png" alt="Mosaic" width="600" height="600">
                </a>
            <div class="desc">{0} in {1} Mosaic</div>
        </div>

""".format(str(obsid),common_name))
        comment = obsids_ex.ix[obs[0],'COMMENT']
#       if isnan(comment): comment = ""

        fhspec.write("</td></tr><tr><td>Comments</td><td><div id=\"comments\">{}</div></td></tr>".format(str(comment).replace("nan","")))
        analysis_html = "<table><tr>"
        img_html = ""
##      print html_files

        for f in js_files["file"].values:
    	    
            try:
##              link = analysis_link(f)
                code = f.split("_")[1]
                code = code.split("_")[0]
                code = code.split(".")[0]
                code = code.lower()
                if code == "stdprod" and len(js_files_snr) > 0:    continue

##              fjs.write('   {\n')

                if f.find("img")>0:
                    letter = img_dict[code]
                    img_html = img_html+'<td class="link"><a href=\"{3}/{0}/html/{1}\"  class="link">{2}</a> '.format(obsid,f.split("/")[-1],letter,config['img_url'])
#                   img_json = img_json+'{code:"'+letter+'",name:"'+code+'"},'
                    pass
                else:

                    letter = spec_dict[code]
                    analysis_html = analysis_html+'<a href=\"{3}/{0}/html/{1}.html\">{2}</a> '.format(obsid,code,letter,config['img_url'])

                    tarball = "ae{0}_{1}.tar.gz".format(obsid,code)

                    if code == "snr":
                        region_description = "Region  generated by SIMBAD database search"
                        region_number =  1

                    if code in ["std","stdprod"]:

                        region_description = "Standard Regions (circular regions centered on the XIS center)"
                        region_number = 5

                    if code == "pie":

                        region_description = "Pie Regions (pie pattern centered on SNR center)"
                        region_number = 9

                    if code in ["ann","annuli"]:

                        region_description = "Annuli Regions (set of annulus centered on SNR center)"
                        region_number = 9

                    if code == "man":

                        region_description = "Regions generated manually by visual inspection"
                        region_number = 1

                    snrs_html = """
                    <tr>
                        <td colspan=2 style="background-color:#99ccff;">{0}</td>
                    </tr>
                    """.format(region_description)



                    for snr_ind in range(region_number):

#                       print code.lower(),snr_ind,sname


                        isnr = snr_ind
                        rcode = "r"+str(isnr)+letter.lower()

                        if code == "snr":
                            try:
                                rcode = regions.loc[ (regions["obsid"] == obsid) & (regions["snr"] == name) ,"region"].values[0]
#                               print name,obsid,rcode
                            except Exception as e:
                                # if search in snr regions.txt fails look into regcode in main table
                                # which is stored in obs_regions list
                                try:
                                    reglist = [r for r in obs_regions if r[2] == "r"]
                                    rcode = reglist[0]
                                except:
                                    rcode = "r1r"
#                               print e.message,name,name.__class__.__name__,obsid.__class__.__name__,rcode
#                                pass
#                       if

                        image_links_html = "<table><tr><td>Images:</td>"
                        spec_links_html = "<table><tr><td>Spectra:</td>"
#                       print rcode

#                       fit_link = ""


                        for xi in inst.split():

                            if xi == "1" and (snr_ind == 1 or rcode in obs_regions):
                                pass
                            else:
                                image_links_html += "<td class=\"link\"><a href=\"{0}/{1}/html/ae{1}xi{2}_{3}_reg.png\" class=\"link\">XIS {2} <a></td>".format(ftp_url,str(obsid),xi,rcode)
                                spec_links_html += "<td class=\"link\"><a href=\"{0}/{1}/html/ae{1}xi{2}_{3}_src_nxb_spec.png\" class=\"link\">XIS {2}</a></td>".format(ftp_url,str(obsid),xi,rcode)


                        if rcode in obs_regions:
#                           print rcode
                            spec_links_html += "<td class=\"link\"><a href=\"{}\">SPECTRAL FIT</a></td>".format(lines_link_url)
#                           print spec_links_html
                        try:
                            reg_desc = reg_desc_dict[code][snr_ind]
#                           print reg_desc
                        except Exception as e:
#                           print e.message
                            reg_desc = ""

#                       if snr_ind == 1 or rcode in obs_regions:

                        sdesc = ""
                        if code == "snr":
                            sdesc  =" - {0}{1}".format(name+"/",common_name)
                        snrs_html += """

                        <tr>
                            <td colspan=2> Region {2} - {4}</td>
                        </tr><tr>
                            <td>
                                <div class="img">
                                    <a target="_blank" href="ae{3}xi0_{2}_reg.png">
                                        <img src="ae{3}xi0_{2}_reg.png" alt="XIS 0 Image" width="400" height="300">
                                    </a>
                                    {0}</tr></table>
                                </div>
                            </td>
                            <td>
                                <div class="img">
                                    <a target="_blank" href="ae{3}xi0_{2}_src_nxb_spec.png">
                                        <img src="ae{3}xi0_{2}_src_nxb_spec.png" alt="XIS 0 Image" width="400" height="300">
                                    </a>

                                    {1}</tr></table>
                                </div>

                            </td>
                        </tr>
""".format(image_links_html,spec_links_html,rcode,str(obsid),sdesc)

#                       else:
#                           snrs_html += """
#                       <tr>
#                           <td colspan=2> Region r{4}{2} - {5}</td>
#                       </tr><tr>
#                           <td>
#                                   {0} </tr></table>
#                           </td>
#                           <td>
#                                   {1} </tr></table>
#                           </td>
#                       </tr>
#""".format(image_links_html,spec_links_html,letter.lower(),obsid,isnr,reg_desc)


                    fhspec.write("""
{3}
    <tr>
        <td>Data Products</td>
        <td><a href="{1}/obsids/{2}/{0}">{0}</a></td>
    </tr>
""".format(tarball,ftp_url,obsid,snrs_html))
#                   fhspec.write("<a href=\"../../../{0}\"><- Back to observations list for {1}</a>".format(link,common_name))
#                   fhspec.close()

#                   analysis_json = analysis_json+'{code:"'+letter+'",name:"'+code+'"},'
#               fjs.write('   },\n')
            except Exception as e:
                print(e)
                pass

        bkgrec = bkgs[bkgs["oid"] == obsid]

        try:

            fhspec.write("""
                <tr>
                    <td>Background</td>
                    <td colspan=2><a href="{1}/obsids/{0}/ae{0}_bkg.tar.gz">ae{0}_bkg.tar.gz</a></td>
                </tr>""".format(bkgrec.iloc[0,2],ftp_url))
        except:
            pass

        fhspec.write("""</table>\n""")
#       fhspec.write("<script src=\"comments.js\"></script>")

        fhspec.close()


        if not os.path.isdir('obsids/'+str(obsid)+"/html"): os.mkdir('obsids/'+str(obsid)+"/html")
        badmark = ""
        if bad in [0,2,3]: badmark = 'class="badobs"'
        obs_html = """

   <tr {8}>
       <td class="link"><a href=\"../{7}\"  class="link">{0}</a></td>
       <td>{1}</td>
       <td>{2}</td>
       <td>{3}</td>
       <td>{4}</td>
       <td>{5}</td>
       <td>{6}</td>
    </tr>
""".format(obsid,szsrc,ra_obs,dec_obs,exp,inst,len(df_l),products_url,badmark)

#       if lnum > 0:

#           fhl = open('obsids/'+str(obsid)+"/html/lines_"+sname+".html","w")
#           fhl.write(table_css)
#           fhl.write("<center><h2>Spectral lines observed in {0}.</h2></center>\n".format(common_name,))
#           fhl.write("<center><h2>Observation sequence  <a href=\"ae{0}.html\">{0}</a>, Region {1}</h2></center>\n".format(obsid,regcode))
#           fhl.write("<br><a href=\"../../../{0}\"><- Back to observations list for {1}</a>".format(link,common_name))
#           fhl.write(df_l.loc[:,["id","linee","linee_min","linee_max","sigma","sigma_min","sigma_max","norm","norm_min","norm_max","ew","ew_min","ew_max"]].sort("linee").to_html(index=False))
#           fhl.write("<center><h2>Continuum charactristics</h2></center>\n".format(common_name,obsid))
#           fhl.write(df_l.loc[:,["pl_model","pl_index","pl_indmin","pl_indmax","nh","nh_min","nh_max"]][0:1].to_html(index=False))
#       fhl.write("""    <tr>
#        <th>Line ID</th>
#        <th>LineE</th>
#        <th>E Min.</th>
#        <th>E Max.</th>
#        <th>Sigma</th>
#        <th>Sigma Min.</th>
#        <th>Sigma Max.</th>
#        <th>Norm</th>
#        <th>Norm Min.</th>
#        <th>Norm Max.</th>
#    </tr>""")

#           fhl.write("<a href=\"../../../{}\"><- Back to the main SNR table</a><br>".format(main_catalog_table))#
#           fhl.write("<br><a href=\"../../../{0}\"><- Back to observations list for {1}</a>".format(link,common_name))
#           fhl.write("<a href=\"../../../{}\"><- Back to the main SNR table</a><br>".format(main_catalog_table))

#           fhl.close()
#       fjs.write('        {\n')
#       fjs.write('            obsid:\"{}\",\n'.format(str(obs[1])))
#       fit_code = ""
#       if str(obs[1]) in obsid_list: fit_code = "Y"
#       fjs.write('            fit:\"{}\",\n'.format(fit_code))
#
##      print obs[1]
#
#
#       stype = "point_src"
#       try:
#           if s[7]>1.0: stype = "ext_src"
#           if s[7]>20.0: stype = "vext_src"
#       except:
#           pass
#       fjs.write('            stype:{},\n'.format(stype))
#       fjs.write('            object:\"{}\",\n'.format(szsrc))
#
#       analysis_json =  analysis_json+"],\n"
#       img_json =  img_json+"],\n"
#       fjs.write(analysis_json)
#       fjs.write(img_json)

#       print html_files
#       obs_list.append([obs[1],obs[5],src_ind])
#       szsrc = master[master["OBSID"] == obs[1]]["NAME"].values[0]
        fhobs.write(obs_html)

        obs_ind = obs_ind + 1
#       print obs


    fhobs.write("</table><br><br>")

    fhobs.close()
    contam = ""
    if obs_snrs_num > 1: contam = str(obs_snrs_num-1)

    src_html = """
   <tr>
       <td><a href=\"{6}\">{0}</td>
       <td>{1}</td>
       <td style="white-space: nowrap;">{2}</td>
       <td>{3}</td>
       <td>{4}</td>
       <td>{7}</td>
       <td>{8}</td>
       <td>{10}</td>
       <td>{11}</td>
       <td>{5}</td>
       <td>{9}</td>
    </tr>
""".format(name,common_name,ra,dec,rad,num_obs,link,age,dist,contam,sn_type,morph)


        #fjs.write('        },\n')


    if obs_ind > 0:
#       src_list.append([s[2],s[3],s[4]," ",s[8],s[9],s[5]])
        src_ind = src_ind + 1
#       src_html = "<tr><td rowspan=\"{}\">{}</td>\n".format(str(obs_ind),s[1])+src_html[4:]



    fh.write(src_html)

#    fjs.write('    ]\n')
#    fjs.write('},\n')
fh.write("</table>")
#fjs.write("]")
fh.close()
spec_list.close()
#fjs.close()
