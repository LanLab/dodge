from time import sleep as sl
import argparse
from sklearn.cluster import AgglomerativeClustering
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import time
import itertools
from datetime import datetime,timedelta
from csv import reader
import sys
import os
from collections import OrderedDict
from collections import Counter
import glob
import calendar

from multiprocessing import Pool


class isolateOb:
    def __init__(self,name,date,year,month,new,strainmgtdict,investigation=False):
        self.name = name
        try:
            self.date = datetime.strptime(date, '%Y-%m-%d')
        except:
            self.date = None
        if year == '':
            self.year = None
        else:
            self.year = int(year) ## format '%Y'
        if month == '':
            self.month = None
        else:
            self.month = month ## format '%m'
        if month != "" and year != "":
            self.monthdate = datetime.strptime('{}-{}'.format(year, month), '%Y-%m')
        if new == 'False' or new == False:
            self.new = False
        elif new == 'True' or new == True:
            self.new = True
        self.partof = {}
        self.mgtid_dict = strainmgtdict
        self.investigation = investigation

    def __str__(self):
        return "Name:{}, Date:{}, Month: {}, Year: {}, new?: {}, investigation: {}, partof: {}, mgtids: {}".format(self.name,self.date,self.month,self.year,self.new,self.investigation,str(self.partof),str(self.mgtid_dict))

class clusterOb:
    def __init__(self, level,id,mgtid="",investigation=False):
        self.id = id
        if mgtid == "":
            self.mgtid = ""
        else:
            self.mgtid = mgtid
        self.level = level
        self.inheritedid = ""
        self.strains = []
        self.size = len(self.strains)
        self.contains = {}
        self.partof = {}
        self.status = ""
        self.investigation = investigation
        self.mgtlist = []
        self.distmatrix = ""
        self.pairwisemaxdistance = ""
        self.prevclusters = []

    def __str__(self):
        if self.size == 0:
            return "empty cluster object"
        else:
            return "ID:{}, Level: {}, Size: {}, Days: {}, Months: {}, Years: {}, Strains: {}, Investigation: {}".format(self.id,self.level,self.size,self.dayspan,self.monthspan,self.yearspan,",".join(self.strainnames),str(self.investigation))

    def addstrain(self, strain,strainmgtdict,distmatrix,strain_to_st):
        self.strains.append(strain)
        self.size = len(self.strains)
        self.dayspan = ""
        self.monthspan = ""
        self.yearspan = ""
        self.mindate = ""
        self.maxdate = ""
        self.minmonth = ""
        self.maxmonth = ""
        self.mgtlist.append(strainmgtdict)
        self.strainnames = [x.name for x in self.strains]

        # calculate timespan of cluster
        dates = [x.date for x in self.strains if x.date != None]
        if len(dates) > 1:
            self.mindate = min(dates)
            self.maxdate = max(dates)
            timespan = self.maxdate - self.mindate
            self.dayspan = timespan.days

        months = [str(x.year) + "-" + str(x.month)+"-01" for x in self.strains if x.month not in [None,""] and x.year not in [None,""]]
        if len(months) > 1:
            dateobjs = [datetime.strptime(x, '%Y-%m-%d') for x in months]
            minmonth = min(dateobjs)
            self.minmonth = minmonth
            maxmonth = max(dateobjs)
            self.maxmonth = maxmonth
            num_months = (maxmonth.year - minmonth.year) * 12 + (maxmonth.month - minmonth.month+1)
            self.monthspan = int(num_months)
        else:
            self.monthspan = 1

        # calculate yearspan of cluster
        years = [int(x.year) for x in self.strains if x.year != None]
        if len(years) > 1:
            earliest = min(years)
            latest = max(years)
            yearspan = latest - earliest+1
            self.yearspan = yearspan
        else:
            self.yearspan = 1

        # get pairwise distance matrix and maximum pairwise distance
        if isinstance(distmatrix, pd.DataFrame):
            # distance matrix is made up of MGT9 STs not strains need to convert
            strain_namelstmp = [x.name for x in self.strains]
            strain_namels = []
            [strain_namels.append(x) for x in strain_namelstmp if x not in strain_namels]
            if strain_to_st != {}:

                stlist = list(strain_namels)
            else:
                stlist = list(strain_namels)
            stlistint = list(map(str,stlist))
            try:
                self.distmatrix = distmatrix.loc[stlistint,stlist]
            except Exception as e:
                print(e)
                print("continuing because common type error, if not type error then could be real issue ")
                stlistint = list(map(int, stlist))
                self.distmatrix = distmatrix.loc[stlistint, stlist]
            self.distmatrix['newindex'] = self.strains
            self.distmatrix.set_index('newindex',inplace=True)
            self.distmatrix.columns = self.strains
            self.pairwisemaxdistance = self.distmatrix.values.max()

    def make_mgtid(self,nextint):

        # calculate cluster ID from largest MGT level with all the same ST
        mgtcounts = OrderedDict()
        for strain in self.strains:
            for mgtlevel in strain.mgtid_dict:
                st = strain.mgtid_dict[mgtlevel]
                if mgtlevel not in mgtcounts:
                    mgtcounts[mgtlevel] = [st]
                else:
                    mgtcounts[mgtlevel].append(st)
        if mgtcounts:
            newid = str(self.mgtid)
            for level in mgtcounts:
                lev_sts = mgtcounts[level]
                num_sts = len(list(set(lev_sts)))
                if num_sts == 1:
                    newid = "{} ST{}".format(level,mgtcounts[level][0])
                elif num_sts > 1:
                    topst,topst_freq = most_frequent(lev_sts)
                    if topst_freq > 0.70:
                        newid = "{} ST{}".format(level, topst)

            self.mgtid = str(newid)
        else:
            self.mgtid = str(nextint)
        return self.mgtid

    def make_hierccid(self,nextint):
        # calculate cluster ID from largest MGT level with all the same ST
        hcccounts = OrderedDict()
        for strain in self.strains:
            for hcclevel in strain.mgtid_dict:
                st = strain.mgtid_dict[hcclevel]
                if hcclevel not in hcccounts:
                    hcccounts[hcclevel] = [st]
                else:
                    hcccounts[hcclevel].append(st)
        hiercc_levels = list(hcccounts.keys())
        hiercc_levels = sorted(hiercc_levels,key=lambda x:int(x.replace("HC","")),reverse=True)
        if hcccounts:
            newid = str(self.mgtid)
            for level in hiercc_levels:
                lev_cs = hcccounts[level]
                num_cs = len(list(set(lev_cs)))
                if num_cs == 1:
                    newid = "{} {}".format(level,hcccounts[level][0])
                elif num_cs> 1:
                    topst,topst_freq = most_frequent(lev_cs)
                    if topst_freq > 0.70:
                        newid = "{} {}".format(level, topst)

            self.mgtid = str(newid)
        else:
            self.mgtid = str(nextint)
        return self.mgtid



    def write(self, path,args):
        if not os.path.isfile(path):
            outf = open(path, "w")
            outf.write("ID\tmgtid\tLevel\tSize\tMax distance\tTimespan\tMindate\tMaxdate\tStrains\tstatus\tInvestigation\tcontains\tpartof\n")
        else:
            outf = open(path, "a+")

        if args.timesegment == "week":
            dateobjmin = self.mindate
            dateobjmax = self.maxdate
            dateformat = '%Y-%m-%d'
            timespan = self.dayspan
        elif args.timesegment == "month":
            dateobjmin = self.minmonth
            dateobjmax = self.maxmonth
            dateformat = '%Y-%m'
            timespan = self.monthspan

        if isinstance(dateobjmin,datetime):
            mind = dateobjmin.strftime(dateformat)
        else:
            mind = dateobjmin

        if isinstance(dateobjmax,datetime):
            maxd = dateobjmax.strftime(dateformat)
        else:
            maxd = dateobjmax

        strainids = [x.name for x in self.strains]
        partof = [str(self.partof[level])+":"+str(level) for level in self.partof]
        contains = [str(self.contains[level]) + ":" + str(level) for level in self.contains]
        outf.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
            self.id,
            self.mgtid,
            self.level,
            self.size,
            self.pairwisemaxdistance,
            timespan,
            # ",".join(map(str, dates)),
            str(mind),
            str(maxd),
            ",".join(strainids),
            self.status,
            str(self.investigation),
            contains,
            partof)
        )
        outf.close()

    def writesummary(self, path,args):
        strainids = [x.name for x in self.strains]
        if not os.path.isfile(path):
            outf = open(path, "w")
            outf.write("ID\tMGTID\tLevel\tSize\tMax distance\tTimespan\tMindate\tMaxdate\tStrains\tstatus\tInvestigation\n")
        else:
            outf = open(path, "a+")

        if args.timesegment == "week":
            dateobjmin = self.mindate
            dateobjmax = self.maxdate
            dateformat = '%Y-%m-%d'
            timespan = self.dayspan
        elif args.timesegment == "month":
            dateobjmin = self.minmonth
            dateobjmax = self.maxmonth
            dateformat = '%Y-%m'
            timespan = self.monthspan

        if isinstance(dateobjmin,datetime):
            mind = dateobjmin.strftime(dateformat)
        else:
            mind = dateobjmin

        if isinstance(dateobjmax,datetime):
            maxd = dateobjmax.strftime(dateformat)
        else:
            maxd = dateobjmax
        # if self.mgtid == "MGT9 ST50288":
        #     sl(1)
        outf.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(self.id,self.mgtid, self.level, self.size,
                                                                 self.pairwisemaxdistance, timespan,
                                                                 mind,
                                                                 maxd,
                                                                 ",".join(strainids),
                                                                 self.status,str(self.investigation)))#, ",".join(olds)))
        outf.close()


def most_frequent(l):
    counter = Counter(l)
    most_common = counter.most_common(1)
    commonid = most_common[0][0]
    commoncount = most_common[0][1]
    freq = float(commoncount)/len(l)
    return commonid,freq

def unneg(a):
    if "-" in a:
        return int(a.split("_")[0][1:])
    else:
        return int(a)

def get_strain_to_st(args,strain,st):
    st_to_strain = {}
    strain_to_st = {}
    inprofiles = open(args.variant_data, "r").read().splitlines()
    for line in inprofiles[1:]:
        col = line.split("\t")
        if col[st] not in st_to_strain:
            st_to_strain[col[st]] = [str(col[strain])]
        else:
            st_to_strain[col[st]].append(str(col[strain]))
        strain_to_st[col[strain]] = col[st]
    return st_to_strain,strain_to_st

def import_allele_data(args,strains,backgroundstrains,hasdistance):
    """

    :param args:
    :param strains:
    :param backgroundstrains:
    :return:
    idlist = list of strains with data successfully gathered
    profs = {strain:allele profile for strain as list of numbers as strings}
    st_to_strain = {ST:[list of strains assigned ST]}
    strain_to_st = {strain:ST}
    """

    if args.enterobase_data:
        strain = 0
        st = 1
        profstart = 2
    else:
        strain = 0
        st = 1
        profstart = 3

    strainstmp = strains + backgroundstrains
    strains = []
    [strains.append(x) for x in strainstmp if x not in strains]

    if list(set(strains).difference(set(hasdistance))) == [] and len(hasdistance) > 0:
        print("\nAll isolates found in input pairwise distances file")
        st_to_strain,strain_to_st = get_strain_to_st(args, strain, st)
        return strains, {}, st_to_strain, strain_to_st, []


    inprofiles = open(args.variant_data, "r").read().splitlines()

    profs = {}
    st_to_strain = {}
    strain_to_st = {}

    ##per strain
    for line in inprofiles[1:]:
        col = line.split("\t")
        if col[strain] in strains:
            if col[strain] not in profs and col[st] != "":
                if args.enterobase_data:
                    noneg = ["0" if x in ['0','',"-"] else x for x in col[profstart:]]
                    noneg = [unneg(x) for x in noneg]
                else:
                    noneg = ["0" if x in ['0', '', "-"] else x for x in col[profstart:]]
                    noneg = [unneg(x) for x in noneg]
                profs[col[strain]] = noneg
                st_to_strain[col[st]] = [str(col[strain])]
            else:
                if col[strain] in profs and col[st] not in st_to_strain:
                    print(f"possible duplication of {col[strain]} in allele profiles")
                    continue
                st_to_strain[col[st]].append(str(col[strain]))

            strain_to_st[col[strain]] = col[st]

    idlist = []
    [idlist.append(x) for x in profs.keys() if x not in idlist]

    missing_profiles = [x for x in strains if x not in idlist]
    print("{} strains missing allele profiles: {}".format(len(missing_profiles),",".join(missing_profiles)))
    return idlist, profs, st_to_strain, strain_to_st,missing_profiles

def parse_mask_bed(args):
    bed = open(args.mask).read().splitlines()
    maskls = []
    for i in bed[1:]:
        col = i.split("\t")
        st = int(col[1])
        en = int(col[2])
        maskls += list(range(st,en+1))
    return maskls

def collect_SNPs(args,strains):
    """
    extract snps from snippy input files that correspond to list of strains in strains input
    :param args:
    :param strains:
    :return:
    snps = snps[position][(reference nucleotide,alternate nucleotide)] = [list of strains with SNP]
    refpos = dictionary of reference nucleotide for each position where snp is called
    namels = list of strains where input data was found
    """
    inp = args.variant_data + "/*.subs.vcf"
    filels = []
    if "*" in inp:
        filels = glob.glob(inp)

    if len(filels) == 0:
        sys.exit("input path {} found no vcf files".format(inp))

    maskls = []
    if args.mask:
        maskls = parse_mask_bed(args)

    snps = {}
    refpos = {}
    namels = []
    for file in filels:
        name = os.path.basename(file).replace(".subs.vcf","")
        if name in strains:
            namels.append(name)
            f = open(file,"r").read().splitlines()
            for line in f:
                if line[0] != "#" and "TYPE=snp" in line and "OLDVAR=" not in line:
                    col = line.split("\t")
                    pos = col[1]
                    ref = col[3]
                    mut = col[4]
                    qual = col[5]
                    if float(qual) < args.snpqual:
                        mut = "N"
                    if pos not in maskls:
                        if pos not in snps:
                            snps[pos] = {(ref,mut):[name]}
                            refpos[pos] = ref
                        else:
                            if (ref,mut) in snps[pos]:
                                snps[pos][(ref,mut)].append(name)
                            else:
                                snps[pos][(ref,mut)] = [name]
    return snps,refpos,namels

def get_genomes(strainls,args):
    """
    Get genome for each strain, .consensus.subs.fa files have "-" and "N" at correct positions relative to ref - allows
    checking for missing data
    :param strainls:
    :param args:
    :return: genomedict = {strain:string of subbed reference genome}
    """
    inp = args.variant_data + "/*.consensus.subs.fa"
    filels = []
    if "*" in inp:
        filels = glob.glob(inp)

    if len(filels) == 0:
        sys.exit("input path {} found no *.consensus.subs.fa files".format(inp))

    genomedict = {}
    for i in filels:
        strain = i.split("/")[-1].replace(".consensus.subs.fa","")
        if strain in strainls:
            ingenome = open(i,"r").read().splitlines()
            genomeseq = "".join(ingenome[1:])
            genomedict[strain] = genomeseq
    return genomedict

def make_alignment(snpdict,genomedict,strains,args,refpos):
    """
    Make dict of alignment  outalign[strain] = "stringOfAlignment"
    :param snpdict:
    :param genomedict:
    :param strains:
    :param args:
    :param refpos:
    :return:
    """
    outalign = {}
    if args.useref:
        outalign["Reference"] = ""
    for strain in strains:
        outalign[strain] = ""
    posls = sorted(map(int,list(snpdict.keys())))
    c=0
    for pos in posls:
        for strain in strains:
            toadd = ""
            for mut in snpdict[str(pos)]:
                if strain in snpdict[str(pos)][mut]:
                    toadd = mut[1]
            if toadd == "":
                if strain in genomedict:
                    call = genomedict[strain][pos-1]
                    toadd = call
            outalign[strain]+=toadd
        if args.useref:
            outalign["Reference"] += refpos[str(pos)]
        if c%100 == 0:
            print("{}/{} snp positions processed".format(c,len(posls)), end="\r", flush=True)
        c+=1
    print("{}/{} snp positions processed".format(c, len(posls)))
    return outalign,posls

def get_identical_alignments(align):
    outdict = {}
    id_to_strain = {}
    strain_to_id = {}
    idls = []
    for strain in align:
        alignment = align[strain]
        outdict[strain] = alignment
        idls.append(strain)
        id_to_strain[strain] = [strain]
        strain_to_id[strain] = strain


    return idls,outdict,id_to_strain,strain_to_id

def check_files_present(args,allstrains):
    inp = args.variant_data + "/*.subs.vcf"
    filels = []
    if "*" in inp:
        filels = glob.glob(inp)

    if len(filels) == 0:
        sys.exit("input path {} found no vcf files".format(inp))

    namels = []
    for file in filels:
        name = os.path.basename(file).replace(".subs.vcf", "")
        if name in allstrains:
            namels.append(name)
    return namels

def import_snp_data(args,strains,background_strains,hasdistance):
    allstrains = list(strains+background_strains)

    print("import_snpdata: import snp start")

    missing_distance = list(set(allstrains).difference(set(hasdistance)))
    strainswdata = check_files_present(args, allstrains)
    missing_distance_with_data = list(set(missing_distance).intersection(set(strainswdata)))
    if len(missing_distance_with_data) == 0:
        print("\nAll isolates with data are found in input pairwise distances file")
        return allstrains, {}, []
    else:
        print(f"to run: {missing_distance_with_data}")
        strainswdata = check_files_present(args, allstrains)
        missing_distance = list(set(strainswdata).difference(set(hasdistance)))
        if list(missing_distance) == []:
            print("\nAll isolates found in input pairwise distances file")
            return strainswdata, {}, []

        snpdict,refpos,snpfile_present = collect_SNPs(args,strainswdata)
        print("import_snpdata: snps collected")

        genomedict = get_genomes(strainswdata,args)

        usablestrains = []
        c=[]
        for strain in strainswdata:
            if strain in snpfile_present and strain in genomedict.keys() and strain not in usablestrains:
                usablestrains.append(strain)
            elif strain in usablestrains:
                continue
            else:
                c.append(strain)
        print("import_snpdata: {} strains without genome or vcf file ({})".format(len(c),",".join(c)))
        print("import_snpdata: genomes with subs collected")
        outalign, posls = make_alignment(snpdict, genomedict, usablestrains, args, refpos)
        print("import_snpdata: alignment generated")

        exist = []
        for i in usablestrains:
            if i not in exist:
                exist.append(i)
            else:
                print("{} duplicated in usablestrains".format(i))

        return usablestrains,outalign,c

def import_clusters(args,strainobjdict):
    """

    :param args: input args from argparse
    :param strainobjdict: dictionary of {strain ID: strain Object}
    :return:
    clusters: dictionary clusters[cluster level(cutoff)][cluster identity] = cluster object containing sets of isolate objects and other functions
    prevstrains: list of strains previously assigned to outbreak clusters
    """
    clusterids = []
    inclusters = open(args.inclusters).read().splitlines()
    clusters = {}
    prevstrains = []
    for line in inclusters[1:]:
        col = line.replace('"','').split("\t")

        mgtid = col[1]
        lev = int(col[2])
        size = int(col[3])
        ident = "old_"+str(col[0])+"_"+str(lev)

        strains = col[8].split(',')
        prevstrains += strains
        if col[10] == 'False':
            invest = False
        elif col[10] == 'True':
            invest = True
        else:
            invest = False

        if lev not in clusters:
            clusters[lev] = {}

        init = True
        for pos in range(0,size):
            strainobj = strainobjdict[strains[pos]]
            if init:
                clusters[lev][ident] = clusterOb(lev, ident, mgtid=mgtid, investigation=invest)
                clusters[lev][ident].addstrain(strainobj,{},"",{})
                init=False
            else:
                clusters[lev][ident].addstrain(strainobj,{},"",{})

    return clusters,prevstrains

def ident_expanded_old(oldclusters,newclusters):

    """
    at each level
    for each new cluster
    get list of old clusters within new
        to do this:
        get list of old cluster isolates + list of new cluster isolates
        if old = investigation:
            if new>old and intersect = old - >old expanded
            OR
            if new = old and intersect = old old unchanged

            record old in new

            if at end no old in new -> new cluster

    """

    new_clust_in_old_invest = {}
    unchanged_oldclust = {}
    if oldclusters == {}:
        return new_clust_in_old_invest,unchanged_oldclust
    # Check that clustering levels used in previous clustering run are the same as current run
    missingnew = [x for x in oldclusters.keys() if x not in newclusters.keys()]
    missingold = [x for x in newclusters.keys() if x not in oldclusters.keys()]
    if len(missingnew) > 0:
        sys.exit("Cluster levels {} between previous runs/background and this run are inconsistent\n"
                 "Make sure you specify matching -l values for background and all subsequent"
                 " identification runs".format(",".join(missingnew)))
    elif len(missingold) > 0:
        sys.exit("Cluster levels {} between previous runs/background and this run are inconsistent\n"
                 "Make sure you specify matching -l values for background and all subsequent"
                 " identification runs".format(",".join(missingold)))

    # for each cutoff level
    for level in oldclusters:
        # get old and new clusters
        levoldclusters = oldclusters[level]
        levnewclusters = newclusters[level]
        # get pairs of old and new clusters at the same level
        oldnew_pairs = itertools.product(levoldclusters.keys(),levnewclusters.keys())
        # for each pair of clusters
        for pair in oldnew_pairs:
            old = oldclusters[level][pair[0]]
            new = newclusters[level][pair[1]]
            if old.investigation: # if old cluster was assigned as an investigation cluster
                if old.strains == new.strains: # if the set of strains in old and new are the same
                    new.investigation = True
                    new.mgtid = old.mgtid
                    new.status = "unchanged"
                    if level not in unchanged_oldclust:
                        unchanged_oldclust[level] = [new]
                    else:
                        unchanged_oldclust[level].append(new)

                    # set new cluster to investigation, new cluster adopts old cluster id, status set to unchanged and add new cluster to unchanged list
                elif set(new.strains).intersection(set(old.strains)) == set(old.strains):# if the set of strains in old is a subset of new
                    new.investigation = True
                    new.mgtid = old.mgtid
                    new.status = "expanded"
                    # new_clust_in_old_invest.append(new)
                    if level not in new_clust_in_old_invest:
                        new_clust_in_old_invest[level] = [new]
                    else:
                        new_clust_in_old_invest[level].append(new)
                    # set new cluster to investigation, new cluster adopts old cluster id, status set to expanded and add new cluster to new_clust_in_old_invest list
    return new_clust_in_old_invest,unchanged_oldclust

def allele_dist_metric(a,b,args):
    match = 0
    missmatch = 0
    for i in range(len(a)):
        aAllele = a[i]
        bAllele = b[i]
        if aAllele == 0 or bAllele == 0 or aAllele == bAllele:
            match += 1
        else:
            missmatch += 1
            if missmatch >= args.max_missmatch:
                return missmatch
    return missmatch

def snp_dist_metric(a,b,args):
    match = 0
    missmatch = 0
    for pos in range(len(a)):
        anuc = a[pos]
        bnuc = b[pos]
        missing = ["N","n","X","x","-"]
        if anuc in missing or bnuc in missing or anuc == bnuc:
            match +=1
        else:
            missmatch +=1
            if missmatch >= args.max_missmatch:
                return missmatch
    return missmatch

def run_multiprocessing(func, i, n_processors):
    with Pool(processes=n_processors) as pool:
        return pool.starmap(func, i)

def pairdist(args,pair,existids,initiald,profs):
    id1 = pair[0]
    id2 = pair[1]
    if id1 in existids and id2 in existids:
        dist = int(initiald.at[str(id1), str(id2)])


    else:
        ap1 = profs[id1]
        ap2 = profs[id2]

        if args.inputtype == "snp":
            dist = snp_dist_metric(ap1, ap2, args)
        elif args.inputtype == "allele":
            dist = allele_dist_metric(ap1, ap2, args)

    return dist

def run_dist(args,profs,idlist,hasdistance,initiald,odc10_to_strains):
    # get existing distances if present, also get list of strain ids
    # make empty dataframe with all current strains

    initiald.index = initiald.index.astype(str)
    initiald.columns = initiald.columns.astype(str)

    if list(set(idlist).difference(set(hasdistance))) == []:
        print("\nAll isolates found in input pairwise distances file1")
        ## subset/reorder input df to match isolate list (idlist)
        newdf = initiald[idlist]
        newdf = newdf.T[idlist].T
        return newdf

    useodc10 = False

    exist = []
    for i in idlist:
        if i not in exist:
            exist.append(i)
        else:
            print("{} duplicated in idlist".format(i))


    ######################################################

    start_time = time.time()
    newdf = pd.DataFrame(index=idlist,columns=idlist,dtype=str)
    newdf.index = newdf.index.astype(str)
    newdf.columns = newdf.columns.astype(str)
    newdf.fillna(args.max_missmatch,inplace=True)

    exist = []
    for i in newdf.columns.tolist():
        if i not in exist:
            exist.append(i)
        else:
            print("{} duplicated in new df".format(i))


    pairs = itertools.combinations(idlist,r=2)
    pairs = list(pairs)
    tot = (len(idlist)*len(idlist)-len(idlist))/2
    frac = int(tot*0.01)
    if frac == 0:
        frac = 1
    c=0

    print("\nCalculating pairwise distances\nDone\t\t% done\t\ttime")

    inps = [(args,x,hasdistance,initiald,profs) for x in pairs]

    distlist = run_multiprocessing(pairdist,inps,args.no_cores)

    for pos in range(len(pairs)):
        pair = pairs[pos]
        dist = distlist[pos]
        id1 = pair[0]
        id2 = pair[1]
        newdf.at[id1,id2] = dist
        newdf.at[id2,id1] = dist
        c+=1

        if c%frac == 0:
            print("{}\t\t{}%\t\t{:.3f} seconds".format(c,int((c/tot)*100),time.time() - start_time), end="\r", flush=True)

    print("{}\t\t{}%\t\t{:.3f} seconds".format(c, 100, time.time() - start_time),
          end="\r", flush=True)
    for i in idlist:
        newdf.at[i,i] = 0
    newdf.to_csv(args.distancesout,sep="\t")

    return newdf

def run_agglom_cluster(args,idlist,distancedf,id_to_strain,strain_to_st):
    """
    get cluster assignments for each strain at different cutoffs after running agglomerative clustering
    :param args:
    :param idlist:
    :param distancedf:
    :param id_to_strain:
    :param strain_to_st:
    :return:
    strain_to_cluster = {strain ID:{cluster distance:cluster ID}}
        for snps distance 0 == strain name
        for alleles distance 0 == ST id
    """

    distances = get_distances_frm_args(args)

    clusterlists = {}

    print("\n\nCalculating cluster distance")

    print("distances: ",distances)
    ## run agglomerative clustering (single linkage)
    for dist in distances:
        clusters = AgglomerativeClustering(n_clusters=None,distance_threshold=dist,metric="precomputed", linkage="single").fit_predict(distancedf)
        clusterls = list(clusters)
        clusterlists[dist] = clusterls

    ## assign strains to cluster
    strain_to_cluster = {}

    for i in range(len(idlist)):
        if args.inputtype == "snp":
            strain = idlist[i]
            if strain not in strain_to_cluster:
                strain_to_cluster[strain] = {}
            for d in distances:
                dreal = d-1
                strain_to_cluster[strain][dreal] = int(clusterlists[d][i])
            strain_to_cluster[strain][0] = strain
        elif args.inputtype == "allele":
            strain = idlist[i]
            if strain not in strain_to_cluster:
                strain_to_cluster[strain] = {}
            for d in distances:
                dreal = d-1
                strain_to_cluster[strain][dreal] = int(clusterlists[d][i])
            strain_to_cluster[strain][0] = strain_to_st[strain]

    return strain_to_cluster

def flatten(l):
    return [item for sublist in l for item in sublist]

def make_strains(args):
    """

    :param args: input arguments
    :return:
    strainobjdict - dictionary of {strainID:strainObject}
    odc10_to_strains - dict of {odc10 cluster: strain list}
    """
    straininfodict = {}
    strainobjdict = {}
    odc10_to_strains = {}
    with open(args.straininfo, 'r') as read_obj:
        csv_reader = reader(read_obj,delimiter="\t")
        straininfo = list(csv_reader)
        header = straininfo[0]
        mgtd = OrderedDict()
        odc10col = ""
        if args.enterobase_data:
            for col in range(len(header)):
                colname = header[col]
                if colname == "Name":
                    isolatecol = col
                if colname.startswith("HC"):
                    name = colname.split(" ")[0]
                    # name = "MGT" + colname[2:].replace(" 0 st", "")
                    mgtd[name] = int(col)
                if "Collection Day" in colname:
                    datecol = col
                if "Collection Month" in colname:
                    monthcol = col
                if "Collection Year" in colname:
                    yearcol = col
        else:
            for col in range(len(header)):
                colname = header[col]
                if colname in ["Strain","Isolate"]:
                    isolatecol = col
                if colname == "MGT 1":
                    name = "MGT1"
                    mgtd[name] = col
                if colname == "ODC10":
                    odc10col = col
                if colname.startswith("MGT"):
                    mgtd[colname] = col
                if "Date".lower() in colname.lower():
                    datecol = col
                if "Month".lower() in colname.lower():
                    monthcol = col
                if "Year".lower() in colname.lower():
                    yearcol = col

        for row in straininfo[1:]:
            strainmgtdict = OrderedDict()
            strain = row[isolatecol]
            year = row[yearcol]
            month = row[monthcol]
            if args.enterobase_data:
                day = row[datecol]
                date = "{}-{}-{}".format(year,month,day) ## 2014-02-25
            else:
                date = row[datecol]
            if odc10col != "":
                odc10 = row[odc10col]
                if odc10 not in odc10_to_strains:
                    odc10_to_strains[odc10] = [strain]
                else:
                    odc10_to_strains[odc10].append(strain)
            straininfodict[strain] = row
            missing = False
            for l in mgtd:
                mgtst = straininfodict[strain][mgtd[l]]
                noneg = str(mgtst).split(".")[0]
                if mgtst == "":
                    missing = True
                strainmgtdict[l] = noneg
            timepresent = True
            if args.background_data:
                if year == '' and date == '':
                    timepresent = False
                elif year != '' and date == '' and month == '':
                    month = "01"
                    date = "{}-{}-{}".format(year,month,"01") ## 2014-02-25
                elif year != '' and date == '' and month != '':
                    date = "{}-{}-{}".format(year, month, "01")  ## 2014-02-25

                if not missing and timepresent:
                    strainobjdict[strain] = isolateOb(strain, date, year, month, False, strainmgtdict)
            else:
                if args.timesegment == "week":
                    if year == '' and date == '':
                        timepresent = False
                    elif year != '' and date == '' and month == '':
                        month = "01"
                        date = "{}-{}-{}".format(year, month, "01")  ## 2014-02-25
                        if datetime.strptime(date,'%Y-%m-%d') >= datetime.strptime(args.startdate,'%Y-%m-%d'):
                            timepresent = False
                            sys.exit(
                                f"strain {strain} is in non-background timeframe ({date}) but lacks complete date metadata")
                    elif year != '' and date == '' and month != '':
                        date = "{}-{}-{}".format(year, month, "01")  ## 2014-02-25
                        if datetime.strptime(date,'%Y-%m-%d') >= datetime.strptime(args.startdate,'%Y-%m-%d'):
                            timepresent = False
                            sys.exit(
                                f"strain {strain} is in non-background timeframe ({date}) but lacks complete date metadata")
                elif  args.timesegment == "month":
                    if date == '' and month == '' and year != '':
                        date = "{}-{}".format(year, "01")  ## 2014-02-25
                        if datetime.strptime(date,'%Y-%m') >= datetime.strptime(args.startdate,'%Y-%m'):
                            timepresent = False
                            sys.exit(f"strain {strain} is in non-background timeframe ({date}) but lacks sufficient month and year metadata")
                    elif date == '' and month == '' and year == '':
                        timepresent = False
                        sys.exit(
                            f"strain {strain} is in non-background timeframe ({date}) but lacks sufficient month and year metadata")

                if not missing and timepresent:
                    strainobjdict[strain] = isolateOb(strain,date,year,month,False,strainmgtdict)
    return strainobjdict,odc10_to_strains

def get_newest_strain_date(straindict):
    newest = datetime.strptime('1000-01-01', '%Y-%m-%d')
    for strain in straindict:
        strainob = straindict[strain]
        if strainob.date:
            if strainob.date > newest:
                newest = strainob.date
    print("newest date {}".format(newest))
    return newest

def get_newest_strain_month(straindict):
    newest = datetime.strptime('1000-01', '%Y-%m')
    for strain in straindict:
        strainob = straindict[strain]
        if strainob.month:
            strainmonth = strainob.monthdate
            if strainmonth > newest:
                newest = strainmonth
    print("newest date {}".format(newest))
    return newest

def get_oldest_strain_month(straindict,prevcluster):
    prevstrains = []
    if len(prevcluster.keys()) > 0:
        mincluster = min(map(int,list(prevcluster.keys())))
        for cluster in prevcluster[mincluster]:
            for strain in prevcluster[mincluster][cluster].strains:
                if strain.name not in prevstrains:
                    prevstrains.append(strain.name)
    oldest = datetime.strptime('3000-01', '%Y-%m')
    for strain in straindict:
        if strain not in prevstrains:
            strainob = straindict[strain]
            if strainob.month:
                strainmonth = strainob.monthdate
                if strainmonth < oldest:
                    oldest = strainmonth
    print("oldest date {}".format(oldest))
    return oldest

def get_oldest_strain_date(straindict,prevcluster):
    prevstrains = []
    if len(prevcluster.keys()) > 0:
        mincluster = min(map(int,list(prevcluster.keys())))
        for cluster in prevcluster[mincluster]:
            for strain in prevcluster[mincluster][cluster].strains:
                if strain.name not in prevstrains:
                    prevstrains.append(strain.name)
    oldest = datetime.strptime('3000-01-01', '%Y-%m-%d')
    for strain in straindict:
        if strain not in prevstrains:
            strainob = straindict[strain]
            if strainob.date:
                if strainob.date < oldest:
                    oldest = strainob.date
    print("oldest date {}".format(oldest))
    return oldest

def get_start_end(args,strainobjdict,prevcluster):
    if args.timesegment == "week":
        if args.startdate and args.enddate:
            startdate = datetime.strptime(args.startdate, '%Y-%m-%d')
            enddate = datetime.strptime(args.enddate, '%Y-%m-%d')
        elif args.startdate:
            startdate = datetime.strptime(args.startdate, '%Y-%m-%d')
            enddate = get_newest_strain_date(strainobjdict)
        elif args.enddate:
            enddate = datetime.strptime(args.enddate, '%Y-%m-%d')
            startdate = get_oldest_strain_date(strainobjdict,prevcluster)
        else:
            enddate = get_newest_strain_date(strainobjdict)
            startdate = get_oldest_strain_date(strainobjdict, prevcluster)
    elif args.timesegment == "month":
        if args.startdate and args.enddate:
            startdate = datetime.strptime(args.startdate, '%Y-%m')
            enddate = datetime.strptime(args.enddate, '%Y-%m')
        elif args.startdate:
            startdate = datetime.strptime(args.startdate, '%Y-%m')
            enddate = get_newest_strain_month(strainobjdict)
        elif args.enddate:
            enddate = datetime.strptime(args.enddate, '%Y-%m')
            startdate = get_oldest_strain_month(strainobjdict,prevcluster)
        else:
            enddate = get_newest_strain_month(strainobjdict)
            startdate = get_oldest_strain_month(strainobjdict, prevcluster)
    else:
        sys.exit("timesegment argument should be either week or month")
    return startdate,enddate

def get_background_strains(args,strainobjdict,prevcluster={}):
    startdate, enddate = get_start_end(args, strainobjdict, prevcluster)
    print(startdate,enddate)
    strainls = []
    for strain in strainobjdict:
        strainob = strainobjdict[strain]
        if args.timesegment == "week":
            if strainob.date:
                if strainob.date >= startdate and strainob.date <= enddate:
                    strainls.append(strainob.name)
            else:
                strainls.append(strainob.name)
        elif args.timesegment == "month":
            try:
                if strainob.monthdate:
                    if strainob.monthdate >= startdate and strainob.monthdate <= enddate:
                        strainls.append(strainob.name)
                else:
                    strainls.append(strainob.name)
            except Exception as e:
                print(e)
                print(strainob)
                sys.exit()

    return [("background",strainls)],[]

def add_months(sourcedate, months):
    month = sourcedate.month - 1 + months
    year = sourcedate.year + month // 12
    month = month % 12 + 1
    day = min(sourcedate.day, calendar.monthrange(year,month)[1])
    dt = datetime.strptime("{}-{}-{}".format(year,month,day), "%Y-%m-%d")
    return dt

def get_analysis_groups(args,strainobjdict,prevcluster={}):
    """
    1- get date range of isolates to analyse, from input args or from cluster + strain objects
    2 - get previous strains/strains without date data
    3 - get analysis period and separate isolates into sets
    strainobjdict - dictionary of {strainID:strainObject}
    prevcluster - clusters from previous run, used to determine start time for this analysis
    out
    timesegment_lists = list of tuples with (name of group i.e. date, list of strains included)
    prev_strains = list of strains from before analysis period start
    """
    # 1

    startdate, enddate = get_start_end(args,strainobjdict,prevcluster)

    #2
    prev_strains = []
    for strain in strainobjdict:
        strainob = strainobjdict[strain]
        if args.timesegment == "week":
            if strainob.date:
                if strainob.date < startdate:
                    prev_strains.append(strainob.name)
            else:
                prev_strains.append(strainob.name)
        elif args.timesegment == "month":
            if strainob.monthdate:
                if strainob.monthdate < startdate:
                    prev_strains.append(strainob.name)
            else:
                prev_strains.append(strainob.name)

    #3
    timesegment_lists = []
    start = startdate
    if args.timesegment == "week":
        while start <= enddate:
            subls = []
            end = start + timedelta(days=6)
            for strain in strainobjdict:
                strainob = strainobjdict[strain]
                if strainob.date:
                    if strainob.date >= start and strainob.date <= end:
                        subls.append(strain)
            groupname = "{}_{}".format(start.strftime('%Y-%m-%d'), end.strftime('%Y-%m-%d'))
            timesegment_lists.append((groupname,subls))
            start = start + timedelta(days=7)
    elif args.timesegment == "month":
        while start <= enddate:
            subls = []
            end = add_months(start,1)
            for strain in strainobjdict:
                strainob = strainobjdict[strain]
                if strainob.monthdate:
                    if strainob.monthdate >= start and strainob.monthdate < end:
                        subls.append(strain)
            groupname = start.strftime('%Y-%m')
            timesegment_lists.append((groupname,subls))
            start = add_months(start, 1)
    return timesegment_lists,prev_strains



def make_clusters(args,strain_to_cluster,distance_df,strainobjdict,strain_to_st):
    """

    :param args:
    :param strain_to_cluster:
    :param distance_df:
    :param strainobjdict:
    :param strain_to_st:
    :return:
    clusters = {clusterlev:{clusterid:cluster_instance}}
    """

    distances = get_distances_frm_args(args)
    distances = [0]+[x-1 for x in distances]
    clusters = {}
    for clusterlev in distances:
        clusters[clusterlev] = {}
        for strain in strainobjdict:
            if strain in strain_to_cluster:
                clusterid = strain_to_cluster[strain][clusterlev]
                # if clusterid.isdigit():
                strainobj = strainobjdict[strain]
                if args.timesegment == "week":
                    date_var_to_check = strainobj.date
                elif args.timesegment == "month":
                    date_var_to_check = strainobj.monthdate
                if date_var_to_check != "" and date_var_to_check != "X":
                    if clusterid not in clusters[clusterlev]:

                        c = clusterOb(clusterlev,clusterid)
                        clusters[clusterlev][clusterid] = c

                        c.addstrain(strainobj, strainobj.mgtid_dict, distance_df,strain_to_st)
                        strainobj.partof[clusterlev] = c

                    else:
                        c = clusters[clusterlev][clusterid]

                        c.addstrain(strainobj, strainobj.mgtid_dict, distance_df, strain_to_st)
                        strainobj.partof[clusterlev] = c


    # uncomment below to check for clustering problems
    for clusterlev in distances:
        for clusterid in clusters[clusterlev]:
            clust = clusters[clusterlev][clusterid]

            for strainobj in clust.strains:
                for l in strainobj.partof:
                    if l > clusterlev:
                        if l not in clust.partof:
                            clust.partof[l] = strainobj.partof[l]
                            if len(clust.strains) > len(strainobj.partof[l].strains):
                                print(f"clustering problem1 - strains in cluster {clust} are not recorded in clusters with larger cutoff {strainobj.partof[l]}")
                                sl(1)

                    elif l < clusterlev:
                        strainclust = strainobj.partof[l]
                        if l not in clust.contains:
                            clust.contains[l] = [strainclust]
                            if len(clust.strains) < len(strainclust.strains):
                                print(f"clustering problem2 subcluster {strainclust} within current cluster {clust} has more strains")
                                print(clust)
                                print(strainclust)
                                print(l)
                                print(clusterlev)
                                sl(1)
                        else:
                            if strainclust not in clust.contains[l]:
                                clust.contains[l].append(strainclust)
                                if len(clust.strains) < len(strainclust.strains):
                                    print(f"clustering problem3  sub cluster{strainclust} within current cluster {clust} has more strains in second or later subcluster in same cluster")
                                    sl(1)


    return clusters

def get_distances_frm_args(args):
    """
    parse clustering distances from inputs in args.dist_limits
    :param args:
    :return:
    """
    diststring = args.dist_limits
    dists = diststring.split(",")
    distances = []
    for i in dists:
        if "-" in i:
            n = i.split("-")
            nlist = list(range(int(n[0])+1,int(n[1])+2))
            # distance cutoffs seems to be non inclusive i.e. cutoff of 3 means max distance is 2
            # therefore need to add 1 to all values
            # caused by agglom cluster cutoff setting being pythonic (i.e. 3 == up to 3 == 2. I think)
        else:
            nlist = [int(i)+1]
        distances += nlist
    return distances

def get_smallest_gendist(cluster,smallest):
    if smallest:
        return cluster
    else:
        current_strainls = cluster.strains
        prevclusters = cluster.prevclusters
        if len(prevclusters) > 1:
            return get_smallest_gendist(cluster,True)
        elif prevclusters == []:

            return get_smallest_gendist(cluster, True)
        else:
            prev = prevclusters[0]
            prevstrainls = prev.strains
            if set(current_strainls) == set(prevstrainls):
                return get_smallest_gendist(prev, False)
            else:
                return get_smallest_gendist(cluster, True)

def get_detection_clusters(args,clusters,detectlevel):

    outc = []
    for clusterid in clusters[detectlevel]:
        cluster = clusters[detectlevel][clusterid]
        if cluster.dayspan != "":
            if cluster.dayspan <= args.timewindow and cluster.size >= 2:
                outc.append(cluster)

    return outc


def get_init_clusters(args,clusters,initlevel):
    outc = []

    for clusterid in clusters[initlevel]:
        cluster = clusters[initlevel][clusterid]
        if args.timesegment == "week":
            if cluster.dayspan != "" and cluster.size >= args.minsize:
                    outc.append(cluster)
        elif args.timesegment == "month":
            if cluster.monthspan != "" and cluster.size >= args.minsize:
                    outc.append(cluster)
    return outc

def get_topdown_invest_clusters_week(args,prevclusters,clusters,distances,currpos,out):
    """

    :param args:
    :param prevclusters: initially clusters from the largest distance allowed, subsequently the clusters
    :param clusters: all clusters dictionary
    :param distances: list of distance cutoffs defined in args
    :param currpos: which position in the distance list the current recursion is on
    :param out:
    :param new_clust_in_old_invest:
    :param unchanged_oldclust:
    :return: list of cluster instances assigned as investigation clusters
    """
    print("Processing {} distance clusters".format(currpos))
    currlevclusters = list(clusters[currpos].values())
    passtonext = []
    #for clusters in the currpos-1 distance set
    for prevcluster in prevclusters:
        # if cluster is already an outbreak cluster add to output
        #TODO need to write expanded and unchanged to summary separately
        # for clusters in the currpos distance set
        for cluster in currlevclusters:
            #if the current currpos-1 cluster contains the currpos cluster
            if cluster in prevcluster.contains[currpos]:
                # get all clusters as any distance that contain cluster
                withinclusters = [cluster.partof[x] for x in cluster.partof]
                # if none of the clusters that the current cluster is part of are in investigation clusters output (out)
                if not any([x for x in withinclusters if x in out]):
                    #if prevcluster not stored in prevclusters and distance is not the largest in list
                    if prevcluster not in cluster.prevclusters and currpos != distances[-1]:
                        cluster.prevclusters.append(prevcluster)
                    # if cluster not already in list to process at smaller distance and neither cluster not prevcluster is assigned as investigation
                    if cluster not in passtonext and cluster not in out and prevcluster not in out:
                        # if cluster has an intact timespan
                        if cluster.dayspan != "" and cluster.yearspan != "" and cluster.monthspan != "":
                            # if cluster yearspan and monthspan is <= 2 AND dayspan is <= timewindow in args #TODO paramaterise
                            if cluster.yearspan <= 2 and cluster.monthspan <= 2 and cluster.dayspan <= args.timewindow:
                                #if strain set of current cluster is the same as previous cluster pass to n+1 distance calc
                                if cluster.strains == prevcluster.strains:
                                    passtonext.append(cluster)
                                else:
                                    # if current cluster has a reduced number of strains and prevcluster passes time threshold then prevcluster is investigation cluster #TODO paramaterise
                                    if prevcluster.yearspan < 3 and prevcluster.dayspan <= args.timewindow:
                                        out.append(prevcluster)

                                    else:
                                        # if current distance is the minimum in args then return cluster as investigation
                                        if currpos == distances[0]:
                                            out.append(cluster)
                                        # if current distance is not the minimum in args pass to n+1 distance to check for reduced number of strains
                                        else:
                                            passtonext.append(cluster)
                            # If cluster does not pass time cutoffs pass to next distance cutoff
                            else:
                                passtonext.append(cluster)

    # If currently at minimum distance
    if currpos == distances[0]:
        #confirm cluster in out are above min size # TODO check necessary
        out = [x for x in out if x.size >= args.minsize]
        out2 = []
        for i in out:
            if i not in out2:
                out2.append(i)
        #return investigation clusters
        return out2
    else:
        # run next distance
        return get_topdown_invest_clusters_week(args,passtonext,clusters,distances,currpos-1,out)

def get_topdown_invest_clusters_month(args,prevclusters,clusters,distances,currpos,out):
    print("Processing {} distance clusters".format(currpos))
    currlevclusters = list(clusters[currpos].values())
    cluster_to_check = []
    passtonext = []
    for c in prevclusters:
        containsclusters = c.contains[currpos]
        for cluster in containsclusters:
            partofclusters = [cluster.partof[x] for x in cluster.partof]
            if not any([x in out for x in partofclusters]):
                if cluster in c.contains[currpos]:
                    if c not in cluster.prevclusters and currpos != distances[-1]:
                        cluster.prevclusters.append(c)
                    if cluster.yearspan != "" and cluster.monthspan != "":
                        if cluster.yearspan <= 2 and cluster.monthspan <= args.timewindow:
                            if cluster.strains == c.strains:
                                passtonext.append(cluster)
                            else:
                                if c.yearspan <= 2 and c.monthspan <= args.timewindow and c.size >= args.minsize:
                                    out.append(c)
                                else:
                                    if currpos == distances[0] and cluster.size >= args.minsize:
                                        out.append(cluster)
                                    else:
                                        passtonext.append(cluster)
                        else:
                            passtonext.append(cluster)

    if currpos == distances[0]:
        out = [x for x in out if x.size >= args.minsize]
        out2 = []
        for i in out:
            if i not in out2:
                out2.append(i)
        return out2
    else:
        return get_topdown_invest_clusters_month(args,passtonext,clusters,distances,currpos-1,out)


def ident_additional_isolates(prevcluster,new_strains,idlist):

    investigation_cluster_strains = []
    for distance in prevcluster:
        for cluster in prevcluster[distance]:
            clusterobj = prevcluster[distance][cluster]
            if clusterobj.investigation == True:
                clusterstrainids = [x.name for x in clusterobj.strains]
                investigation_cluster_strains += clusterstrainids
    strainstocluster = investigation_cluster_strains + new_strains
    strainstocluster = [x for x in strainstocluster if x not in idlist]
    strainstoclusterout = []
    [strainstoclusterout.append(x) for x in strainstocluster if x not in strainstoclusterout]
    return strainstoclusterout

def get_remaining_isolates(new_clust_in_old_invest,unchanged_oldclust,toprocess,prevstrains,strains_missing_inputs):
    remove = []
    for level in new_clust_in_old_invest:
        for clusterobj in new_clust_in_old_invest[level]:
            remove += clusterobj.strainnames
    for level in unchanged_oldclust:
        for clusterobj in unchanged_oldclust[level]:
            remove += clusterobj.strainnames
    outids = [x for x in toprocess if x not in remove]
    outids += [x for x in prevstrains if x not in remove]
    outids = [x for x in outids if x not in strains_missing_inputs]
    outids2 = []
    [outids2.append(x) for x in outids if x not in outids2]
    return outids2

def topdown_clust(args,clusters,new_clust_in_old_invest,unchanged_oldclust):
    """

    :param args:
    :param clusters: all clusters from combined dataset
    :param new_clust_in_old_invest:
    :param unchanged_oldclust:
    :return:
    """
    print("Running topdown clustering")
    # get distance settings from args
    distances = get_distances_frm_args(args)
    distances = [0]+[x-1 for x in distances]
    # get initClusters (all clusters at maximum genetic distance above minimum size)
    initClusters = get_init_clusters(args, clusters, distances[-1])

    dateformat = ""
    if args.timesegment == "week":
        dateformat = "%Y-%m-%d"
        # run topdown clustering with number of days as cutoff
        invest_clusters = get_topdown_invest_clusters_week(args, initClusters, clusters, distances, distances[-2], [])
    elif args.timesegment == "month":
        dateformat = "%Y-%m"
        invest_clusters = get_topdown_invest_clusters_month(args, initClusters, clusters, distances, distances[-2], [])

    for i in invest_clusters:
        # check that no clusters within investigation cluster are already marked as investigation
        contains_investigation = False
        for lev in i.contains:
            for containsclust in i.contains[lev]:
                if containsclust.investigation:
                    contains_investigation = True

        if args.timesegment == "week":
            mindate = i.mindate
        elif args.timesegment == "month":
            mindate = i.minmonth
        # mark cluster as investigation if:
        #   doesn't contain investigation cluster
        #   this run is not within the background dataset
        #   The earliest date in this cluster is after the start of the analysis period
        if not contains_investigation and not args.background_data and mindate >= datetime.strptime(args.startdate,dateformat):
            i.investigation = True

    print("Compare to existing clusters")

    invest_out = []
    mgtids = []
    strains_in_old_clusters = []

    new_clust_in_old_invest_clusters = flatten(list(new_clust_in_old_invest.values()))
    unchanged_oldclust_clusters = flatten(list(unchanged_oldclust.values()))
    oldclusters = new_clust_in_old_invest_clusters + unchanged_oldclust_clusters
    for cluster in oldclusters:
        mgtids.append(str(cluster.mgtid))
        strains_in_old_clusters += cluster.strains
    strains_in_old_clusters = set(strains_in_old_clusters)

    for i in invest_clusters:
        #collect any strains that are in existing clusters but are also in new investigation clusters (should be none)
        no_already_clustered = len(set(i.strains).intersection(strains_in_old_clusters))
        if no_already_clustered == 0:
            # assign cluster status as new investigation
            i.status = "new"
            # below section identifies if clusters have an MGT identifier or a simple number (where no MGT assignments are given)
            mgtidsn = []
            for x in mgtids:
                if isinstance(x,int):
                    mgtidsn.append(x)
                else:
                    if x.isdigit():
                        mgtidsn.append(int(x))
            mgtids = list(mgtidsn)
            if mgtids != []:
                nextintid = str(max(mgtids)+1)
            else:
                nextintid = str(1)
            if args.enterobase_data:
                mgtid = i.make_hierccid(nextintid)
            else:
                #if MGT data is present then nextintid is not used, if MGt data is not present the cluster is assigned the name nextintid
                mgtid = i.make_mgtid(nextintid)
            mgtids.append(mgtid)
            # add investigation cluster to list for output
            invest_out.append(i)
    return invest_out

def static_clusters(args, clusters, new_clust_in_old_invest, unchanged_oldclust):
    cutoff = args.static_cutoff
    clusters = clusters[cutoff]
    invest_clusters = [clusters[x] for x in clusters if clusters[x].size >= args.minsize]
    invest_out = []
    mgtids = []
    strains_in_old_clusters = []
    for cluster in list(new_clust_in_old_invest + unchanged_oldclust):
        mgtids.append(cluster.mgtid)
        strains_in_old_clusters += cluster.strains
    strains_in_old_clusters = set(strains_in_old_clusters)

    for i in invest_clusters:
        no_already_clustered = len(set(i.strains).intersection(strains_in_old_clusters))
        if no_already_clustered == 0:
            i.status = "new"
            i.investigation = True
            mgtids = [int(x) for x in mgtids if x.isdigit()]
            if mgtids != []:
                nextintid = max(mgtids) + 1
            else:
                nextintid = 1
            if args.enterobase_data:
                mgtid = i.make_hierccid(nextintid)
            else:
                mgtid = i.make_mgtid(nextintid)
            mgtids.append(mgtid)
            invest_out.append(i)
    return invest_out


def unsupervised_clust():
    sl(1)
    # todo later
def add_uniqueclusterid(strain_to_cluster):
    outres = {}
    for strain in strain_to_cluster:
        outres[strain] = {}
        for dist in strain_to_cluster[strain]:
            if dist == 0:
                outres[strain][dist] = strain_to_cluster[strain][dist]
            else:
                outres[strain][dist] = str(strain_to_cluster[strain][dist])+"_o"
    return outres


def writeout_isolates(args,clusters,strain_2_cluster):
    # output should be "STRAIN   METADATA(time, location etc)    MGT levels  CC-maybe??  11-ODC-levels   investigate?\n
    # strain, ODCs and investigation from clusters
    # MGT and CC from MGT output

    distances = get_distances_frm_args(args)

    distances = [0]+[x-1 for x in distances if x != 1]
    odcslist = [str(x)+'cluster' for x in reversed(distances)]

    straininfodict = {}
    header = ''

    mgtd = OrderedDict()
    ccd = OrderedDict()
    countrycol = ""
    continentcol = ""
    statecol = ""
    postcodecol = ""
    datecol = ""
    monthcol = ""
    yearcol = ""
    sourcecol = ""
    typecol = ""
    hostcol = ""
    diseasecol = ""


    with open(args.straininfo, 'r') as read_obj:
        csv_reader = reader(read_obj,delimiter="\t")
        straininfo = list(csv_reader)
        header = straininfo[0]
        for col in range(len(header)):
            colname = header[col]
            if args.enterobase_data:
                if colname == "Name":
                    isolatecol = col
                if "Collection Day" in colname:
                    datecol = col
                if "Collection Month" in colname:
                    monthcol = col
                if "Collection Year" in colname:
                    yearcol = col
                if colname.startswith("HC"):
                    mgtd[colname] = col
            else:
                if colname in ["Name","Strain","Isolate"]:
                    isolatecol = col
                if "Date".lower() in colname.lower():
                    datecol = col
                if "Month".lower() in colname.lower():
                    monthcol = col
                if "Year".lower() in colname.lower():
                    yearcol = col
                if colname == "MGT 1":
                    name = "MGT1"
                    mgtd[name] = col
                if colname.startswith("MGT"):
                    mgtd[colname] = col
                if colname.startswith("CC-") and "ODC" not in colname:
                    name = colname.split("-")[-1]
                    ccd[name] = col
            if "Country" in colname:
                countrycol = col
            if "Continent" in colname:
                continentcol = col
            if "State" in colname:
                statecol = col
            if "postcode" in colname:
                postcodecol = col

            if "source" in colname:
                sourcecol = col
            if "type" in colname:
                typecol = col
            if "host" in colname:
                hostcol = col
            if "Host disease" in colname:
                diseasecol = col
        metacols = [countrycol,continentcol,statecol,postcodecol,datecol,monthcol,yearcol,sourcecol,typecol,hostcol,diseasecol]
        print(metacols)
        if isolatecol != "":
            for row in straininfo[1:]:
                straininfodict[row[isolatecol]] = row
        else:
            sys.exit("please name column with strain name as one of 'Isolate','Strain','Name'")

    outheader = ["Isolate"]
    for col in metacols:
        if col != "":
            outheader.append(header[col])
    for col in mgtd:
        outheader.append(col)
    for col in ccd:
        outheader.append(col)
    for level in distances:
        outheader.append(str(level)+"cluster")
    outheader.append("investigation cluster")
    # outf.write("\t".join(outheader))
    outrows = [outheader]
    for strain in straininfodict:
        outrow = []
        outrow.append(strain)
        for col in metacols:
            if col != "":
                outrow.append(straininfodict[strain][col])
        for col in mgtd:
            outrow.append(straininfodict[strain][mgtd[col]])
        for col in ccd:
            outrow.append(straininfodict[strain][ccd[col]])
        investigate = ""
        if strain in strain_2_cluster:
            for level in distances:
                id = strain_2_cluster[strain][level]
                if id in clusters[level]:
                    cluster = clusters[level][id]
                    clusterid = str(cluster.id)
                    clustermgtid = str(cluster.mgtid)
                    if cluster.size > 1:
                        outrow.append(clusterid)
                    else:
                        outrow.append("")
                    if cluster.investigation:
                        investigate = "id " + clustermgtid + ":" + str(level)
                else:
                    outrow.append("")
        else:
            for i in distances:
                outrow.append("")
        outrow.append(investigate)
        outrows.append(outrow)

    outdf = pd.DataFrame.from_records(outrows[1:],columns=outrows[0])
    outdf.sort_values(odcslist,inplace=True)
    outdf.to_csv(args.isolates_out,sep="\t",index=False)


def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    ## input files

    required_args = parser.add_argument_group('Required input/output')

    required_args.add_argument("-i", "--variant_data",
                        help="file containing allele profiles (tab delimited table) "
                             "or snp data (wildcard path to snippy outputs "
                             "e.g. /folder/*_snippy = /folder/straina_snippy + /folder/strainb_snippy + ...)"
                             'If using wildcards in path make sure to add ""', required=True)
    required_args.add_argument("--inputtype", help="is input data alleles or snps", choices=["snp","allele"], required=True)


    required_args.add_argument("-s", "--strainmetadata",
                        help="file containing isolate information (downloaded from mgtdb, enterobase or created for SNPs)", required=True)

    required_args.add_argument("--outputPrefix", help="output path and prefix for output file generation",
                        default="output_files", required=True)

    ## optional

    optional_args = parser.add_argument_group('Optional input/output')

    optional_args.add_argument("-d", "--distances",
                        help="file containing pairwise distances corresponding to the alleleprofiles file (from previous run of this script if applicable)")
    optional_args.add_argument("-c", "--inclusters",
                        help="existing clusters to be imported")
    optional_args.add_argument("--background_data",
                        help="data in this input set / time window to be used for background (no outbreak predictions)",action='store_true')
    optional_args.add_argument("-n", "--no_cores",
                               help="number cores to increase pairwise distance speed",default=8,type=int)

    ## snp specific

    snp_args = parser.add_argument_group('SNP input specific')

    snp_args.add_argument("--useref",
                        help="include reference in distances/clusters for snp inputtype",
                        action='store_true')
    snp_args.add_argument("--mask",
                        help="bed file for reference used to generate SNPs with regions to ignore SNPs (i.e. phages etc)")
    snp_args.add_argument("--snpqual",
                        help="minimum allowable SNP quality score",default=1000,type=int)

    ## allele specific

    allele_args = parser.add_argument_group('Allele input specific')

    allele_args.add_argument("--enterobase_data",
                        help="metadata and allele profiles downloaded from enterobase, if hierCC in metadata table hierCC will be used for outbreak naming (i.e. column named HCXXX)",
                        action='store_true')

    ## Time info

    dateranges = parser.add_argument_group('Date / time options')

    dateranges.add_argument("--startdate",
                        help="start date for new cluster analysis (format YYYY-MM-DD) if left blank earliest date not in inclusters will be identified from strain metadata")
    dateranges.add_argument("--enddate",
                        help="end date for new cluster analysis (format YYYY-MM-DD) if left blank latest date in input metadata will be used")
    dateranges.add_argument("--timesegment",
                        help="time segment to perform analysis. every month or every week", choices=["week","month"], default="month")
    dateranges.add_argument("-t", "--timewindow", help="time period a cluster must fall into to be called as investigation"
                                                       " --outbreakmethod topdown or --outbreakmethod twostage only"
                                                       " --timesegment week default 28"
                                                       " --timesegment month default 2", type=int)

    ##cluster options

    cluster_options = parser.add_argument_group('Clustering options')

    cluster_options.add_argument("-l", "--dist_limits",
                        help="comma separated list of cluster cutoffs or range or both i.e 1,2,5 or 1-8 or 1,2,5-10",
                        default="1-5")
    cluster_options.add_argument("-m", "--max_missmatch",
                        help="maximum number of missmatches reported between 2 isolates (will default to max of --dist_limits + 1 if not set)", default=1, type=int)

    ## detection method options

    detection_options = parser.add_argument_group('Outbreak detection algorithm options')

    detection_options.add_argument("--minsize", help="smallest cluster size for outbreak detection", default=5, type=int)
    detection_options.add_argument("--outbreakmethod",
                        help="algorithm for outbreak detection dodge or static",choices=["dodge","static",False],default="dodge")

    detection_options.add_argument("--static_cutoff",
                        help="cutoff for static genetic cutoff method, must be used with '--outbreakmethod static' ",
                        type=int, default=5)


    args = parser.parse_args()

    #################################

    if args.outbreakmethod == "static":
        args.dist_limits = str(args.static_cutoff)
    # print(args.dist_limits)
    distances = get_distances_frm_args(args)
    # print(args.dist_limits)
    # print(distances)
    if args.max_missmatch <= max(distances):
        args.max_missmatch = max(distances)+1

    if args.timesegment == "week":
        if not args.timewindow:
            args.timewindow = 28
    elif args.timesegment == "month":
        if not args.timewindow:
            args.timewindow = 2

    return args



def subsequent_occurrence_indexes(lst):
    indexes = {}
    for index, item in enumerate(lst):
        if item in indexes:
            indexes[item].append(index)
        else:
            indexes[item] = [index]
    return indexes

def main():
    starttime = time.time()

    args = parseargs()
    # Make strain objects that contain MGT type, metadata
    strainobjdict,odc10_to_strains = make_strains(args)
    print(f"make strains --- {time.time() - starttime} seconds ---")
    starttime = time.time()
    prevcluster = {}
    hasdistance = []
    incluster = []
    if args.inclusters:
        print("importing existing clusters")
        # add strain objects to cluster objects
        prevcluster,incluster = import_clusters(args,strainobjdict)

        # get strains that are not already in clusters
        # work out end of previous analysis and start of next
        analysis_groups,prev_strains = get_analysis_groups(args, strainobjdict,prevcluster)
    elif args.background_data:
        #if --background_data set then store all strains in background group
        analysis_groups,prev_strains = get_background_strains(args, strainobjdict)
    else:
        #get sets of strains for each analysis period
        #prevstrains are a list of isolates before the start period of the analysis (background or not)
        analysis_groups,prev_strains = get_analysis_groups(args, strainobjdict)
        #    analysis_groups = list of tuples with (name of group i.e. date, list of strains included)
        #    prev_strains = list of strains from before analysis period start
    print(f"import clusters --- {time.time() - starttime} seconds ---")
    starttime = time.time()
    # import pairwise allele profile difference  matrix (if present) for clustering analysis
    initiald = pd.DataFrame()
    if args.distances:
        initiald = pd.read_csv(args.distances,sep="\t",index_col=0)
        if hasdistance == []:
            hasdistance = list(initiald.head())

        print(f"import distances --- {time.time() - starttime} seconds ---")
    starttime = time.time()

    strains_missing_inputs = [] # list to exclude isolates that have missing input data

    # for each time period defined in analysis groups run full clustering and file output production
    for pos,group in enumerate(analysis_groups):
        print(f"start group {group} --- {time.time() - starttime} seconds ---")
        starttime = time.time()
        groupname = group[0]
        strainls = group[1]

        ## in case where only a single isolate in period, no clustering and add to next periods existing strains
        if len(strainls) == 1 and len(incluster) == 0:
            existstrain = strainls[0]
            next_group = analysis_groups[pos+1]
            next_group[1].append(existstrain)
            continue

        ## Get list of strains that are not in clusters or missing data from existing and new strain lists for period
        inputlist = set(strainls + prev_strains)
        if args.inputtype == "snp":
            strainswdata = check_files_present(args, inputlist)
            strains_missing_inputs = strains_missing_inputs + list(set(inputlist).difference(set(strainswdata)))

        done_and_missing = set(incluster + strains_missing_inputs)
        toprocess = list(inputlist.difference(done_and_missing))

        # if none to process move to next period
        if len(toprocess) == 0:
            continue

        print("\n##############\nStarting analysis for {} on {} isolates".format(groupname,len(toprocess)))

        # output files for current period
        args.distancesout = args.outputPrefix + "_" + groupname + "_pairwise_distances.txt"
        args.clusters_out = args.outputPrefix + "_" + groupname + "_all_clusters.txt"
        args.investigation_out = args.outputPrefix + "_" + groupname + "_investigation_clusters.txt"
        args.isolates_out = args.outputPrefix + "_" + groupname + "_isolate_information.txt"
        print(f"group {group} setup --- {time.time() - starttime} seconds ---")
        starttime = time.time()

        if args.inputtype == "snp":
            # import snippy snp data and make alignment
            idlist, diffdata,missing_inputs = import_snp_data(args,strainls,prev_strains,hasdistance)
            #idlist = list of strains that data was imported for
            #diffdata = dict of {strain:snp_alignment_string}
            # missing_inputs = list of strains not able to be used
            strains_missing_inputs += missing_inputs
            st_to_strain = {}
            strain_to_st = {}
        elif args.inputtype == "allele":
            #import allele data as ...
            idlist, diffdata, st_to_strain, strain_to_st,missing_inputs = import_allele_data(args,strainls,prev_strains,hasdistance)
            strains_missing_inputs += missing_inputs
            #   idlist = list of strains with data successfully gathered
            #   profs = {strain:allele profile for strain as list of numbers as strings}
            #   st_to_strain = {ST:[list of strains assigned ST]}
            #   strain_to_st = {strain:ST}
        else:
            sys.exit("input type is not one of 'snp' or 'allele'")
        print(f"group {group} diffdata import --- {time.time() - starttime} seconds ---")
        starttime = time.time()

        if args.outbreakmethod == "dodge" and not args.background_data:

            starttime = time.time()
            strains_left_to_cluster = ident_additional_isolates(prevcluster,toprocess,strains_missing_inputs)
            # distance_df1 = run_dist(args, diffdata, idlist, hasdistance, initiald, odc10_to_strains,writeout=False)
            distance_df1 = run_dist(args, diffdata, strains_left_to_cluster, hasdistance, initiald, odc10_to_strains)
            print(f"\t run_dist make df subset --- {time.time() - starttime} seconds ---")
            starttime = time.time()
            strain_to_cluster1 = run_agglom_cluster(args, strains_left_to_cluster, distance_df1, st_to_strain, strain_to_st)
            print(f"\t agglom cluster --- {time.time() - starttime} seconds ---")
            starttime = time.time()
            strain_to_cluster2 = add_uniqueclusterid(strain_to_cluster1)
            clusters_existing = make_clusters(args, strain_to_cluster2, distance_df1, strainobjdict, strain_to_st)
            print(f"\t make clusters --- {time.time() - starttime} seconds ---")
            starttime = time.time()
            new_clust_in_old_invest, unchanged_oldclust = ident_expanded_old(prevcluster, clusters_existing,)
            print(f"\t ident_expanded and old --- {time.time() - starttime} seconds ---")
            starttime = time.time()
            # sl(0.1)
            #run above again but with only non investigation isolates and new isolates not already assigned
            remaining_toprocess = get_remaining_isolates(new_clust_in_old_invest,unchanged_oldclust,toprocess,prev_strains,strains_missing_inputs)
            distance_df = run_dist(args, diffdata, remaining_toprocess, hasdistance, initiald, odc10_to_strains)
            print(f"\t run secondary dist --- {time.time() - starttime} seconds ---")
            starttime = time.time()
            strain_to_cluster = run_agglom_cluster(args, remaining_toprocess, distance_df, st_to_strain, strain_to_st)
            print(f"\t run secondary agglom cluster --- {time.time() - starttime} seconds ---")
            starttime = time.time()
            clusters = make_clusters(args, strain_to_cluster, distance_df, strainobjdict, strain_to_st)
            print(f"\t run secondary make clusters --- {time.time() - starttime} seconds ---")
            starttime = time.time()


        else:
            idlist_nomissing = list(set(idlist).difference(set(strains_missing_inputs)))
            # Use diffdata to generate or update a pairwise distance matrix (depending if one was supplied)
            distance_df = run_dist(args, diffdata, idlist_nomissing, hasdistance, initiald, odc10_to_strains)

            # Use distance to generate single linkage clusters
            strain_to_cluster = run_agglom_cluster(args,idlist_nomissing,distance_df,st_to_strain,strain_to_st)
            #     strain_to_cluster = {strain ID:{cluster distance:cluster ID}}
            #         for snps, distance 0 == strain name
            #         for alleles, distance 0 == ST id

            # generates cluster instances and adds strain instances to them returns them as clusters dict
            clusters = make_clusters(args,strain_to_cluster,distance_df,strainobjdict,strain_to_st)
            # clusters = {clusterlev: {clusterid: cluster_instance}}

            # compare new clusters to old ones to identify those that have expanded
            new_clust_in_old_invest,unchanged_oldclust = ident_expanded_old(prevcluster,clusters)
            #new_clust_in_old_invest = list of investigation cluster instances that have gained new isolates
            #unchanged_oldclust = list of unchanged investigation clusters

        # Run Cluster detection
        if args.background_data:
            pass

        elif args.outbreakmethod == "dodge": ## USE THIS
            # Run top down clustering
            investigationClusters = topdown_clust(args,clusters,new_clust_in_old_invest,unchanged_oldclust)
            # investigationClusters = list of clusters assigned as investigation
            print("Writing outputs for {}".format(groupname))
            #output investigation cluster summaries for new, expanded and unchanged investigation clusters
            for invest_clust in investigationClusters:
                if invest_clust.investigation == True:
                    invest_clust.writesummary(args.investigation_out,args)
            for level in new_clust_in_old_invest:
                for expandinvest in new_clust_in_old_invest[level]:
                    investigationClusters.append(expandinvest)
                    expandinvest.writesummary(args.investigation_out,args)
            for level in unchanged_oldclust:
                for unchanged in unchanged_oldclust[level]:
                    investigationClusters.append(unchanged)
                    unchanged.writesummary(args.investigation_out,args)

        elif args.outbreakmethod == "static": ## FOR COMPARISON TO MORE NUANCED METHODS
            # identify new, expanded and unchanged investigation clusters using a static distance cutoff
            investigationClusters = static_clusters(args,clusters,new_clust_in_old_invest,unchanged_oldclust)
            print("Writing outputs for {}".format(groupname))
            for invest_clust in investigationClusters:
                invest_clust.writesummary(args.investigation_out,args)
            for expandinvest in new_clust_in_old_invest:
                expandinvest.writesummary(args.investigation_out,args)
            for unchanged in unchanged_oldclust:
                unchanged.writesummary(args.investigation_out,args)
        #write out all clusters at all levels
        incluster = []


        if args.outbreakmethod == "dodge" and not args.background_data:
            for lev in new_clust_in_old_invest:
                for cluster in new_clust_in_old_invest[lev]:
                    cluster.write(args.clusters_out, args)
                    strains = [x.name for x in cluster.strains if x.name not in incluster]
                    incluster += strains
            for lev in unchanged_oldclust:
                for cluster in unchanged_oldclust[lev]:
                    cluster.write(args.clusters_out, args)
                    strains = [x.name for x in cluster.strains if x.name not in incluster]
                    incluster += strains
            for lev in clusters:
                for id in clusters[lev]:
                    cluster = clusters[lev][id]
                    cluster.write(args.clusters_out, args)
                    strains = [x.name for x in cluster.strains if x.name not in incluster]
                    incluster += strains
        else:
            for lev in clusters:
                for id in clusters[lev]:
                    cluster = clusters[lev][id]
                    cluster.write(args.clusters_out, args)
                    strains = [x.name for x in cluster.strains if x.name not in incluster]
                    incluster += strains

        if args.outbreakmethod == "dodge" and not args.background_data:
            writeout_isolates(args, clusters_existing, strain_to_cluster2)
        else:
            # write out all isolates
            writeout_isolates(args, clusters, strain_to_cluster)
        #update variables for next analysis group based on results from this analysis group
        initiald_head = set(initiald.head())
        distance_head = set(distance_df.head())
        if list(distance_head.difference(initiald_head)) != []:
            initiald = pd.DataFrame(distance_df)

        hasdistance = list(initiald.head())
        prev_strains += idlist
        if args.background_data:
            prevcluster = clusters
        else:
            merged_clusters = {}
            for level in clusters:
                merged_clusters[level] = {}
                for cluster in clusters[level]:
                    cobj = clusters[level][cluster]
                    merged_clusters[level][cluster] = cobj
            if args.outbreakmethod == "dodge":
                for cluster in investigationClusters:
                    merged_clusters[cluster.level][cluster.id] = cluster
                    #TODO utilise/ test below code to report merges better
                # if level in clusters_existing:
                #     for oldcluster in clusters_existing[level]:
                #         if oldcluster not in merged_clusters[level]:
                #             cobj = clusters_existing[level][oldcluster]
                #             merged_clusters[level][oldcluster] = cobj
                #         else:
                #             cobj = clusters_existing[level][oldcluster]
                #             merged_clusters[level][oldcluster + "_2"] = cobj
            prevcluster = merged_clusters

        print("Analysis finished for for {}\n##############\n".format(groupname))


if __name__ == '__main__':
    main()