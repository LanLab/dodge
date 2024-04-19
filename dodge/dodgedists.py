from time import sleep as sl
import sys
import os
import glob
from multiprocessing import Pool
import pandas as pd
import itertools
import time
import argparse

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

def import_allele_data(args,hasdistance):
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

    inprofiles = open(args.variant_data, "r").read().splitlines()


    profs = {}
    st_to_strain = {}
    strain_to_st = {}
    strains = []
    ##per strain
    for line in inprofiles[1:]:
        col = line.split("\t")
        strainid = col[strain]
        if strainid not in profs and col[st] != "":
            noneg = ["0" if x in ['0', "", "-",'""',"''"] else x for x in col[profstart:]]
            noneg = [unneg(x) for x in noneg]
            profs[strainid] = noneg
            st_to_strain[col[st]] = [str(strainid)]
            strains.append(str(strainid))
        else:
            if strainid in profs and col[st] not in st_to_strain:
                print(f"possible duplication of {strainid} in allele profiles")
                continue
            st_to_strain[col[st]].append(str(strainid))
            strains.append(str(strainid))

        strain_to_st[strainid] = col[st]

    if list(set(strains).difference(set(hasdistance))) == [] and len(hasdistance) > 0:
        print("\nAll isolates found in input pairwise distances file")
        st_to_strain, strain_to_st = get_strain_to_st(args, strain, st)
        return strains, {}, st_to_strain, strain_to_st, []

    idlist = []
    [idlist.append(x) for x in profs.keys() if x not in idlist]

    return idlist, profs, st_to_strain, strain_to_st

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

def check_files_present(args):
    inp = args.variant_data + "/*.subs.vcf"
    filels = []
    if "*" in inp:
        filels = glob.glob(inp)

    if len(filels) == 0:
        sys.exit("input path {} found no vcf files".format(inp))

    namels = []
    for file in filels:
        name = os.path.basename(file).replace(".subs.vcf", "")
        namels.append(name)
    return namels

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

def import_snp_data(args,hasdistance):

    print("import_snpdata: import snp start")


    strainswdata = check_files_present(args)
    missing_distance = list(set(strainswdata).difference(set(hasdistance)))
    missing_distance_with_data = list(set(missing_distance).intersection(set(strainswdata)))
    if len(missing_distance_with_data) == 0:
        print("\nAll isolates with data are found in input pairwise distances file")
        return strainswdata, {}, []
    else:
        print(f"to run: {missing_distance_with_data}")
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

def run_multiprocessing(func, i, n_processors):
    with Pool(processes=n_processors) as pool:
        return pool.starmap(func, i)



def pairdist(args,pair,existids,initiald,ap1,ap2):
    id1 = pair[0]
    id2 = pair[1]
    if id1 in existids and id2 in existids:
        dist = int(initiald.at[str(id1), str(id2)])


    else:

        if args.inputtype == "snp":
            dist = snp_dist_metric(ap1, ap2, args)
        elif args.inputtype == "allele":
            dist = allele_dist_metric(ap1, ap2, args)

    return dist

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

def run_dist(args,profs,idlist,hasdistance,initiald):
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

    inps = [(args, x, hasdistance, initiald, profs[x[0]], profs[x[1]]) for x in pairs]

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

    return newdf

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

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # input files

    required_args = parser.add_argument_group('Required input/output')

    required_args.add_argument("-i", "--variant_data",
                        help="file containing allele profiles (tab delimited table) "
                             "or snp data (wildcard path to snippy outputs "
                             "e.g. /folder/*_snippy = /folder/straina_snippy + /folder/strainb_snippy + ...)"
                             'If using wildcards in path make sure to add ""', required=True)
    required_args.add_argument("--inputtype", help="is input data alleles or snps", choices=["snp","allele"], required=True)

    required_args.add_argument("--output", help="output filepath for distance matrix",
                        default="dodge_distance_matrix.txt", required=True)

    ## distances in out

    distance_args = parser.add_argument_group('Optional inputs')

    distance_args.add_argument("-d", "--distances",
                        help="file containing pairwise distances corresponding to the alleleprofiles file (from previous run of this script if applicable)")



    opt_args = parser.add_argument_group('Run options')

    opt_args.add_argument("-n", "--no_cores",
                               help="number cores to increase pairwise distance speed",default=8,type=int)

    opt_args.add_argument("-m", "--max_missmatch",
                        help="maximum number of missmatches reported between 2 isolates ", default=25, type=int)

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
                        help="metadata and allele profiles downloaded from enterobase, if hierCC in metadata table hierCC will be used for outbreak naming (i.e. columns named HCXXX)",
                        action='store_true')




    args = parser.parse_args()

    return args

def main():


    args = parseargs()

    hasdistance = []

    starttime = time.time()
    # import pairwise allele profile difference  matrix (if present) for clustering analysis
    initiald = pd.DataFrame()
    if args.distances:
        initiald = pd.read_csv(args.distances,sep="\t",index_col=0)
        if hasdistance == []:
            hasdistance = list(initiald.head())

        print(f"import distances --- {time.time() - starttime} seconds ---")


    if args.inputtype == "snp":
        # import snippy snp data and make alignment
        idlist, diffdata, missing_inputs = import_snp_data(args, hasdistance)
        # idlist = list of strains that data was imported for
        # diffdata = dict of {strain:snp_alignment_string}
        # missing_inputs = list of strains not able to be used
        st_to_strain = {}
        strain_to_st = {}
    elif args.inputtype == "allele":
        # import allele data as ...
        idlist, diffdata, st_to_strain, strain_to_st = import_allele_data(args, hasdistance)

        #   idlist = list of strains with data successfully gathered
        #   profs = {strain:allele profile for strain as list of numbers as strings}
        #   st_to_strain = {ST:[list of strains assigned ST]}
        #   strain_to_st = {strain:ST}
    else:
        sys.exit("input type is not one of 'snp' or 'allele'")
    print(f"diffdata import --- {time.time() - starttime} seconds ---")
    starttime = time.time()

    distance_df = run_dist(args, diffdata, idlist, hasdistance, initiald)
    
    distance_df.to_csv(args.output,sep="\t")
    print(f"distance matrix complete --- {time.time() - starttime} seconds ---")


if __name__ == '__main__':
    main()