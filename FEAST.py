#!/usr/bin/python

# last update sept 12 2016   @abigailc
# some changes to Shorten/ShortenKeep to remove + = and \
# some comments added.

pfamdatabasedir = "/home/abigail/Downloads/PfamScan/PfamData/"

# edit the above to point to the directory on your machine that contains
# pfam database files (eg pfam-A.hmm) if you plan on running pfam locally
# (not required)


# this script has been cobbled together from around 15 seperate scripts that were created as I was learning python...
# sorry for the huge amount of spaghetti.


# if add taxonomy errors, be sure you are connected to the internet and also have XMLLINT.
# if entire program errors, check you are running python 2.7
# if you are running on 3. it might work if you remove all instances of """ such as the print information bit.
# this script will likely be replaced soon-ish by things using classes and
# persistant storage of seqs.

# hahahaha so many functions_lol why did i not use a reasonable structure.___


# this runs MultiDatasetSubSampling.. all it does it wrap MDSevery inside
# an open file.

def MultiDataSub(inputstr, inlist, indirectory):
        # setting a bunch of variables here, instead of below the parser, for
        # ease of understanding that section
    ins = inputstr.split(",")
    inrank = ins[0]
    innum = ins[1]
    if "strain" in inputstr:
        strain = "yes"
    else:
        strain = "no"
    if "NAdrop" in inputstr:
        NAdrop = "yes"
    else:
        NAdrop = "no"
        # creating and opening the information file
    informationfile = "MultiSub_Info.txt"
    with open(informationfile, "w") as information:

        flist = inlist
        # calling the actual program that will run the multi-dataset
        # subsampler.
        nfl = MDSevery(flist, inrank, innum, strain, NAdrop, information)
    print("MDS done. Info at: " + informationfile)
    return (nfl)


# pfam in looks like PfamOnline("info", fasta) where info == "L, C_term N_term"

def PfamSelect(info, fasta, directory):
    if "," in info:
        pass
    else:
        print("You need to specify the Method of Pfam (O/L/F) then a comma ',' then your domain(s).\neg 'L, FAD_binding' or 'F SQMO_Pfam.txt, SE'")
        print("Try again. exiting now.")
        raise SystemExit
    mode, domains = info.split(",")
    domainlist = domains.split()
    pfop = len(domainlist)
    # run locally
    if mode == "L":
        # avoid fatal error
        print("Hopefully you remembered to set the Pfam Variable as described in the readme/info file, or have the pfamdatabase in your path.\nIf not local pfam will likely not work.")
        if os.path.isfile(str(fileout)) == True:
            print("Cannot run local Pfam -  there is already a file: " + fileout)
            print(
                "Consider using --pfile [your pfam file] instead of running a new scan")
            print("Or delete the offending file, and try again")
            raise SystemExit
        if pfamdatabasedir == "/home/abigail/Downloads/PfamScan/PfamData/":
            print("pfamdatabasedir was not changed; attempting to run without database path specified... outdated!! warning")
            if pfop == 1:
                a = everythinglocal1nodir(fasta, domainlist[0], directory)
            if pfop == 2:
                a = everythinglocal2nodir(
                    fasta, domainlist[0], domainlist[1], directory)
        print("attempting to use pfamdatabasedir variable:" + pfamdatabasedir)
        if pfop == 1:
            a = everythinglocal1(fasta, domainlist[0], directory)
        if pfop == 2:
            a = everythinglocal2(fasta, domainlist[0], domainlist[
                                 1], args.directory)
        # run online
        # you have to upload and save the files that are emailed to you
        # yourself.
    if mode == "O":
        print("Pfam-vetting online!\nThis requires manual input - upload the file(s), and then save the generated pfam files in the indicated manner.\nSorry!")
        if pfop == 1:
            a = everythingnotlocal1(fasta, 4000, domainlist[0], directory)
        if pfop == 2:
            a = everythingnotlocal2(fasta, 4000, domainlist[
                                    0], domainlist[1], directory)
        # run from already generated pfam file.
    if "F" in mode:
        mode, pfile = mode.split()
        if pfop == 1:
            a = GivenPfam1(fasta, pfile, domainlist[0], directory)
        if pfop == 2:
            a = GivenPfam2(fasta, pfile, domainlist[
                           0], domainlist[1], directory)
    print("Pfam exiting")
    return (a)


# tiny wrappers. I don't remember why they exist.
def ShortenKeep(fasta):
    ski = ShortenKeepInfo(fasta)
    return ski


def ShortenNoKeep(fasta):
    skii = Shorten(fasta)
    return skii


# wrapper. this calles extract (no capital) with keepone set to True or
# False. Just a wrapper to get it to work merged in with everything else.
def Extract(directory, fasta, info):
    print(info)
    qualifier = info.split()
    if "keepone" in qualifier:
        qualifier.remove("keepone")
        keepone = True
    else:
        keepone = False
    return(extract(directory, fasta, "Extracted", qualifier, keepone))

# wrapper. this calls seperate with keepone set to True or False. just
# compliance.


def Seperate(directory, fasta, info):
    qualifier = info.split()
    if "keepone" in qualifier:
        qualifier.remove("keepone")
        keepone = True
    else:
        keepone = False
    return (sep(directory, fasta, "Chosen", qualifier, keepone))


# SubSamp
# so input will look like SubSampling(fasta, "inrank innum nadrop")
# im arbitrarily setting NAnum to 1 b/c otherwise you might get a stupid
# amount of things, and nobody got time to set that variable
# if need be, add another variable to be passed through to change NAnum.
# this is also basically a wrapper for the subsampleALL function SSAll.
def SubSampling(fasta, info):
    inf = info.split()
    inrank = inf[0]
    innum = inf[1]
    if "nadrop" in inf:
        nanum = 1
    else:
        nanum = 0
        # period error
    fa, ex = fasta.split(".")

    outputname = fa + "SubSamp" + innum + "p" + inrank + ".fasta"
    return(SSAll(fasta, innum, nanum, inrank, outputname, "no", "$$$$$$$$$no"))

# Merge
   # i want to turn this into it's sub functions of
# 1. merge
# 2. removedupgis
# 3. remove duo taxa
# 4. remove dup AA sequence
##
# could turn last three into one fxn with optional switch eg the pfam function
# like -remove AA/TX/GI
# i like this idea
# *consider making something that removed seqs of very different lengths from the rest#TODO

# methods should be a string input consisting of "AA" "GI" "TX" or a combo
# like "AA GI" or "AA GI TX" for all three operations.

# this calls removeal functions based on what you passed it. wrapper.


def Remove(fasta, methods):
    infile = fasta
    if "GI" in methods:
        infile = RemoveDupGIs(infile)
    if "AA" in methods:
        infile = RemoveDupSeqs(infile)
    if "TX" in methods:
        infile = RemoveDupTaxa(infile)
    if "SI" in methods:
        infile = RemoveDupSeqIDs(infile)
    return (infile)

#... i really couldn't tell you.
# probably this just exists because it used to to something.
# anyways, calls for merge with no sequence removal.


def MergeMerge(inputlist, output):
    return(MergeNoRemove(output, inputlist))


# Summarizes the taxa in a fasta file by rank, eg. "you have 3 bacteria
# and one archaea"
def SumTaxa(fasta, rank):
    countingtaxa = []
    number = int(rank) - 1
    print("got")
    with open(fasta) as old:
        print("here")
        for line in old:
            # print("her")
            if ">" in line:
                # print("gothere")
                alltax = re.sub(
                    "(>)(.*)(\|)(gi#?\|?)([0-9]*)(.*)", "\\2", line)
# print gitax
                if "NT|" in alltax:
                    thetax = "NT"
                else:
                    taxlist = alltax.split("|")
                    thetax = taxlist[number]
                countingtaxa.append(thetax)
##                print (thetax)

    from collections import Counter
    counta = Counter(countingtaxa).most_common()
    print("The following sequences were found:")
    for thing in counta:
        print("%s %s \n" % thing)

# multidatasubsampling

# Looking for shared taxa across several sets of gene sequences.
##
# 1) make a dic of key(rank) value(species name) for each gene eg "fish(scopthalamus rhombus)"
# 2) make list of all ranks for each dataset eg [fish, cat, bird] for SQMO-dataset
# 3) make list of all ranks in all four genes ( rank "fish" present in ubi sqmo shc and osc) 1 option
# 4) make list of all ranks in three (four options)
# 5) make list of all ranks in two (six options)
# 6) make list of all ranks in one (four options)
# 7) print each in nice format
# Orders in all four genes:
# Fish
# scop
# poec
# cicl
# ALSO run sub-sampling on each largest set
##
# when writing this up: we downloaded the full data set for each gene
# we decided to keep one sequence per [order?].
# within each [order?], we chose the taxa that was represented across the largest number of datasets.
# if multiple taxa, we randomly chose one.
##
# >gene|sk|phylum|class|order|species_name|gi|8984932849|blah
# >1|2|3|4|5|sp_name|gi|num|etc
# NOTE: this will keep 1 per order in each dataset, maximizing for
# "shared" species (but always keeping representative breadth)

# makes dictionaries of all the taxa in each dataset.


def MDSfastatodic(fasta, rank, strain="no", na="no", information="na"):
    NAranks = []
    rspec = {}
    information.write("Creating dic for fasta: " + fasta + "\n")
    r = str(int(rank) + 1)
##    replace = "\\"+str(r)+"~\\6"
    cd = os.getcwd()
# print(cd)
    with open(fasta) as original:
        for line in original:
            if ">" in line:
                example = re.sub(
                    "(>)([A-Za-z]*\|)([A-Za-z]*\|)([A-Za-z]*\|)([A-Za-z]*\|)([A-Za-z]*\|)([A-Za-z_0-9]*|)([A-Za-z#|]*)([0-9]*)(.*)", "\\" + r + "~\\7", line)
# print(example)
                try:
                    ra, ss = example.split("~")
                except:
                    print(example)
                    raise SystemExit
                if ra == "NA|":
                    if na == "yes":
                        rn = int(r) + 1
                        srn = str(rn)
                        example = re.sub(
                            "(>)([A-Za-z]*\|)([A-Za-z]*\|)([A-Za-z]*\|)([A-Za-z]*\|)([A-Za-z]*\|)([A-Za-z_0-9]*|)([A-Za-z#|]*)([0-9]*)(.*)", "\\" + srn + "~\\7", line)
                        ra, ss = example.split("~")
                        if ra not in NAranks:
                            NAranks.append(ra)
                ran = ra[:-1]
                sp = ss[:-1]
                if strain == "yes":
                    sp = re.sub("([A-Z][a-z]*_[a-z]*)(.*)", "\\1", sp)
# print(ran)
# print(sp)
                if ran not in rspec:
                    rspec[ran] = []
                rspec[ran].append(sp)
    information.write("NA substitution occurred on " + str(len(NAranks)) +
                      " sequences, resulting in the addition of sub-ranks: \n")
    for naran in NAranks:
        information.write(naran + "\n")
    return rspec

# compares those dictionaries


def MDScomparedics(diclist, sampnum, information):
    listofkeys = []
    dictionaryRtoS = {}
    kn = {}
    import random
    kn = {}
    subslist = []
    # list of all classes
    for d in diclist:
        for key in d:
            if key not in listofkeys:
                listofkeys.append(key)
    for k in listofkeys:
        n = 0
        for d in diclist:
            if k in d:
                n += 1
        if n not in kn:
            kn[n] = []
        kn[n].append(k)
    for number in kn:
        information.write(
            "The following rank(s) are found in " + str(number) + " datasets:\n")
        for rank in kn[number]:
            information.write(rank)
            information.write("\n")
    for k in listofkeys:
        counts = {}
        listofspec = []
        countlist = []
        for d in diclist:
            if k in d:
                for spec in d[k]:
                    if spec in listofspec:
                        pass
                    else:
                        listofspec.append(spec)
        for species in listofspec:
            count = 0
# print(species)

            for d in diclist:
                if k in d:

                    if species in d[k]:
                        count += 1
            if count not in counts:
                counts[count] = []
            counts[count].append(species)
# print(counts)
            countlist.append(count)
        most = max(countlist)
# print(countlist)
        information.write(
            "In rank " + k + ", species in the largest number of datasets (" + str(most) + ") are as follows:")
        information.write("\n")
        randomc = []
        for thing in counts[most]:
            information.write(thing)
            information.write("\n")
            randomc.append(thing)
        if int(sampnum) > len(randomc):

            for each in randomc:
                subslist.append(each)
                dictionaryRtoS[k] = []
                dictionaryRtoS[k].append(each)
                information.write(
                    "The chosen representative(s): " + each + "\n")
        else:
            outp = random.sample(randomc, int(sampnum))
            if "_" not in outp[0]:

                print(outp)
                print("Attempting redraw due to lack of underscore")

                outp = random.sample(randomc, int(sampnum))
                print(outp)
            for t in outp:
                subslist.append(t)
                dictionaryRtoS[k] = []
                dictionaryRtoS[k].append(t)
                information.write("The chosen representative(s): " + t + "\n")
    information.write("Taking " + str(sampnum) +
                      " species per rank.\n\n\n~~~~~~~\n Optimizing for species which are represented in the most datasets.\n Selected species returned, and printed.")
    information.write("\n")
    for s in subslist:
        information.write(s)
        information.write("\n")
##    print (subslist)
    return dictionaryRtoS

# extracts the chosen sequences from each dataset


def MDSextract(filein, outext, qualifier, keepone):
    filename, e = filein.split(".")
    e = "no"
    with open(filein) as old:
        with open(filename + outext + ".fasta", "w") as newe:
            print(filename + outext + ".fasta")
            if type(qualifier) == list:
                for line in old:
                    if ">" in line:
                        # print(line)
                        e = "no"
                        for thin in qualifier:

                            if thin in line:
                                newe.write(line)
                                e = "yes"
                                if keepone == True:
                                    qualifier.remove(thin)
                                continue
                    else:
                        if e == "yes":
                            newe.write(line)
            else:
                for line in old:
                    if ">" in line:
                        e = "no"
                        if qualifier in line:
                            newe.write(line)
                            e = "yes"
                    else:
                        if e == "yes":
                            newe.write(line)

    print("Finished!")

# organizes all of the above MDS functions into a cohesive
# multi-dataset-subsample.


def MDSevery(fastalist, rank, num, strain, nastate, information):
    print("Beginning everything")
    dlist = []
    listoffastas = []
    for fasta in fastalist:
        print("making dictionary from " + fasta)
        a = MDSfastatodic(fasta, rank, strain, nastate, information)
        dlist.append(a)
        listoffastas.append(str(fasta))
    print("making rank to species dictionary")
    dictionaryRtoS = MDScomparedics(dlist, num, information)
    num = 0
    alts = 0
    seqs = 0
    print(listoffastas)
    for fastafile in listoffastas:
        finalspecies = []
        altlist = []
        newflist = []
        print("Creating subsampling from " + fastafile)
        for key in dictionaryRtoS:
            if key in dlist[num]:

                # print(dlist[num][key])
                for speciesitem in dictionaryRtoS[key]:
                    if speciesitem in dlist[num][key]:
                        finalspecies.append(speciesitem)
                        seqs += 1
                    else:

                        finalspecies.append(dlist[num][key][0])
                        seqs += 1
                        alts += 1
                        altlist.append(key + "-" + dlist[num][key][0])
        print("There were " + str(alts) + "replacements made.")
        information.write("\n")
        information.write("There were " + str(alts) + "replacements made: \n")
        for alt in altlist:
            information.write(alt)
            information.write("\n")
        information.write("The overall list for " +
                          fastafile + " selection is: \n")
        towrite = ""
        for aspecin in finalspecies:
            towrite = aspecin + ", " + towrite
        information.write(towrite)
        information.write("\n")
        print("A total of " + str(seqs) +
              " sequences will be written to file...")
        MDSextract(fastafile, "MDS", finalspecies, True)
        newf = fastafile[:-6] + "MDS.fasta"
        newflist.append(newf)
        information.write("A total of " + str(seqs) +
                          " sequences has been be written to file as " + fastafile[:-6] + "MDS.fasta")
        num += 1
        alts = 0
        seqs = 0
    return(newflist)

# vetatree
#!/usr/bin/python

# SHORTEN: takes an input .fasta file (straight from NCBI)
#(filein), and spits out a new
# file named (fileout) with shortened names
#in format >Species_name|gi#|#########


def Shorten(filein):
    print("Shortening sequence IDs from file " + filein)
    import re
    import sys
    try:
        fin, ext = filein.split(".")
        writeto = fin + "Sh.fasta"
        unknownerror = fin + "TaxaError.fasta"
    except:
        writeto = filein + "Sh.fasta"
        unknownerror = filein + "TaxaError.fasta"
    orig = open(filein)
    unk = "no"
    normal = 0
    ucount = 0
    with open(unknownerror, "w") as unknownsp:
        with open(writeto, "w") as new:
            for line in orig:
                if "gi#" in line:
                    print("Fasta IDs already Shortened.")
                    os.system("rm " + writeto)
                    os.system("rm " + unknownerror)
                    return filein
                if ">" in line:
                    # this removes words in brackets that aren't Species_name
                    # and then changes NCBI's default naming scheme to be
                    #>Species_name|gi#|#########
                    # and makes a list of all gi nums and all duplicates
                    number = re.sub(
                        "(>)(gi)(\|)([0-9]*)(\|)([A-Za-z]*)(\|)(.*)(\[\'?[A-Z]?[a-z]* ?.*\])(.*)", "\\4", line)
                    num = number[:-1]
                    edit1 = re.sub(
                        "(>)(gi)(\|)([0-9]*)(\|)([A-Za-z]*)(\|)(.*)(\[\'?[A-Z]?[a-z]* ?.*\])(.*)", "\\1\\9\\3\\2#|\\4", line)
                    if "[" in edit1:
                        unk = "no"
                        normal += 1
                    elif ">gi" in line:
                        edit1 = line
                        unk = "yes"
                        ucount += 1
                    edit2 = re.sub("[\[\]]", "", edit1)
                    edit3 = re.sub("[:;=,/\+'\.\(\)]", "_", edit2)
                    edit4 = re.sub(" ", "_", edit3)
                    edit4 = re.sub("__", "_", edit4)
                    if unk == "no":
                        new.write(edit4)
                    else:
                        unknownsp.write(edit4)
                elif unk == "no":
                    edit = re.sub("[JZ]", "X", line)
                    new.write(edit)
                else:
                    edit = re.sub("[JZ]", "X", line)
                    unknownsp.write(line)
    print("Sequence IDs were shortened, saved as: " + writeto)
    if ucount == 0:
        os.system("rm " + unknownerror)
    if ucount > 0:
        tot = ucount + normal
        print("Of " + str(tot) +
              " sequences, unable to find Species_name in " + str(ucount))
        print("They were removed and saved to " + unknownerror + ".")
    return (writeto)

# same as shorten but keeps more information


def ShortenKeepInfo(filein):
    print("Modifying sequence IDs from file " + filein)
    import re
    import sys
    try:
        fin, ext = filein.split(".")
        writeto = fin + "Mo.fasta"
        unknownerror = fin + "TaxaError.fasta"
    except:
        writeto = filein + "Mo.fasta"
        unknownerror = filein + "TaxaError.fasta"

    print(filein)
    a = os.getcwd()
    print(a)
    orig = open(filein)

    unk = "no"
    normal = 0
    ucount = 0
    with open(unknownerror, "w") as unknownsp:
        with open(writeto, "w") as new:
            for line in orig:
                if "gi#" in line:
                    print("Fasta IDs already Shortened.")
                    os.system("rm " + writeto)
                    os.system("rm " + unknownerror)
                    return filein
                if ">" in line:
                    # this removes words in brackets that aren't Species_name
                    # and then changes NCBI's default naming scheme to be
                    #[Species_name]|gi#"#########" and makes a list of all gi nums and all duplicates
                    number = re.sub(
                        "(>)(gi)(\|)([0-9]*)(\|)([A-Za-z]*)(\|)(.*)(\[\'?[A-Z]?[a-z]* ?.*\])(.*)", "\\4", line)
                    num = number[:-1]
                    edit1 = re.sub(
                        "(>)(gi)(\|)([0-9]*)(\|)([A-Za-z]*)(\|)(.*)(\|)(.*)(\[\'?[A-Z]?[a-z]* ?.*\])(.*)", "\\1\\11\\3\\2#\\3\\4\\3\\10", line)
                    if "[" in edit1:
                        unk = "no"
                        normal += 1
                    elif ">gi" in line:
                        edit1 = line
                        unk = "yes"
                        ucount += 1
                    edit2 = re.sub("[\[\]]", "", edit1)
                    edit3 = re.sub("[:;=,/\+'\.\(\)]", "_", edit2)
                    edit4 = re.sub(" ", "_", edit3)
                    edit4 = re.sub("__", "_", edit4)
                    if unk == "no":
                        new.write(edit4)
                    else:
                        unknownsp.write(edit4)
                elif unk == "no":
                    edit = re.sub("[JZ]", "X", line)
                    new.write(edit)
                else:
                    edit = re.sub("[JZ]", "X", line)
                    unknownsp.write(line)
    print("Sequence IDs were shortened, saved as: " + writeto)
    if ucount == 0:
        os.system("rm " + unknownerror)
    if ucount > 0:
        tot = ucount + normal
        print("Of " + str(tot) +
              " sequences, unable to find Species_name in " + str(ucount))

        print("They were removed and saved to " + unknownerror + ".")
    return (writeto)

# this option keeps only sequences that have a unique GI within the file.
# any that share GI numbers with at least one other sequence are removed
# because they are likely two partially aligned sequences, that can mess
# up the tree.


def RemoveDupGIs(filein):
    import re
    print("Removing all GI-sharing sequences from " + filein)
    try:
        fin, ext = filein.split(".")
        info = fin + "_Info.txt"
        writeto = fin + "_NoDupGIs.fasta"
        writedups = fin + "_ExcludedDupGIs.fasta"
    except:
        info = filein + "_Info.txt"
        writeto = filein + "_NoDupGIs.fasta"
        writedups = filein + "_ExcludedDupGIs.fasta"
    # initializing lists to remove dulplicates
    listofnums = []
    listofdups = []
    listofok = []
    kept = 0
    removed = 0
    dup = "no"
    # finding the gi numbers of all sequences with at least one duplicate
    with open(filein) as original:
        for line in original:
            if ">" in line:
                number = re.sub(
                    "(>)(gi)(\|)([0-9]*)(\|)([A-Za-z]*)(\|)(.*)", "\\4", line)
                num = number[:-1]
                if num in listofnums:
                    ##                    print("Found a duplicate! "+num)
                    listofdups.append(num)
                else:
                    listofnums.append(num)
    # if no duplicates, setting list of all to be ok
    if listofdups == []:
        listofok = listofnums
        print("No Duplicates!")
        return(filein)
    # otherwise making a list of all ok gi numbers.
    else:
        ##        print("Duplicate sequences will be removed and saved in "+writedups)
        for gi in listofnums:
            if gi in listofdups:
                pass
            else:
                listofok.append(gi)
    # dup is on/off switch to see if sequence should be written to good or dup
    # file
    dup = "no"
    orig = open(filein)
    # opening a file to save duplicate, and file to save everything but
    # duplicates
    with open(writedups, "w") as duplicates:
        with open(writeto, "w") as new:
            for line in orig:
                if ">" in line:
                    # finds gi number in NCBI format
                    try:
                        number = re.sub(
                            "(>)(gi)(\|)([0-9]*)(\|)([A-Za-z]*)(\|)(.*)", "\\4", line)
                        num = number[:-1]
                    # finds gi numbers in Shortened format
                    except:
                        try:
                            junk, numb = line.split("gi#|")
                        except:
                            junk, numb = line.split("gi#")
                        num = numb[:-1]
                    if num in listofok:
                        new.write(line)
                        dup = "no"
                        kept += 1
                    elif num in listofdups:
                        dup = "yes"
                        duplicates.write(line)
                        removed += 1
                    else:
                        print("num")
                        print("gi is nowhere error")
                        raise SystemExit
                else:
                    if dup == "yes":
                        duplicates.write(line)
                    elif dup == "no":
                        new.write(line)
                    else:
                        print("Error in writing sequence")
                        raise SystemExit
    import os
    # this ensures that the function passes on a full directory-level name of the file in question
    # probably not necessary inside everything() functions that already go to the dir in question,
    # but it doesn't hurt to keep it.
    # currently OOO b/c was annoying. change the return key to change it.
    filelocation = os.getcwd() + "/" + writeto
    comb = kept + removed
    with open(info, "a") as inform:
        inform.write("\n\nPerformed Identical GI# Sequence Removal.\nOf " + str(comb) + " sequences, " + str(removed) +
                     " were removed for sharing GI#s.\nNew file at " + writeto + ". \nRemoved sequences at " + writedups)
    print("Of " + str(comb) + " sequences, " + str(removed) +
          " had shared GI numbers and were removed.")
    print("New file at " + filelocation +
          ".\n Duplicates (if any) visible at " + writedups)
    # returns the name of the newly created no duplicates file
    return writeto


# RUNPFAM should take an input file in .fasta format, run it through pfam, and create a pfam-output file
# should except the return from shorten(): runpfam(shorten(filein))
# if filein in a different directory, make sure to add
# /home/abigail/Documents/Summons/ or whatever to filein's name
def runpfam(filein):
    import os
    # creates a name for output file
    try:
        fin, ext = filein.split(".")
        fileout = fin + "_Pfam.txt"
    except:
        fileout = filein + "_Pfam.txt"
    print("Attempting to run local Pfam with input = " + filein)
    # this runs the pfam_scan script from terminal, using perl.
    os.system("pfam_scan.pl -fasta " + filein + " -dir " +
              pfamdatabasedir + " -outfile " + fileout)
    # something to test if the output was actually created correctly
    if os.path.isfile(fileout) == True:
        print("Success! The export file is " + fileout)
    else:
        print("looked in: " + fileout)
        print("no file here?")
        raise SystemExit
    return fileout


def runpfamnodir(filein):
    import os
    # creates a name for output file
    try:
        fin, ext = filein.split(".")
        fileout = fin + "_Pfam.txt"
    except:
        fileout = filein + "_Pfam.txt"
    print("Attempting to run local Pfam with input = " + filein)
    # this runs the pfam_scan script from terminal, using perl.
    os.system("pfam_scan.pl -fasta " + filein + " -outfile " + fileout)
    # something to test if the output was actually created correctly
    if os.path.isfile(fileout) == True:
        print("Success! The export file is " + fileout)
    else:
        print("looked in: " + fileout)
        print("no file here?")
        raise SystemExit
    return fileout

# PFAMSPLT takes a file of shortened names, and splits it into chunks manageable by pfam server online
# for now, just splits (FILEIN) into chunks of (SIZE) sequences, saved in sequential files
# eg test.fasta splits and saves as test1.fasta, test2.fasta etc.
# NOT included in "everything"


def pfamsplit(filein, size):
    # opens the file, lets you know its running
    print("Attempting to make files of length " + str(size) + " sequences")
    # tick counts how many sequences have been written to a single file
    import re
    tick = 0
    # num represents the number of files that have been created, starting at 1
    num = 1
    # this ensures that whatever you put it, will put out .fasta files
    try:
        fin, ext = filein.split(".")
        writeto = fin + "1.fasta"
    except:
        writeto = filein + "1.fasta"

    # So that you know how many files you will need to upload to pfam server
    # online, and what they are called.
    print("Writing file to: " + writeto)
    with open(filein) as original:
        with open(writeto, "w") as new:
            for line in original:
                if ">" in line:
                    tick += 1
                # specifies when and how to make a new file
                if int(tick) == int(size):
                    num += 1
                    try:
                        writeto = str(fin) + str(num) + ".fasta"
                    except:
                        writeto = filein + str(num) + ".fasta"
                    print("Writing file to: " + writeto)
                    new = open(writeto, "w")
                    tick = 0
                nodash = re.sub("-", "", line)
                new.write(nodash)
            new.flush()
            new.close()
# print("Done")
    # instructions for running the split files through pfam online
    print("Please take the above files and run through Pfam(http://pfam.xfam.org/search#tabview=tab1)")
    print("When finished (may take a while), copy the information (in gmail, use \"see original\") into a .txt document")
    print("Save it as filein#Pfam.fasta.")
    # only say "y" if it is done and you saved it correctly.
    a = raw_input("Have you saved all Pfam files? Type y if so ")
    print (a)
    if a != "y":
        print("That's not valid input")
        a = raw_input("Type y to continue, or anything else to exit.")
        print(a)
    if a == "y":
        try:
            a = fin + "Pfam.fasta"
        except:
            a = filein + "Pfam.fasta"
        return a
    else:
        # tells you what to do if you messed up/decided to exit
        print ("Run Pfam from a pre-gen file using the F option.\nIf your file was split, you will now have to manually combine.")
        quit()


# COMBINEPFAM to be used after pfam split/manual save. NOT included in "everything".
# combines multiple files of type name1Pfam.txt into one large file
# returns the name of the large file
# as input, give it "namePfam.fasta" (no numbers, though all your
# pfam-output files SHOULD have numbers (name1Pfam.fasta)

def combinepfam(pfamfile):
    print("Starting combine")
    import os
    i = 1
    # creates name to open first pfam-output file.
    try:
        fin, ext = pfamfile.split(".")
        checkfile = fin[:-4] + str(i) + "Pfam.fasta"
    except:
        checkfile = pfamfile[:-4] + str(i) + "Pfam.fasta"
    # creates output combined pfam file name.
    newfile = pfamfile[:-10] + "PfamCombined.fasta"
    # opens the new file to write to
    with open(newfile, "w") as new:
        print("looking for: " + checkfile)
        # this will be false when it tries to open, say, the fifth file but you
        # only have four.
        while os.path.isfile(checkfile) == True:
            print("Copying " + checkfile)
            check = open(checkfile)
            for line in check:
                new.write(line)
            # this makes it create a new file to look for of increment one
            # higher than the last it found.
            i += 1
            try:
                checkfile = fin[:-4] + str(i) + "Pfam.fasta"
            except:
                checkfile = pfamfile[:-4] + str(i) + "Pfam.fasta"
    print("The pfam files have been combined and saved as " + newfile)
    return newfile

# SIMPLIFYPFAM takes a pfam file, removes all the junk at the top/bottom, and uses regular
# expressions to reformat it to keep just the sequence ID and identified
# domain.


def simplifypfam(pfamfile):
    import re
    original = open(pfamfile)
    try:
        fin, ext = pfamfile.split(".")
        outfile = fin + "Simp.fasta"
    except:
        outfile = pfamfile + "Simp.fasta"
    print("Simplifying Pfam file " + pfamfile)
# editing the pfam document with regex to leave only id and domain name
    with open(outfile, "w") as new:
        ok = "no"
        for line in original:

            if ok == "yes":

                edit = re.sub(
                    "(.*)(\t)([0-9][0-9]*\t)([0-9][0-9]*\t)([0-9][0-9]*\t)([0-9][0-9]*\t)([A-Z.0-9]*\t)(.*)(\t)([A-Z][a-z]*\t)(.*)", "\\1 \\8", line)
            else:
                edit = ""
            # changed this to match the proper regex to leave only ID and
            # Domain
            new.write(edit)
            if "<clan>" in line:
                ok = "yes"
##    print ("Simplified and saved as "+outfile)
    return outfile


# VERIFYDOMAIN takes input a simplified pfam file + domain to look for, and gives a
# list of GI numbers that
# have that domain.
def verifydomain(inputfile, inputdomain):
    # opens file and initiates list
    print ("Now making a list of all sequences in " +
           inputfile + " that have domain " + inputdomain)
    orig = open(inputfile)
    listofinputdomains = inputdomain.split()
    # makes empty list
    nameswithdomain = []
    # line-by-line search - if has proper domain, adds gi num to list
    for line in orig:
        name, domain = line.split(" ")

        for d in listofinputdomains:

            if str(domain[:-1]) in str(d):
                if "gi#|" in name:

                    junk, number = name.split("gi#|")
                if "gi#" in name:
                    junk, number = name.split("gi#")
                else:
                    number = name

                nameswithdomain.append(number)
    # returns list of gi nums that have the domain
    if nameswithdomain == []:
        print("Error, no sequences found with domain: " + inputdomain)
    else:
        return nameswithdomain

# gets a list of all GI, used to generate "vetted-out"
# def getallgi(filein):
##    a = open(filein)
# for line in filein:
# if ">" in line:
# junk, number = name.split("gi#")
##            num = number[:-1]
# numlist.append(num)
# return numlist


# VETSEQS : combines several functions above:
    # give it domains to look for and a pfam file, and it will run SIMPIFY, run VERIFYDOMAIN
    # for each domain provided, and make a list of all gi#'s that have both
    # domains.
def vetseqs(spfam, domain1, domain2):
    # runs the simplify function, and sets a to the name of the new file
    hasD1 = verifydomain(spfam, domain1)
    hasD2 = verifydomain(spfam, domain2)
    hasBoth = []
    for number in hasD1:
        if number in hasD2:
            hasBoth.append(number)
    return hasBoth


def pfamdomains(simpfam):
    print("Creating a list of domain frequencies")
    listofdomains = []
    with open(simpfam) as a:
        for line in a:
            name, dom = line.split(" ")
            domain = dom[:-1]
            listofdomains.append(domain)
    from collections import Counter
    return Counter(listofdomains).most_common()


def VetWrite(shfasta, hasDom, spfam, indomains):
    print("debug: vetwrite")
    domainnums = pfamdomains(spfam)
    if str(hasDom) == "None":
        print ("No sequence has the indicated domains. Domains found are:")
        for thing in domainnums:
            print("%s %s \n" % thing)
        raise SystemExit
    try:
        fin, ext = shfasta.split(".")
        outputfile = fin + "_Vet.fasta"
        outbad = fin + "_VetRemoved.fasta"
        info = fin + "_Info.txt"
    except:
        filein = shfasta
        outputfile = filein + "_Vet.fasta"
        outbad = filein + "_VetRemoved.fasta"
        info = filein + "_Info.txt"
    origfasta = open(shfasta)
    vout = 0
    keep = 0
    removeddomains = []
##    print (hasDom)
    with open(outbad, "w") as bad:
        with open(outputfile, "w") as new:
            write = 0
            for line in origfasta:
                if ">" in line:
                    write = 0
                    if "gi#" in line:
                        # print("gi")
                        try:
                            junk, num = line.split("gi#|")
                            if "|" in num:
                                newn, blah = num.split("|")
                            else:
                                newn = num

                        except:
                            junk, num = line.split("gi#")
                            newn = num[:-1]
                    else:
                        newn = line[1:]

                    if "\n" in newn:
                        newn = newn[:-1]
                    a = []
                    a.append(newn)
# print(a)
                    if newn in hasDom:
                        ##                        print ("yes")
                        keep += 1
                        write = 1
                        new.write(line)
                    else:
                        vout += 1
                        bad.write(line)
# print("wrote"+line)
                        # adds removed domain-names to removed document
                        with open(spfam) as shortpfam:
                            for seq in shortpfam:
                                if newn in seq:
                                    # print("found")
                                    name, dom = seq.split(" ")
                                    domain = dom[:-1]
    # print("adding"+domain)
                                    bad.write(domain + "\n")
                                    removeddomains.append(domain)
                elif write == 1:
                    new.write(line)
                elif write == 0:
                    bad.write(line)
            new.flush()
            new.close()
    totseq = vout + keep
    percent = (float(vout) / float(totseq)) * 100
    perc = '%.2f' % percent
    with open(info, "w") as information:
        information.write("You Vetted for any of: " + indomains + "\n\n")
        information.write("Of " + str(totseq) + " total sequences, " +
                          str(vout) + " of them were removed during vetting. \n")
        information.write("That is " + perc + "% of all sequences. \n\n")
        information.write(
            "Here is a list of all domains that were found, and their frequency: \n")
        for thing in domainnums:
            information.write("%s %s \n" % thing)
        information.write(
            "\n Here is a list of all domains found on sequences that were removed: \n")
        from collections import Counter
        count = Counter(removeddomains).most_common()
        for item in count:
            information.write("%s %s \n" % item)
    print("You Vetted for any of: " + indomains)
    print("Of " + str(totseq) + " total sequences, " +
          str(vout) + " of them were removed during vetting.")
    print("That is " + perc + "% of all sequences.")
    print("Here is a list of top five domains found, and their frequency:")
    c = 0
    for thing in domainnums:
        if c > 4:
            pass
        else:
            print("%s %s" % thing)
            c += 1
    print("Here is a list of the top five domains on all removed sequences: ")
    from collections import Counter
    count = Counter(removeddomains).most_common(5)
    for item in count:
        print("%s %s" % item)
    import os
    print ("Done, Vetted sequences exported in " + outputfile)
    print ("Anything removed exported as " + outbad)
    print ("Information regarding this Vet available at " + info)
    return outputfile


#######WAY MORE WRAPPERS THAT AT SOME POINT DEALT WITH SOME SHIT LOL#####

# SHORTPFAMVET will shorten, runpfam, vetseqs, make .fasta with two
# verified domains
def ShortPfamVet(shfasta, domain1, domain2):

    bothdom = domain1 + " AND any of: " + domain2
    spfam = simplifypfam(runpfam(shfasta))
    hasBoth = vetseqs(spfam, domain1, domain2)
    return VetWrite(shfasta, hasBoth, spfam, bothdom)

# SHORTPFAMVET will shorten, runpfam, vetseqs, make .fasta with two
# verified domains


def ShortPfamVetnodir(shfasta, domain1, domain2):

    bothdom = domain1 + " AND any of: " + domain2
    spfam = simplifypfam(runpfamnodir(shfasta))
    hasBoth = vetseqs(spfam, domain1, domain2)
    return VetWrite(shfasta, hasBoth, spfam, bothdom)


# SINGLEVET will shorten, run pfam, vetseqs, make .fasta with one verified
# domain.
def SingleVet(shfasta, domain1):
    spfam = simplifypfam(runpfam(shfasta))
    hasD1 = verifydomain(spfam, domain1)
    return VetWrite(shfasta, hasD1, spfam, domain1)


def SingleVetnodir(shfasta, domain1):
    spfam = simplifypfam(runpfamnodir(shfasta))
    hasD1 = verifydomain(spfam, domain1)
    return VetWrite(shfasta, hasD1, spfam, domain1)

# EverythingLocal will run ShortPfamVet, but from the correct specified
# directory to keep all your generated files together
# def everythinglocal(filein, domain1, domain2, directory):
# print "Running Everything with Local Pfam"
##    import os
# os.chdir(directory)
# return ShortPfamVet(filein, domain1, domain2)

# EVERYTHINGLOCAL2: does everything, vets for 2 domains.


def everythinglocal2(shfasta, domain1, domain2, directory):
    # print "Running Everything with Local Pfam, two domains"
    import os
    os.chdir(directory)
    return ShortPfamVet(shfasta, domain1, domain2)

# EVERYTHINGLOCAL1: does everything, vets for 1 domain.


def everythinglocal1(shfasta, domain1, directory):
    # print "Running Everything with Local Pfam, one domain"
    import os
    os.chdir(directory)
    return SingleVet(shfasta, domain1)


def everythinglocal2nodir(shfasta, domain1, domain2, directory):
    # print "Running Everything with Local Pfam, two domains"
    import os
    os.chdir(directory)
    return ShortPfamVetnodir(shfasta, domain1, domain2)

# EVERYTHINGLOCAL1: does everything, vets for 1 domain.


def everythinglocal1nodir(shfasta, domain1, directory):
    # print "Running Everything with Local Pfam, one domain"
    import os
    os.chdir(directory)
    return SingleVetnodir(shfasta, domain1)

# SHORTSPLIT goes to directory, and runs shorten and split on a fasta file.


def shortsplit(filein, num, directory):
    print ("Running Shorten and Split")
    import os
    os.chdir(directory)
    pfamsplit(Shorten(filein), num)


# SHORTEN2 specifies a directory, changes into that directory, and then
# runs shorten.
def Shorten2(filein, directory):
    import os
    os.chdir(directory)
    return Shorten(filein)
# CombVet is what you do after you have a short fasta file and saved pfam files from online,
# and need to run all the vetting stuff.
# Idea is 1. SHORTSPLIT 2.manual upload/saving 3.COMBVET (will be combined in later fxn)
# fasta is the shortened fasta file you should have saved
# pfam is the pfam file you saved, with no numbers.
# domains are the strings of domains you want to match
# directory is the directory in which you fasta and pfam files are saved in.


def CombVet2(shfasta, pfam, domain1, domain2):
    print ("Running Combine and Vet 2D")
    spfam = simplifypfam(combinepfam(pfam))
    hasBoth = vetseqs(spfam, domain1, domain2)
    bothdom = domain1 + " AND any of: " + domain2
    return VetWrite(shfasta, hasBoth, spfam, bothdom)


def CombVet1(shfasta, pfam, domain1):
    print ("Running Combine and Vet 1 Domain")
    spfam = simplifypfam(combinepfam(pfam))
    hasD1 = verifydomain(spfam, domain1)
    return VetWrite(shfasta, hasD1, spfam, domain1)

# EVERYTHINGNOTLOCAL runs everything. goes to directory to keep stuff together, shortens,
# splits, allows for manual pfam, joins pfam files, vets for two domains,
# spits out a vetted .fasta file
# only use if you want to have vetting, but do not have pfam installed on
# your computer.


def everythingnotlocal2(shfasta, splitnum, domain1, domain2, directory):
    import os
    os.chdir(directory)
    return CombVet2(shfasta, pfamsplit(shfasta, splitnum), domain1, domain2)


def everythingnotlocal1(shfasta, splitnum, domain1, directory):
    import os
    os.chdir(directory)
    return CombVet1(shfasta, pfamsplit(shfasta, splitnum), domain1)


def GivenPfam2(shfasta, pfam, domain1, domain2, directory):
    import os
    os.chdir(directory)
    spfam = simplifypfam(pfam)
    hasBoth = vetseqs(spfam, domain1, domain2)
    bothdom = domain1 + " AND any of: " + domain2
    return VetWrite(shfasta, hasBoth, spfam, bothdom)


def GivenPfam1(shfasta, pfam, domain1, directory):
    import os
    os.chdir(directory)
    spfam = simplifypfam(pfam)
    hasBoth = verifydomain(spfam, domain1)
    return VetWrite(shfasta, hasBoth, spfam, domain1)

# AppendRank
# gi to taxid
##


def CountSeqs(fasta):
    inc = 0
    with open(fasta) as tocount:
        for line in tocount:
            if ">" in line:
                inc += 1
    return inc

# this is what fetches taxonomy from NCBI. requires XMLLINT and an
# internet connection.


def AddRank(fasta, ranklist, dire):
    print ("Starting on file: " + fasta)
    import re
    import subprocess

    nam, ext = fasta.split(".")
    namer = ""
    for r in ranklist:
        R = str.capitalize(r)
        namer = namer + R[:1]
    import os
    print(os.getcwd())
    newfasta = nam + namer + ".fasta"
    totallines = CountSeqs(fasta)
    with open(fasta) as old:
        with open(newfasta, "w") as new:
            inc = 0
            perklist = []
            for line in old:
                perc = int((inc % totallines) * 100)
                if perc % 5 == 0:
                    if str(perc) in perklist:
                        pass
                    else:
                        print(str(perc) + " percent done... working")
                        perklist.append(str(perc))
                if "gi#|" in line:
                    ginum = re.sub(
                        "(>)(.*)(\|)(gi#\|)([0-9]*)(\|?)(.*)", "\\5", line)
                    do = "yes"
                elif "gi" in line:
                    ginum = re.sub(
                        "(>)(.*)(\|)(gi#?[\|]?)([0-9]*)(\|?)(.*)", "\\5", line)
                    do = "yes"
##                    print (ginum)
# print("GOTHERE")
                else:
                    new.write(line)
                    do = "no"
                if do == "yes":
                    fullrank = ""
# print(ginum)
# print(ginum[:-1])
                    futuretaxid = "error"
                    ginum = ginum[:-1]
                    GItoTAXID = "xmllint --xpath '/GBSet/GBSeq/GBSeq_feature-table/GBFeature/GBFeature_quals/GBQualifier[GBQualifier_name=\"db_xref\"]/GBQualifier_value/text()' \"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=" + \
                        ginum + "&retmode=xml\""
# print(GItoTAXID)
##                    GItoTAXID = "xmllint --xpath '/GBSet/GBSeq/GBSeq_feature-table/GFeature_quals/GBQualifier[GBQualifier_name=\"db_xref\"]/GBQualifier_value/text()' \"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id="+ginum+"&retmode=xml\""
                    try:
                        futuretaxid = subprocess.check_output(
                            GItoTAXID, shell=True)
                    except:
                        print("Error: Cannot find Taxid for GI num: " +
                              ginum + ", skipping")
                        tosub = "\\1NT\\3\\2\\3\\4\\5\\6\\7"
                        newline = re.sub(
                            "(>)(.*)(\|)(gi#?[\|]?)([0-9]*)(\|?)(.*)", tosub, line)
                        new.write(newline)
                        continue
                    taxid = re.sub("(taxon:)([0-9]*)(.*)", "\\2", futuretaxid)
# print(taxid)
                    for r in ranklist:
                        if r == "CommonName":
                            TAXIDtoRANKNAME = "xmllint --xpath '/TaxaSet/Taxon/OtherNames/GenbankCommonName/text()'  \"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=" + taxid + "\""
                        else:
                            TAXIDtoRANKNAME = "xmllint --xpath '/TaxaSet/Taxon/LineageEx/Taxon[Rank=\"" + r + \
                                "\"]/ScientificName/text()'  \"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=" + taxid + "\""
                        try:
                            rankname = subprocess.check_output(
                                TAXIDtoRANKNAME, shell=True)
                        except:
                            rankname = "NA"
    ##                    print (rankname)
                        rankname = re.sub(" ", "_", rankname)
                        fullrank = rankname + "|" + fullrank

                    tosub = "\\1" + fullrank + "\\2\\3\\4\\5\\6\\7"
                    newline = re.sub(
                        "(>)(.*)(\|)(gi#?[\|]?)([0-9]*)(\|?)(.*)", tosub, line)
                    new.write(newline)
    print("Done! Added ranks. File at " + newfasta)
    return newfasta


def extract(direct, filein, outext, qualifier, keepone):
    import os
    os.chdir(direct)
    print("running in" + direct)
    filename, e = filein.split(".")
    e = "no"
    with open(filein) as old:
        with open(filename + outext + ".fasta", "w") as newe:
            print(filename + outext + ".fasta")
            if type(qualifier) == list:
                for line in old:
                    if ">" in line:
                        # print(line)
                        e = "no"
                        for thin in qualifier:

                            if thin in line:
                                newe.write(line)
                                e = "yes"
                                if keepone == True:
                                    qualifier.remove(thin)
                                    print(thin)

                                continue
                    else:
                        if e == "yes":
                            newe.write(line)
            else:
                for line in old:
                    if ">" in line:
                        e = "no"
                        if qualifier in line:
                            newe.write(line)
                            e = "yes"
                    else:
                        if e == "yes":
                            newe.write(line)

    print("Finished!")
    return(filename + outext + ".fasta")

# seperated sequences with "Qualifier" in seqID from those without it.


def sep(direct, filein, outext, qualifier, keepone):
    import os
    os.chdir(direct)
    filename, e = filein.split(".")
    e = "no"
    with open(filein) as old:
        with open(filename + outext + ".fasta", "w") as newe:
            with open(filename + "Not" + outext + ".fasta", "w") as newn:
                if type(qualifier) == list:

                    for line in old:
                        if ">" in line:
                            e = "no"

                            for thin in qualifier:

                                if thin in line:
                                    newe.write(line)
                                    e = "yes"
                                    if keepone == True:
                                        qualifier.remove(thin)
                                        print("removed" + thin)

                                    continue
                            if e == "no":
                                newn.write(line)

                        else:
                            if e == "yes":
                                newe.write(line)
                            if e == "no":
                                newn.write(line)
                else:
                    for line in old:
                        if ">" in line:
                            e = "no"
                            n = "no"
                            if qualifier in line:
                                newe.write(line)
                                e = "yes"
                                n = "no"
                            else:
                                newn.write(line)
                                e = "no"
                                n = "yes"
                        else:
                            if e == "yes":
                                newe.write(line)
                            if e == "no":
                                if n == "yes":
                                    newn.write(line)

        print("Finished!")
        return(filename + outext + ".fasta")
# SubSample

# uhh this is used in something. probably.


def MakeLists(fasta, rank):
        # period error
    ##    n = fasta
    n, ex = fasta.split(".")
    import re
    number = int(rank)
    diTRL = {}
    diNA = {}
    listofERROR = []
    with open(fasta) as old:
        for line in old:
            if ">" in line:
                if "gi#" in line:
                    gitax = re.sub(
                        "(>)(.*)(\|)(gi#?\|?)([0-9]*)(.*)", "\\5~\\2", line)
        # print gitax
                    gi, tax = gitax.split("~")
                    tn = tax[:-1]
                    taxlist = tn.split("|")
        ##                print (taxlist)
                    if taxlist[0] == "NT":
                        listofERROR.append(gi)
                        continue
                    thetax = taxlist[number - 1]
                else:
                    taxlist = line.split("|")
        ##                print (taxlist)
                    if taxlist[0] == "NT":
                        listofERROR.append(gi)
                        continue
                    thetax = taxlist[number - 1]
                    gi = line
                if thetax == "NA":
                    try:
                        thetax = taxlist[number]
                    except:
                        pass
                    if thetax == "NA":
                        try:
                            thetax = taxlist[number + 1]
                        except:
                            pass
                        if thetax == "NA":
                            listofERROR.append(gi)
                            continue
                        if thetax in diNA:
                            # print(diNA[thetax])
                            diNA[thetax].append(gi)
                        else:
                            diNA[thetax] = [gi]
                elif thetax in diTRL:
                    diTRL[thetax].append(gi)
                else:
                    diTRL[thetax] = [gi]

    print("Created lists")
# print(diNA)
    return diTRL, diNA, listofERROR


def SeqsToKeep(makelistsout, number, NAnum):
    diTRL = makelistsout[0]
    number = int(number)
    diNA = makelistsout[1]
    import random
    keep = []
    for entry in diTRL:
        try:
            r = random.sample(diTRL[entry], number)
        except ValueError:
            le = len(diTRL[entry])
            r = random.sample(diTRL[entry], le)
        for thing in r:
            keep.append(thing)
##    print (keep)
    for entry in diNA:
        try:
            r = random.sample(diNA[entry], NAnum)
        except ValueError:
            le = len(diNA[entry])

            r = random.sample(diNA[entry], le)
        for thing in r:
            keep.append(thing)
    print("Picked a random subset of each group")
    return keep


def GetSeqs(fasta, keep, output, special):
    copy = "no"
    import re

    with open(fasta) as old:
        with open(output, "w") as new:
            for line in old:
                if ">" in line:
                    gi = re.sub(
                        "(>)(.*)(\|)(gi#?\|?)([0-9]*)(.*)(\n)", "\\5", line)
                    if gi in keep:
                        copy = "yes"
                        new.write(line)
                        # to stop it from copying multiple sequences in the
                        # case of shared gi numbers All(whole genome seq etc)
                        keep.remove(gi)
                    elif special in line:
                        copy = "yes"
                        new.write(line)
                    else:
                        copy = "no"
                else:
                    if copy == "yes":
                        new.write(line)
            new.flush()
    print("Created a subset file at " + output)
    return output


def ShowErrors(fasta, output, erlist):
    import re
    n, e = output.split(".")
    num = 0
    errorfa = n + "_SubErr.txt"
    with open(fasta) as old:
        with open(errorfa, "w") as new:
            for line in old:
                if ">" in line:
                    gi = re.sub(
                        "(>)(.*)(\|)(gi#?\|?)([0-9]*)(.*)(\n)", "\\5", line)
##                    print (gi)
                    if gi in erlist:
                        num += 1
                        new.write(line)
    print("There were " + str(num) + " sequences that couldn't be sorted due to lack of taxonomic information. \n Their SeqIDs were saved in " +
          errorfa + " if you'd care to see.")


def ShowErrorsAll(fasta, output, erlist, nadic):
    import re
    # period error
##    n = output
    n, e = output.split(".")
    num = 0
    errorfa = n + "_SubErr.txt"
    with open(fasta) as old:
        with open(errorfa, "w") as new:
            for line in old:
                if ">" in line:
                    gi = re.sub(
                        "(>)(.*)(\|)(gi#?\|?)([0-9]*)(.*)(\n)", "\\5", line)
##                    print (gi)
                    if gi in erlist:
                        num += 1
                        new.write(line)
                    for thing in nadic:
                        if gi in nadic[thing]:
                            num += 1
                            new.write(line)
    print("There were " + str(num) + " sequences that were ignored for lack of taxonomic information. \n Their SeqIDs were saved in " +
          errorfa + " if you'd care to see.")


def SSAll(fasta, number, NAnum, rank, output, SE=0, special="nothingatall"):
    MLO = MakeLists(fasta, rank)
    if SE == "yes":
        if NAnum == 0:
            ShowErrorsAll(fasta, output, MLO[2], MLO[1])
        else:
            ShowErrors(fasta, output, MLO[2])
    KE = SeqsToKeep(MLO, number, NAnum)
    return GetSeqs(fasta, KE, output, special)

# Compare
# i want to open two fasta files, and make a list of the species found in each.
# then make a new .fasta file that contains an entry for each species appears in both original fastas.
# this new entry will look like :
    # >Taxonomy|Info|Species_name|GI#|1232123|seq_id_info
    # AAALIGNEDSEQUENCEFROMGENE1AAALIGNEDSEQUENCEFROMGENE2
# this can then be made into a tree, or subsampled using SubSample.py
import re


# this makes a list of all species_names in a given .fasta in simple format.
def fastatospeclist(fasta, strain):

    allspec = []
    with open(fasta) as original:
        for line in original:
            if ">" in line:
                if "gi#|" in line:

                    spec = re.sub(
                        "([A-Za-z0-9_]*)(\|)(gi#\|)([0-9]*)(\|?)(.*)", "\\1", line)

                    lisbar = spec.split("|")

                    good = lisbar[len(lisbar) - 1]

                    do = "yes"
                elif "gi" in line:

                    spec = re.sub(
                        "([A-Za-z0-9_]*)(\|)(gi#?[\|]?)([0-9]*)(\|?)(.*)", "\\1", line)
                    lisbar = spec.split("|")

                    good = lisbar[len(lisbar) - 1]

                    spec = re.sub(
                        "(>)(.*)(\|)(.*)(gi#?[\|]?)([0-9]*)(\|?)(.*)", "\\4", line)
                    do = "yes"
                if "\n" in good:
                    good = good[:-1]
                if strain == True:
                    good = re.sub(
                        "([A-Za-z]*)(_)([A-Za-z]*)(_)(.*)", "\\1\\2\\3", good)
# print(good)
                if good in allspec:
                    pass
                else:
                    allspec.append(good)
    return allspec

# fastatospeclist("/home/abigail/Documents/Summons/Phylo/sqmo/MostRecentSQMO/August27/ERG7_5pCl_Outg_GIM.fasta")


def comparelists(list1, list2, fasta1, fasta2):
    print("started compare")
    list3 = []
    only2 = []
    only1 = []
    print("the first fasta has " + str(len(list1)) +
          " species. \nthe second fasta has " + str(len(list2)) + " species")
    for item in list1:
        if item in list2:
            list3.append(item)
        else:
            only1.append(item)
    for item in list2:
        if item in list3:
            pass
        else:
            only2.append(item)
    print(str(len(list3)) + " species found in both files.")
    print(str(len(only1)) + " species found in " + fasta1 + " only.")
    print(str(len(only2)) + " species found in " + fasta2 + " only.")
    print("In both files, there are:")
    allth = ""
    for thing in list3:
        allth = thing + " " + allth
    print (allth)
    allth = ""
    print("In " + fasta1 + ", there are:")
    for thing in only1:
        allth = thing + "\n" + allth
    print (allth)
    allth = ""
    print("In " + fasta2 + ", there are:")
    for thing in only2:
        allth = thing + "\n" + allth
    print (allth)


def CompareTwo(f1, f2, strain):
    return (comparelists(fastatospeclist(f1, strain), fastatospeclist(f2, strain), f1, f2))

# Merge(this one is old and janky, with extraneous functions not used in
# this script. sorry)


# createlists takes input files, gets all of their gi numbers, and returns
# 1. a list of all gi (once each) 2. a list of gi that are duplicates)
def CreateLists(inputlist):
    filelist = inputlist
    while "None" in filelist:
        filelist.remove("None")
    import re
    listofall = []
    listofdups = []
    for thing in filelist:
        checkfile = open(thing)
##        print("Scanning "+thing)
        for line in checkfile:
            if ">" in line:
                try:
                    # gets gi num in NCBI format
                    number = re.sub(
                        "(>)(.*)(\|)(gi#\|)([0-9]*)(\|?)(.*)", "\\5", line)

                    if ">" in number:
                        raise TypeError
                    if "\n" in number:
                        number = number[:-1]
##                        print("Yikes, "+number)
# print(numbe)
# print(number)
                    if number in listofall:
                        listofdups.append(number)
                    else:
                        listofall.append(number)
                except:
                    try:
                        # gets gi num in Simple format
                        num = re.sub(
                            "(>)(.*)(\|)(gi#?[\|]?)([0-9]*)(\|?)(.*)", "\\5", line)
# print(num)
                        if "\n" in num:
                            num = num[:-1]
##                            print("Yikes simp"+num)
                        if num in listofall:
                            listofdups.append(num)
                        else:
                            listofall.append(num)
                    except ValueError:
                        print("Error with list construction finding numbers")
                        raise SystemExit
    return(listofall, listofdups, filelist)

# list of OK subtracts listofdups from listofall


def ListOfOk(listofall, listofdups):
    listofok = []
    for num in listofall:
        if num in listofdups:
            pass
        else:
            listofok.append(num)
    return listofok


def MakeFilesKeepOneDups(listofok, listofdups, filelist, outputname, removeseconds="yes"):
    ##    print ("Got here")
    # print(removeseconds)
    try:
        o, fi = outputname.split(".")
        outputname = o
    except:
        pass
# print(removeseconds)
    info = outputname + "_GIM_Info.txt"
    dup = outputname + "_MergeExGI.fasta"
    on = outputname + "_GIM.fasta"

    write = 0
    keep = 0
    dout = 0
    listofextra = []
    listofsequences = []
    import re
    if removeseconds == "no":
        print("no removal of duplicates during merge")
        info = outputname + "_Info.txt"
        dup = outputname + "_dup_shldnt_exist.fasta"
        on = outputname + ".fasta"
        outputname = on
        with open(outputname, "w") as out:

            # opens each input file, gets the gi number, and copies it and the
            # following sequence to the corresponding output (out or dups)
            for efile in filelist:
                with open(efile) as current:
                    for line in current:
                        if ">" in line:
                            listofsequences.append(efile)
                        out.write(line)
            out.flush()
            out.close()
    else:
        outputname = on
        # opens the files to write to: output and duplicates
        with open(outputname, "w") as out:
            with open(dup, "w") as dups:
                # opens each input file, gets the gi number, and copies it and
                # the following sequence to the corresponding output (out or
                # dups)
                for efile in filelist:
                    with open(efile) as current:
                        for line in current:
                            if ">" in line:
                                listofsequences.append(efile)
                                numbe = "unassigned"
                                write = 0
                                try:
                                           # gets gi num in NCBI format
                                    numbe = re.sub(
                                        "(>)(.*)(\|)(gi#\|)([0-9]*)(\|?)(.*)", "\\5", line)
    # print(number)
                                    if ">" in number:
                                        raise TypeError
                                except:
                                                                       # gets
                                                                       # gi num
                                                                       # in
                                                                       # Simple
                                                                       # format
                                    numbe = re.sub(
                                        "(>)(.*)(\|)(gi#?[\|]?)([0-9]*)(\|?)(.*)", "\\5", line)
    # print("Simp"+num)
                                if "\n" in numbe:
                                    numbe = numbe[:-1]
    ##                                    print("Yikes simp"+num)
                                if numbe in listofok:
                                    keep += 1
                                    write = 1
                                    out.write(line)
                                    listofok.remove(numbe)
                                    listofextra.append(numbe)
                                elif numbe in listofextra:
                                    dout += 1
                                    write = 0
                                    dups.write(line)
                                else:
                                    print (
                                        "Error, GI number not found in any list")
    ##                                    print("Dups ")
    # print(listofdups)
    ##                                    print("OK ")
    # print(listofok)
                                    raise SystemExit

                            else:
                                if write == 0:
                                    dups.write(line)
                                if write == 1:
                                    out.write(line)
                out.flush()
                out.close()

    # here i am getting stats to print in the information output file.
    from collections import Counter
    counta = Counter(listofsequences).most_common()
    if removeseconds == "no":
        with open(info, "w") as information:
            information.write("You merged the following files: \n")
            for thing in counta:
                information.write("%s %s \n" % thing)
    if removeseconds == "yes":
        totseq = dout + keep
        percent = (float(dout) / float(totseq)) * 100
        perc = '%.2f' % percent
        with open(info, "w") as information:
            information.write("You merged the following files: \n")
            for thing in counta:
                information.write("%s %s \n" % thing)
            information.write("Performed Identical GI# sequence removal.\nOf " + str(totseq) + " sequences, " +
                              str(dout) + " had identical GI#s to previously copied sequences and were removed.")
            information.write("\nThat is " + perc + " percent.")
            print("Of " + str(totseq) + " sequences, " + str(dout) +
                  " had identical GI#s to previously copied sequences and have been removed. \n That is " + perc + " percent.")
    print("You merged file:number of sequences :")
    for thing in counta:
        print("%s %s" % thing)
    print("\n Output file is " + outputname + ". \n Info at " + info + ".")
    return outputname


# this will merge all files, and keep only one sequence for each distinct
# GI number.
def MergeKeepOneGI(outputname, inputlist):
    print("Beginning Merge and Remove Duplicates (Keeping One Copy)")
    print("GI number accession initiated")
    listofall, listofdups, filelist = CreateLists(inputlist)
    if listofdups == []:
        print("You have no duplicates")
        listofok = listofall
        merge = MakeFilesKeepOneDups(
            listofok, listofdups, filelist, outputname)
    else:
        listofok = listofall
        print("All lists created. Sorting sequences...")
# print(listofok)
        merge = MakeFilesKeepOneDups(
            listofok, listofdups, filelist, outputname, "no")
    return merge

# this merges all files, preforming no removal or conflict resolution


def MergeNoRemove(outputname, inputlist):
    print("Beginning Merge (No Removal)")
    listofall, listofdups, filelist = CreateLists(inputlist)
    listofok = listofall
    emptylist = []
    merge = MakeFilesKeepOneDups(
        listofok, emptylist, filelist, outputname, "no")
    print("Merge Complete")
    return merge


# finds the GI number of all sequences. Keeps only sequences that have unique GI numbers, discarding any that have twins etc
# because they are likely two partial alignments from the same gene, which
# skews results.

def RemoveDupGIs(filein):
    import re
    print("Removing all GI-sharing sequences from " + filein)
    try:
        fin, ext = filein.split(".")
        info = fin + "_Info.txt"
        writeto = fin + "_NDGIs.fasta"
        writedups = fin + "_ExcludedDupGIs.fasta"
    except:
        info = filein + "_Info.txt"
        writeto = filein + "_NDGIs.fasta"
        writedups = filein + "_ExcludedDupGIs.fasta"
    # initializing lists to remove dulplicates
    listofnums = []
    listofdups = []
    listofok = []
    kept = 0
    removed = 0
    dup = "no"
    # finding the gi numbers of all sequences with at least one duplicate
    with open(filein) as original:
        for line in original:
            if ">" in line:
                number = re.sub(
                    "(>)(gi)(\|)([0-9]*)(\|)([A-Za-z]*)(\|)(.*)", "\\4", line)
                num = number[:-1]
                if num in listofnums:
                    ##                    print("Found a duplicate! "+num)
                    listofdups.append(num)
                else:
                    listofnums.append(num)
    # if no duplicates, setting list of all to be ok
    if listofdups == []:
        listofok = listofnums
    # otherwise making a list of all ok gi numbers.
    else:
        ##        print("Duplicate sequences will be removed and saved in "+writedups)
        for gi in listofnums:
            if gi in listofdups:
                pass
            else:
                listofok.append(gi)
    # dup is on/off switch to see if sequence should be written to good or dup
    # file
    dup = "no"
    orig = open(filein)
    # opening a file to save duplicate, and file to save everything but
    # duplicates
    with open(writedups, "w") as duplicates:
        with open(writeto, "w") as new:
            for line in orig:
                if ">" in line:
                    # fins gi number in GI format
                    try:
                        number = re.sub(
                            "(>)(gi)(\|)([0-9]*)(\|)([A-Za-z]*)(\|)(.*)", "\\4", line)
                        num = number[:-1]
                    except:
                        junk, numb = line.split("gi#")
                        num = numb[:-1]
                    if num in listofok:
                        new.write(line)
                        dup = "no"
                        kept += 1
                    elif num in listofdups:
                        dup = "yes"
                        duplicates.write(line)
                        removed += 1
                    else:
                        print("num")
                        print("gi is nowhere error")
                        raise SystemExit
                else:
                    if dup == "yes":
                        duplicates.write(line)
                    elif dup == "no":
                        new.write(line)
                    else:
                        print("Error in writing sequence")
                        raise SystemExit
    import os
    # this ensures that the function passes on a full directory-level name of the file in question
    # probably not necessary inside everything() functions that already go to the dir in question,
    # but it doesn't hurt to keep it.
    # currently OOO b/c was annoying. change the return key to change it.
    filelocation = os.getcwd() + "/" + writeto
    comb = kept + removed
    with open(info, "a") as inform:
        inform.write("\n\nPerformed Identical GI# Sequence Removal.\nOf " + str(comb) + " sequences, " + str(removed) +
                     " were removed for sharing GI#s.\nNew file at " + writeto + ". \nRemoved sequences at " + writedups)
    print("Of " + str(comb) + " sequences, " + str(removed) +
          " had shared GI numbers and were removed.")
    print("New file at " + writeto + ".")
    # returns the name of the newly created no duplicates file
    return writeto


def RemoveDupSeqIDs(filein):
    import re
    print("Removing all but one sequence for each SeqID from: " + filein)
    try:
        fin, ext = filein.split(".")
        info = fin + "_Info.txt"
        writeto = fin + "_NDSID.fasta"
        writedups = fin + "_ExcludedDupSID.fasta"
    except:
        info = filein + "_Info.txt"
        writeto = filein + "_NDSID.fasta"
        writedups = filein + "_ExcludedDupSID.fasta"
    # initializing lists to remove dulplicates
    listoftaxa = []
    kept = 0
    removed = 0
    orig = open(filein)
    dup = "yes"
    # opening a file to save duplicate, and file to save everything but
    # duplicates
    with open(writedups, "w") as duplicates:
        with open(writeto, "w") as new:
            for line in orig:
                if ">" in line:

                    if line in listoftaxa:
                        duplicates.write(line)
                        dup = "yes"
                        removed += 1
                    else:
                        dup = "no"
                        new.write(line)
                        kept += 1
                        listoftaxa.append(line)
                else:
                    if dup == "yes":
                        duplicates.write(line)
                    elif dup == "no":
                        new.write(line)
                    else:
                        print("Error in writing sequence")
                        raise SystemExit
    import os
    filelocation = os.getcwd() + "/" + writeto
    comb = kept + removed
    with open(info, "a") as inform:
        inform.write("\n\nPerformed identical Squence ID sequence Removal.\nOf " + str(comb) + " sequences, " + str(removed) +
                     " were removed for having the same Species as a previous sequence.\nNew file at " + writeto + ". \nRemoved sequences at " + writedups)
    print("Of " + str(comb) + " sequences, " + str(removed) +
          " had identical Species_name to previous sequences and were removed.")
    print("New file at " + writeto + ".")
    # returns the name of the newly created no duplicates file
    return writeto

# this removes extraneous sequences by Taxa_name. It keeps the first sequence it finds of each unique Taxa_name,
# and removes any subsequent sequence with the same name.
# This just cuts down on "clutter" in the final tree


def RemoveDupTaxa(filein):
    import re
    print("Removing all but one sequence for each Species_name from: " + filein)
    try:
        fin, ext = filein.split(".")
        info = fin + "_Info.txt"
        writeto = fin + "_NDTaxa.fasta"
        writedups = fin + "_ExcludedDupTaxa.fasta"
    except:
        info = filein + "_Info.txt"
        writeto = filein + "_NDTaxa.fasta"
        writedups = filein + "_ExcludedDupTaxa.fasta"
    # initializing lists to remove dulplicates
    listoftaxa = []
    kept = 0
    removed = 0
    orig = open(filein)
    dup = "yes"
    # opening a file to save duplicate, and file to save everything but
    # duplicates
    with open(writedups, "w") as duplicates:
        with open(writeto, "w") as new:
            for line in orig:
                if ">" in line:
                    # finds taxa number in Shorten format
                    if "gi#" in line:
                        tax, junk = line.split("|gi#")
                        taxa = tax[1:]
                    # finds taxa number in NCBI .fasta format
                    else:
                        t = re.sub(
                            "(>)(gi)(\|)([0-9]*)(\|)([A-Za-z]*)(\|)(.*)(\[\'?[A-Z]?[a-z]* ?.*\])(.*)", "\\9", line)
                        ta = re.sub("[\[\]]", "", t)
                        tax = re.sub(" ", "_", ta)
                        taxa = tax[:-1]
                    if taxa in listoftaxa:
                        duplicates.write(line)
                        dup = "yes"
                        removed += 1
                    else:
                        dup = "no"
                        new.write(line)
                        kept += 1
                        listoftaxa.append(taxa)
                else:
                    if dup == "yes":
                        duplicates.write(line)
                    elif dup == "no":
                        new.write(line)
                    else:
                        print("Error in writing sequence")
                        raise SystemExit
    import os
    filelocation = os.getcwd() + "/" + writeto
    comb = kept + removed
    with open(info, "a") as inform:
        inform.write("\n\nPerformed identical Taxa Name Sequence Removal.\nOf " + str(comb) + " sequences, " + str(removed) +
                     " were removed for having the same Species as a previous sequence.\nNew file at " + writeto + ". \nRemoved sequences at " + writedups)
    print("Of " + str(comb) + " sequences, " + str(removed) +
          " had identical Species_name to previous sequences and were removed.")
    print("New file at " + writeto + ".")
    # returns the name of the newly created no duplicates file
    return writeto


# this removes extraneous sequences by AA sequence.
# It keeps the first sequence it finds of each unique AA sequence, and removes any subsequent matching seqs.
# This just cuts down on "clutter" in the final tree

def RemoveDupSeqs(filein):
    import re
    print("Removing extraneous identical AA Sequences from " + filein)
    try:
        fin, ext = filein.split(".")
        writeto = fin + "_NDAAs.fasta"
        writedups = fin + "_ExcludedDupAAs.fasta"
        info = fin + "_Info.txt"
    except:
        info = filein + "_Info.txt"
        writeto = filein + "_NDAAs.fasta"
        writedups = filein + "_ExcludedDupAAs.fasta"
    # initializing lists to remove dulplicates
    listofseqs = []
    removed = 0
    kept = 0
    a = 0
    # finding the gi numbers of all sequences with at least one duplicate
    with open(filein) as original:
        with open(writedups, "w") as duplicates:
            with open(writeto, "w") as new:
                for line in original:
                    if ">" in line:
                        if a == 0:
                            prevline = line
                        else:
                            if allline in listofseqs:
                                duplicates.write(prevline)
                                duplicates.write(allline)
                                removed += 1
                                allline = ""
                            else:
                                listofseqs.append(allline)
                                new.write(prevline)
                                new.write(allline)
                                kept += 1
                                allline = ""
                        prevline = line
                    else:
                        if a == 0:
                            a = 1
                            allline = line
                        else:
                            mline = allline + line
                            allline = mline
                if allline in listofseqs:
                    duplicates.write(prevline)
                    duplicates.write(allline)
                    removed += 1
                else:
                    listofseqs.append(allline)
                    new.write(prevline)
                    new.write(allline)
                    kept += 1
    import os
    # this ensures that the function passes on a full directory-level name of the file in question
    # probably not necessary inside everything() functions that already go to the dir in question,
    # but it doesn't hurt to keep it.
    # currently OOO b/c was annoying. change the return key to change it.
    filelocation = os.getcwd() + "/" + writeto
    comb = kept + removed
    with open(info, "a") as inform:
        inform.write("\n\nPerformed Identical AA Sequence Removal.\nOf " + str(comb) + " sequences, " + str(removed) +
                     " were removed for being identical to a previous sequence.\nNew file at " + writeto + ". \nRemoved sequences at " + writedups)
    print("Of " + str(comb) + " sequences, " + str(removed) +
          " were removed for being identical to a previous sequence.")
    print("New file at " + writeto + ".")
    return writeto


def informationallscript():
    import textwrap
    print textwrap.dedent("""\



FEAST (fasta editing and subsampling tool) Info

Welcome to my script. It is useful for modifying and subsampling .fasta files based on seqid and taxonomic information.

___________________________
Commands:
Simplify/SimplifyKeep
Remove
PFAM
AppendRanks
Seperate
Extract
SubSample
Summarize
Compare
MultiDataSubSample
Merge
___________________________
Bugs:
please only use files named foo.fasta as input. not ending the string in .fasta may cause an error. having extra periods in the string may cause an error.
test.fasta = ok
test.example.fasta = may error
test.txt = may error
___________________________
Examples:
	FEAST.py /path/to/directory example.fasta -si -ar \"kingdom class\" -ex \"Eukaryotes\" -ss \"2 1\"
This will look at example.fasta in directory. It will \"simplify\" the SeqIDs from NCBI online download or BLASTP outfmt \"6 sseqid stitle sseq\". Then it will append kingdom and class information to each seqID (>oldseqID ---> >kingdom|class|oldseqID) via \"appendrank\". It will then \"extract\" all sequences with the string \"Eukaryotes\". Finally, it will perform subsampling on the extracted eukarotes-containing sequences. The output will contain one random sequence per unique depth-2 string (in this case depth-2 corresponds to \"class\").
	Generated files will include:
		exampleSh.fasta (simplified)
		exampleShKC.fasta (with kingdom and class appended)
		exampleShKCExtracted.fasta (only containing seqs with str \"Eukaryotes\"
		exampleShKCExtractedSubSamp1p2.fasta (subsampled 1 random sequence per class)

	FEAST.py /path/to/directory example.fasta -sk -rm TX -pf \"L, Test_N_Term Test_C_Term\"
This will look at example.fasta in directory. It will simplify the SeqIDs from NCBI online download or BLASTP outfmt \"6 sseqid stitle sseq\", keeping gene information via \"simplifykeep\". Then it will keep only one sequence per unique species using the \"remove TX\" command. Each resulting sequence will then be vetted for the presence of conserved domains Test_N_Term and Test_C_Term via a local install of Pfam - sequences with both domains will be kept.
	Generated files will include:
		exampleMo.fasta (simplified keeping info)
		exampleMoNDTaxa (removed duplicates by Taxon)
		exampleMoNDTaxaVet (after Pfam vetting for Test_N_Term and Test_C_Term)
___________________________
Command Details
___________________________
Directory:
Usage:
	FEAST.py /path/to/directory
Purpose:
	Specifies a directory to run in - all input files should be here, and all output files should write to here.

Fasta:
Usage:
	FEAST.py /path/to/directory example.fasta
	FEAST.py /path/to/directory \"example1.fasta example2.fasta example3.fasta\"
Purpose:
	Specifies a fasta file or files to run on. If multiple files are listed, all operations will be done to each included file. Note that Compare requires exactly 2 files to be listed, and MultiDataSubsampling requires a minimium of 2.
___________________________
Simplify
Usage:
	-simplify
	-si
Example:
	FEAST.py /path/to/directory example.fasta -si
		This will produce exampleSh.fasta, with all of the sequences of example.fasta but with shortened names.
Purpose:
	Takes output from NCBI Blast webserver download or otherwise formatted like >gi|587572769|gb|AHK06412.1|:1-224 manganese superoxide dismutase [Megalobrama amblycephala]
	and replaces the seqID with just >Species_name|gi#|##### (in this case, >Megalobrama_amblycephala|gi|587572769).
	Removes disallowed characters like spaces, :;,` etc.
	This operation (or simplifykeep) is necessary for to be preformed on a .fasta before AppendRank, SubSample, Compare, Summarize, or MDS will work - as they rely on this format.
How it works:
	This is just regex to get a consistant seqID placement of species and gi#.
	NOTE that sequences that do not fit the NCBI-out input pattern (eg; lack a given species name) will be saved to a TaxaError file.
___________________________
SimplifyKeep
Usage:
	-simplifykeep
	-sk
Example:
	FEAST.py /path/to/directory example.fasta -sk
		This will produce exampleMo.fasta, with all of the sequences of example.fasta but with shortened names.
Purpose:
	Takes output from NCBI Blast webserver download or otherwise formatted like >gi|587572769|gb|AHK06412.1|:1-224 manganese superoxide dismutase [Megalobrama amblycephala]
	and replaces the seqID with just >Species_name|gi#|#####|Gene_Info (in this case, >Megalobrama_amblycephala|gi|587572769|1_224_manganese_superoxide_dismutase).
	Removes disallowed characters like spaces, :;,` etc.
	This operation (or simplify) is necessary for to be preformed on a .fasta before AppendRank, SubSample, Compare, Summarize, or MDS will work - as they rely on this format.
How it works:
	This is just regex to get a consistant seqID placement of species and gi#.
	NOTE that sequences that do not fit the NCBI-out input pattern (eg; lack a given species name) will be saved to a *TaxaError file.
___________________________
Remove
Usage:
	-remove mode
	-rm mode
Example:
	FEAST.py /path/to/directory example.fasta -re SI
		This will produce a new file with duplicates removed based on seqID (keeping one per seqID).
	FEAST.py /path/to/directory example.fasta -re \"SI TX\"
		This will produce a new file with duplicates removed based on seqID (keeping one per seqID), and also based on Taxon_name (keeping one per taxon).
Variables:
	mode - 	chose from one of four modes
		AA - keep one per identical amino acid sequence
		TX - keep one per identical Taxon_name
		SI - keep one per identical SequenceID
		GI - if a GI number is present in more than one sequence, remove them all.
Purpose:
	Removes duplicates in various ways. Reduces number of sequences.
How it works:
	Makes a new file without duplicates based on seqID info or sequence.

___________________________
PFAM
Usage:
	-pfam \"mode, conserved_domain(s)\"
	-pf \"mode, conserved_domain(s)\"
Example:
	FEAST.py /path/to/directory example.fasta -pf \"O, Example_C_term\"
		This will run online pfam, keeping only sequences identified to have an Example_C_term domain.
	FEAST.py /path/to/directory example.fasta -pf \"L, Example_C_term Example_N_Term\"
		This will run local pfam, keeping only sequences identified to have an Example_C_term domain and an Example_N_term domain.
	FEAST.py /path/to/directory example.fasta -pf \"F<examplePfamFile.txt>, Example_C_term\"
		This will take a previously-generated pfam file that corresponds to your given .fasta, and parse both, creating a new fasta including only sequences identified to have an Example_C_term domain.
Variables:
	mode - 	Chose from one of three modes to run in:
		Online - Type an \"O\". This requires you to manually upload a generated file to pfam online database (http://pfam.xfam.org/), wait until it is done running, and then save the output in the directory folder as indicated. Chose this if you want to run pfam, but do not have it installed locally.
		Local - Type an \"L\". Ensure that you have pfam_scan.pl installed, and have the pfamdatabase in your path, or have edited the third line of this script containing the variable pfamdatabasedir to point to your Pfam data files eg \"/home/abigail/Downloads/PfamScan/PfamData/\".
		File - Type an \"F\" followed by the name of your Pfam file within <>. The specified Pfam file will be parsed, instead of re-running pfam at all.
	conserved_domain(s) - Type one or mode space-seperated pfam domains to vet for the presence of. If you use multiples, vetting will require the presence of ALL indicated domains.
Purpose:
	Allows for integrated vetting of sequences based on their conserved domains. This may ensure that any given sequence after vetting is actually a representative of your gene(s) of interest. Type \"mode,conserved_domain\" including quotes and the comma.
How it works:
	Pfam gets the conserved domain information for each sequence. This script parses Pfam output, correlates it with your .fasta file, and makes a new file containing only sequences that were found by Pfam to have your indicated domain(s).
___________________________
AppendRanks
Usage:
	-appendranks rank
	-ar rank
Example:
	FEAST.py /path/to/directory example.fasta -ar \"kingdom\"
		This will produce a new file called exampleK.fasta that contains each sequence from example with seqIDs formatted like >kingdom|old_seqID
	FEAST.py /path/to/directory example.fasta -re \"kingdom phylum order\"
		This will produce a new file called exampleK.fasta that contains each sequence from example with seqIDs formatted like >kingdom|phylum|order|old_seqID
Variables:
	rank - include within quotes a space-seperated list of taxonomic ranks to append to seqID. All data is pulled from NCBI, so you can include things like common_name or subclass etc.
Purpose:
	Adds taxonomic information for ease of understanding produced trees, or for use by other parts of the script such as subsampling.
How it works:
	Parses each sequence for it's GI number, queries NCBI for the associated TAXID, and then queries NCBI taxonomy database for the associated taxonomic rank(s) specified. If there is no given taxonomic information, the script will print \"NA\" in the specified place. If no taxonomic information is available at all (wrong formatting of SeqIDs, sequence not from NCBI, or no internet access) the script will print NT| and move on.

___________________________

Seperate
Usage:
	-seperate string
	-sp string
Example:
	FEAST.py /path/to/directory example.fasta -sp Bacteria
		this will split example.fasta into two files, exampleChosen.fasta (containing all sequences that have \"Bacteria\" in the seqID) and exampleNotChosen.fasta (containing all the rest)
	FEAST.py /path/to/directory example.fasta -sp \"Cyanobacteria Proteobacteria\"
		Same as above, but any sequence containing \"Cyanobacteria\" or \"Proteobacteria\" will be in the exampleChosen.fasta file, the rest will be in exampleNotChosen.fasta.
	FEAST.py /path/to/directory example.fasta -sp \"Cyanobacteria Proteobacteria keepone\"
		The exampleChosen.fasta file will have the first sequence whose ID contains \"Cyanobacteria\" and the first sequence whose ID contains \"Proteobacteria\". The exampleNotChosen file will have all other sequences.
Variables:
	string - this should be the thing you want to check seqIDs for. if using multiples, keep them spece-seperated and within quotes.
Purpose:
	Sometimes you want to seperate your data, this makes it quick and easy to do.
	If any of the input string(s) are found in the seqID of a given sequence, it will be written to a *Chosen file.
	If none of the strings are found, it will be written to a *NotChosen file.
	Including \"keepone\" in the given list will change behaviour to keep just the first sequence containing a given string in the *Chosen file, and send the rest to NotChosen.
How it works:
	Literally just writes to different files based on seqID names.
___________________________
Extract
Usage:
	-extract string
	-ex string
Example:
	FEAST.py /path/to/directory example.fasta -ex Bacteria
		This will create a new file called exampleExtracted.fasta that contains only sequences from example.fasta that had the string \"Bacteria\" in their SeqID.
	FEAST.py /path/to/directory example.fasta -ex \"Cyanobacteria Proteobacteria\"
		Same as above, but any sequence containing \"Cyanobacteria\" or \"Proteobacteria\" will be in the exampleExtracted.fasta file.
	FEAST.py /path/to/directory example.fasta -ex \"Cyanobacteria Proteobacteria keepone\"
		The exampleExtracted.fasta will have the first sequence whose ID contains \"Cyanobacteria\" and the first sequence whose ID contains \"Proteobacteria\".
Variables:
	string - this should be the thing you want to check seqIDs for. if using multiples, keep them space-seperated and within quotes.
Purpose:
	Sometimes you want just a subset of your data, this makes it quick and easy to do.
	If any of the input string(s) are found in the seqID of a given sequence, it will be written to a *Extracted file.
	Including \"keepone\" in the given list will change behaviour to keep just the first sequence containing a given string in the *Extracted file.
How it works:
	Writes to a new file if seqeuenceID includes given terms.

___________________________
SubSample
Usage:
	-subsample \"ranknum, seqnum nadrop(optional)\"
	-ss \"ranknum, seqnum nadrop(optional)\"
Example:
	FEAST.py /path/to/directory example.fasta -ss \"1,1\"
		This will produce a new file called example1per1.fasta that contains a subsampling of sequences from example.fasta. If seqIDs from example look like: >kingdom|phylum|order|old_seqID, then this will write one random sequence per unique kingdom-string to the new file.
	FEAST.py /path/to/directory example.fasta -ss \"1,1 nadrop\"
		This will produce a new file called example1per1.fasta that contains a subsampling of sequences from example.fasta. If seqIDs from example look like: >kingdom|phylum|order|old_seqID, then this will write one random sequence per unique kingdom-string to the new file. If rank = \"NA\", look one rank lower.
	FEAST.py /path/to/directory example.fasta -ss \"2,3\"
	This will produce a new file called example3per2.fasta that contains a subsampling of sequences from example.fasta. If seqIDs from example look like: >kingdom|phylum|order|old_seqID, then this will write three random sequences per unique order-string to the new file.
Variables:
	ranknum - how many bars deep in the appended taxonomy to look. for example, the seqID >kingdom|phylum|order|old_seqID has \"kingdom\" one deep, \"phylum\" two deep, \"order\" three deep, \"old_seqID\" four deep. the script will compare all strings of a specified depth to each other, such that you can do subsampling at whichever taxonomic rank you require. Ensure that all sequences had taxonomy appended in the same way, else your subsampling will not be representative.
	seqnum - the number of sequences per unique string to keep.
	nadrop - type \"nadrop\" to toggle on. if NCBI taxonomy returned \"NA\" for a given sequenceID-rank combo, look one rank lower for that sequence only. this can preserve data that might otherwise have been lost, but could result in oversampling (for example, all cyanobacteria are missing class on NCBI) - use with caution.
Purpose:
	Subsamples evenly across a dataset at any given taxonomic depth.
How it works:
	Parses each sequenceID for the string at specified depth(ranknum), starting at one and incrementing at bars. All strings are compared, and matching-string sequences are pooled together. A given number(seqnum) are randomly selected from each pool and printed to the new document.
  ___________________________
Summarize
Usage:
	-summarize ranknum
	-su ranknum
Example:
	FEAST.py /path/to/directory example.fasta -su 1
		This will print a list of all strings in rank 1 (kingdom if the seqIDs are like >kingdom|phylum|order|old_seqID), and how many taxa have that unique string. eg
		Bacteria 12
		Eukaryota 3
Variables:
	ranknum - how many bars deep in the appended taxonomy to look. for example, the seqID >kingdom|phylum|order|old_seqID has \"kingdom\" one deep, \"phylum\" two deep, \"order\" three deep, \"old_seqID\" four deep. the script will compare all strings of a specified depth to each other, such that you get an accurate summary at whichever taxonomic rank you require. Ensure that all sequences had taxonomy appended in the same way, else your summary will not be representative.
Purpose:
	Quickly summarize large datasets at whatever precision you require.
How it works:
	Parses each sequenceID for the string at specified depth(ranknum), starting at one and incrementing at bars. All strings are compared, and each time a duplicate comes up that string tally is incremented by one. Results are printed highest tally -> lowest
___________________________
Compare
Usage:
	-compare
	-co
Example:
	FEAST.py /path/to/directory \"example.fasta example2.fasta\" -co
		This will compare the species_names in exactly two fasta files, and print 1. all taxa in both example.fasta and example2.fasta, 2. all taxa only in example.fasta, and 3. all taxa only in example2.fasta
Purpose:
	Quickly find what taxa are shared in two files, and which are not.
How it works:
	Uses regex to identify the species name (given seqIDs in shortened / shorten-appended format), and checks between two files.
___________________________
MultiDataSubSample
Usage:
	-multidatasub \"ranknum seqnum strain(optional) nadrop(optional)\"
	-ms \"ranknum seqnum strain(optional) nadrop(optional)\"
Example:
       FEAST.py /path/to/directory \"example.fasta example2.fasta\" -ms \"4 1\"
		This will take one sequence per unique string of depth(ranknum) 4, biased towards taking the same taxon from both example and example2 datasets.
Variables:
	ranknum - how many bars deep in the appended taxonomy to look. for example, the seqID >kingdom|phylum|order|old_seqID has \"kingdom\" one deep, \"phylum\" two deep, \"order\" three deep, \"old_seqID\" four deep. the script will compare all strings of a specified depth to each other, such that you can do subsampling at whichever taxonomic rank you require. Ensure that all sequences had taxonomy appended in the same way, else your subsampling will not be representative.
	seqnum - the number of sequences per unique string to keep.
	nadrop - type \"nadrop\" to toggle on. if NCBI taxonomy returned \"NA\" for a given sequenceID-rank combo, look one rank lower for that sequence only. this can preserve data that might otherwise have been lost, but could result in oversampling (for example, all cyanobacteria are missing class on NCBI) - use with caution.
	strain - type \"strain\" to toggle on. this is an option when working with data that has species name written like Species_name_sample1 or Species_name_strain2x5 etc and you want to ignore everything after Species_name for purposes of determining presence of species across datasets.
Purpose:
   Create subsampled set of sequences across several datasets, preserving taxonomic diversity yet optimizing for shared taxa for purposes of comparison and/or concatenation.
How it works:
  This script parses your input .fasta files. It requires that some taxonomic rank be appended to each sequence in the format outputted by AppendRank. Looking at the level of taxonomic rank specified by you, it looks to see what datasets include sequences sharing each unique rank (eg \"which datasets have at least one sequence with order = \"rodent\"). It will then look at the specific species within that rank across all datasets, and return the single species found in the largest number of datasets (eg if mouse found in 3/4 datasets and rat in 2/4, it will return mouse as the chosen sequence). If there are several equally well represented options, it randomly choses one to return. For any datasets that have a sequence for the chosen species, that sequence will be written to that dataset's output file. For any dataset with at least one sequence of the specified rank, but none of the chosen species, a different species will be randomly chosen as a substitute to preserve taxonomic breadth (eg a dataset with rodent sequences but no mouse will have a random other rodent sequence written to file). This process will continue until the script has examined all strings of the rank specified in all input files, and chosen a representative species for each (and alternates if necessary). This ensures that in all datasets, sampling properly captures the breadth of the taxa represented.
  An information file will be created, including what ranks and species are shared across datasets, what representative species were chosen, and what (if any) substitutions were made.
___________________________
Merge
Usage:
	-merge
	-me
Example:
	FEAST.py /path/to/directory \"example.fasta example2.fasta\" -me
		This will merge all sequences from both example.fasta and example2.fasta in a new merged file.
Purpose:
	Merge many files into one large file - more efficient than manual merging given a large number of .fasta files.
How it works:
	Just copies sequences over into one new merged file. Very simple.
___________________________
Information
Usage:
	-i
	-information
Purpose:
	prints this file.


""")


import textwrap
examp = textwrap.dedent("""\
Examples
\tFEAST.py /path/to/directory example.fasta -si -ar \"kingdom class\" -ex \"Eukaryotes\" -ss \"2 1\"\nThis will look at example.fasta in directory. It will \"simplify\" the SeqIDs from NCBI online download or BLASTP outfmt \"6 sseqid stitle sseq\". Then it will append kingdom and class information to each seqID (>oldseqID ---> >kingdom|class|oldseqID) via \"appendrank\". It will then \"extract\" all sequences with the string \"Eukaryotes\". Finally, it will perform subsampling on the extracted eukarotes-containing sequences.The output will contain one random sequence per unique depth-2 string (in this case depth-2 corresponds to \"class\").\n\n\t FEAST.py /path/to/directory example.fasta -sk -rm TX -pf \"L, Test_N_Term Test_C_Term\"\nThis will look at example.fasta in directory.It will simplify the SeqIDs from NCBI online download or BLASTP outfmt \"6 sseqid stitle sseq\", keeping gene information via \"simplifykeep\". Then it will keep only one sequence per unique species using the \"remove TX\" command.Each resulting sequence will then be vetted for the presence of conserved domains Test_N_Term and Test_C_Term via a local install of Pfam - sequences with both domains will be kept.")
    """)


#
#___________________parser_________

if __name__ == "__main__":

    print("Running in terminal")
    print ("\nThis is a Fasta Editing, Annotation, and Subsampling tool (FEAST) used to modify the contents of fasta files.\nFor extended information type '-i', for simple help: '-h'.")
    print("\nSimplify or SimplifyKeep must be preformed on a .fasta before AppendRank, SubSample, Compare, Summarize, or MDS will work. AppendRank must be preformed for SubSample, Summarize, and MDS to work. This is not an error, just a reminder :) \n\n")

    import sys
    import argparse
    import os
    import re
    from argparse import ArgumentParser
    parser = argparse.ArgumentParser(
        description="All Options", epilog=examp, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("directory", nargs='?', default=os.getcwd(
    ), type=str, help="type name of directory to run in (where .fasta(s) reside)")
    parser.add_argument("FASTA", nargs='?', default="no_fasta_specified", type=str,
                        help="type the name of your .fasta file or a space seperated list of files within quotes")
    parser.add_argument("-si", "--simplify", action="store_true",
                        help="This simplifies fasta file names to >Species_name|gi#|#####")
    parser.add_argument("-sk", "--simplifykeep", action="store_true",
                        help="This simplifies fasta file names to >Species_name|gi#|#####|OtherInfo")
    parser.add_argument("-re", "--remove", action="store", dest="RE",
                        help="Type abbrev(s) you want: AA = keep one per unique AA seq. TX = keep one per species. SI = remove all but one of each SeqID. GI = remove all that share a GI number.")
    parser.add_argument("-pf", "--pfam", action="store", dest="PF",
                        help="Type in quotes O online L local F <file> preexisting file, comma, space-sep domain list eg 'L, C_term N_term' or 'F gene.txt, Domain1'")
    parser.add_argument("-ar", "--appendranks", action="store", dest="AR",
                        help="Type ranks to append eg \"superkingdom phylum order\"")
    parser.add_argument("-sp", "--seperate", action="store", dest="SP",
                        help="Type string (or space-sep list of strs in quotes) to query seqIDs. Creates 2 new fasta, chosen and not-chosen")
    parser.add_argument("-ex", "--extract", action="store", dest="EX",
                        help="Type string (or space-sep list of strs in quotes) to search seqIDs for. Creates new fasta of seqs w str. default: keep all matching. can also keep 1per str by typing in list 'keepone'")
    parser.add_argument("-ss", "--subsample", action="store", dest="SS",
                        help="Type rank # to look at, number to keep. eg \"2 1\" keeps 1 per phylum if seqids = >sk|phy|ord|tax|gi|##. Type '# # nadrop'  to toggle keep one downrank if NA")
    parser.add_argument("-su", "--summarize", action="store", dest="SU",
                        help="Type rank # to look at while summarizing number of sequences in each category")
    parser.add_argument("-co", "--compare", action="store_true",
                        help="Toggle compares two (exactly 2) .fasta files and describes what taxa are shared/in each")
    parser.add_argument("-ms", "--multidatasetsubsample", action="store", dest="MS",
                        help="Type rank, number -eg \"2,1\". toggles creation of list of taxa number per rank, optimizing for highest representation across datasets")
    parser.add_argument("-me", "--merge", action="store", dest="ME",
                        help="This merges all listed fasta files into one new file. Type name for new file")
    parser.add_argument("-i", "--information",
                        action="store_true", help="Prints extended info")


args = parser.parse_args()


if args.information:

    informationallscript()


directory = args.directory
os.chdir(directory)

fastalist = args.FASTA

try:
    a = fastalist.split()
except:
    a = [fastalist]

#
# how to keep threading? working = most recent result. ensure each fxn returns the newly-generated file name
    # do this
    # write cute print-outs to track where the script is
b = []
for fasta in a:
    if "./" in a:
        a = a[:-2]
    working = fasta
    print("Opening file: " + working)
    if args.simplifykeep:
        print("Beginning shorten, keeping information")
        working = ShortenKeep(working)
        print("Done")
    if args.simplify:
        print("Beginning shorten")
        working = ShortenNoKeep(working)
        print("Done")
    if str(args.RE) == "None":
        pass
    else:
        print("Beginning remove duplicates...")
        working = Remove(working, args.RE)
    if str(args.PF) == "None":
        pass
    else:
        print("Beginning Pfam Vetting... This may take some time.")
        working = PfamSelect(args.PF, working, directory)
        print("Done")
    if str(args.AR) == "None":
        pass
    else:
        print("Adding Rank Information... This may take some time.")
        rlist = args.AR.split()
        rlist.reverse()
        working = AddRank(working, rlist, directory)
        print("Done")
    if str(args.SP) == "None":
        pass
    else:
        print("Beginning seperation...")
        working = Seperate(directory, working, args.SP)
        print("Done")
    if str(args.EX) == "None":
        pass
    else:
        print("Beginning extraction...")
        working = Extract(directory, working, args.EX)
        print("Done")
    if str(args.SS) == "None":
        pass
    else:
        print("Beginning SubSampling...")
        working = SubSampling(working, args.SS)
        print("Done")
    if str(args.SU) == "None":
        pass
    else:
        print("Running Summarize on file: " + working)
        SumTaxa(working, args.SU)
        print("Done")
    b.append(working)
if args.compare:
    if len(b) == 2:
        one = b[0]
        two = b[1]
        print("Beginning Compare...")
        whatever = CompareTwo(one, two, "no")
    else:
        print("error can not compare unless there are exactly two fasta files given")
if str(args.MS) == "None":
    pass
else:
    print("Beginning MultiDataset SubSampling Selection")
    working = MultiDataSub(args.MS, b, args.directory)
    print("Done")
if str(args.ME) == "None":
    pass
else:
    print("Beginning file merge")
    working = MergeMerge(b, args.ME)
    print("Done")

print("Finished! Exiting.")
