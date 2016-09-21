#!/usr/bin/python

#version 4 of tree parser, with support for running and comparing each subtree and picking two most divergent non-transfer sequences.
#sept 20 2016 abigailc@actaeon

#added support for printing fast_subtree, fasta of sttips, list of species, superkingdom
#added send-to-make-species -tree-support
#added make-raxml support
#added ranger compare support
#added SS after ranger support.



#jun 9 added r script to get subtrees
#jun 17 fixed repeated last tree error, added min value
#july - added subsampling support
#september: adding support for generation of RAXML and SPECIES trees.

#current plan: 1. generate subtrees list using R script
#2. check output for monophyly of a given string
#3. next check for monophyly with inclusions or with outside artifacts
#by scanning for "best clade" which is clade w highest numberstr*percentstr
#and then counting the number of str outside of best clade, and comparing to a para/artifact cut off
#if best clade does not have 100% strtipspercent, mark it inclusion
#if best clade does not include each strtip, mark it artifacts if under cutoff, else poly

def MakeSubtreesFile(nexus, subtrees, verbose = False):
    print("Making subtreesfile")
    import os
    isfile = os.path.isfile("./"+subtrees)
    if isfile == True:
        print("There is already a subtrees file made from this subtree! \n it is called: "+subtrees+".")
        print("You can indicate it with the -t flag to save time")
        answer = raw_input("Type \"continue\" if you would like to proceed anyways.  :")
        print(answer)
        if answer == "continue":
              pass
        else:
              print("Quitting")
              raise SystemExit
    os.system("subtreegen "+nexus+" "+subtrees)
    return subtrees

##class Subtree(object):
##    def __init__(self, name, tips):
##        self.name = name
##        self.tips = tips
##    def GetName(self):
##        return self.name
##    def GetTips(self):
##        return self.tips
def GetRankFromSubtreesFile(subtrees, ranknum, verbose = False):
    print("Getting ranks from subtrees file seqIDS")
    ranklist = []
    lines = []
    new = "no"
    newtree = "no"
    with open(subtrees) as old:
        for line in old:
            if line == "SUBTREES_BEGIN\n":
                print("beginning...")
            elif new == "yes":
                new = "no"
                newtree = "yes"
            elif newtree == "yes":
                newtree = "no"
            elif line == "subtree#\n":
                new = "yes"
            else:
                if "\n" in line:
                    line = line[:-1]
                
                if line in lines:
                    pass
                else:
                    lines.append(line)
##    print(lines)
    #also stop double tree at the front of paraphyly found trees see streptophytes
    for example in lines:
        import re
        
        gitax = re.sub("(.*)(\|)(gi#?\|?)([0-9]*)(.*)", "\\4~\\1", example)
##  this bit only works if your seqID format is taxon|omic|info|gi#|whateverelse
##if switch to accession nums, replace gi(above) with acn or something.
        
        try:
            gi, tax = gitax.split("~")
        except:
            print ("ERROR WITH:"+gitax)
            raise SystemExit
        
        tn = tax[:-1]
        taxlist = tn.split("|")
        
        therankstr = taxlist[ranknum-1]
        therankstr = therankstr+"|"
        
        if therankstr in ranklist:
            pass
        else:
            if therankstr == "X|":
                pass
            else:
                ranklist.append(therankstr)
##        print(ranklist)
    print("Added "+str(len(ranklist))+" groupnames to parsing queue")
    return ranklist
    
def MakeSubtreesDict(subtrees, verbose = False):
    print("Making Subtrees Dictionary")
    stdict = {}
    sttips = []
    st_trees = {}
    new_tree = "no"
    new = "no"
    first = "yes"
    with open(subtrees) as old:
        for line in old:
            if line == "SUBTREES_BEGIN\n":
                if verbose == True:
                    print("beginning subtree parsing")
            elif new_tree == "yes":
                st_trees[stname] = line[:-1]
                new_tree = "no"
            elif new == "yes":
                if first == "yes":
                    stname = "1"
                    first = "no"
                    new = "no"
                    new_tree = "yes"
                else:
                    #add name to dict with tips as ke
##                    print(stname+" has tips: "+str(len(sttips)))
                    stdict[stname] = sttips
                    #reset tips
                    sttips = []
                    #resetname
                    stname = line
                    new = "no"
                    if "\n" in stname:
                        stname = stname[:-1]
                    new_tree = "yes"
            elif line == "subtree#\n":
                new = "yes"
            else:
                if "\n" in line:
                    line = line[:-1]
                sttips.append(line)
        #catch the final subtre
        stdict[stname] = sttips
    print("Subtree dict finished")
##    print (stdict)
    return stdict, st_trees
    

def CheckMonophyly(string, stdict, fewnum, verbose = False):
    goodst = {}
    alltips = []
    for subtree in stdict:
        tips = stdict[subtree]
        goodtips = 0
        status = "good"
        for tip in tips:
            if string in tip:
                goodtips+=1
                if tip in alltips:
                    pass
                else:
                    alltips.append(tip)
            else:
                status = "bad"
                
        if status == "good":
            goodst[subtree] = goodtips
            #break
            ##to save time - returns fewnum here rather than calculating later.
            #BUT this means we have to loop through all tips adding to alltips here
            #wheras before, we just looped until a bad status and then broke, which was
            #much faster. unclear which way is better in long run, but this works for now.
            #negligible save in >2000 seq datasets
    if len(alltips) < fewnum:
        if verbose == True:
           print(str(len(alltips))+string+"is not enough tips to reliably assign clades")
        return("few")
    if goodst == {}:
        return "NotMono"
    #get largest monophyletic subtree
    ke, val = max(goodst.iteritems(), key=lambda x:x[1])
            #ke = max(ex.keys(), key = lambda k: ex[k])
    #ensure all tips containing string are in largest mono clade
    mono = "yes"
    for item in alltips:
        if item in stdict[ke]:
            pass
        else:
            mono = "no"
            break
    if mono == "yes":
        return (string, "Monophyly_Strict", ke, str(goodst[ke]), "0", str(len(alltips)))
    else:
        return "NotMono"

#    return "type, subtree, goodtips, cladetips, outsidetips"

def CheckMonoSingleSubtree(subtreename, subtreetipslist, alltipslist, verbose = False):
    mono = "yes"
    for item in subtreetipslist:
        if item in alltipslist:
            pass
        else:
            mono = "no"
            break
    if mono == "yes":
        return "Monophyly"
    else:
        return "NotMono"

def CheckInclusions(string, stdict, few = "no", verbose = False):
    valst = {}
    totst = {}
    goodst = {}
    alltips = []
    for subtree in stdict:
##        print("Subtree: "+subtree)
        tips = stdict[subtree]
        
        goodtips = 0
        badtips = 0
        for tip in tips:
            
            if string in tip:
                goodtips+=1
                if tip in alltips:
                    pass
                else:
                    alltips.append(tip)
            else:
                badtips+=1
        #weight subtrees by number strtips * percent of tips that are strtips
        #make dictionaries of subtreename: total tips, goodtips, %goodtips, weightedvalue
        tottips = goodtips+badtips
        totst[subtree] = tottips
        goodst[subtree] = goodtips
   
        perctips = float(goodtips)/float(tottips)
       
        valtips = goodtips*perctips
        
        valst[subtree] = valtips
##        if str(goodtips) == "0":
##            pass
##        else:
##            print("good"+str(goodtips))
##            print("bad"+str(badtips))
##            print("total"+str(tottips))
##            print("perc"+str(perctips))
##            print("Val"+str(valtips))
    #get largest value subtree.
    
    ke, val = max(valst.iteritems(), key=lambda x:x[1])
    if verbose == True:
        print("Largest subtree: "+str(ke)+" "+str(val))
            #ke = max(ex.keys(), key = lambda k: ex[k])
    #ensure all tips containing string are in largest mono(w inclusions) clade
    mono = "yes"
    for item in alltips:
        if item in stdict[ke]:
            pass
        else:
            mono = "no"
            break
    if few == "few":
        return ("few", ke, goodst, totst, alltips, valst, stdict)
    if mono == "yes":
##        print ("Subtree number "+ke+" is monophyletic with inclusions with respect to "+string+" .\nIt has "+totst[ke]+" tips total, of which "+goodst[ke]+" are "+string+".")
##set a limit on how big inclusions can be?? else we get the case where it returns subtree 1: bacteria(?)

        return (string, "Paraphyly_Strict", ke, str(goodst[ke])+"/"+str(totst[ke]), "0", str(len(alltips))) 
    else:
        return ("NotParaAlone", ke, goodst, totst, alltips, valst, stdict)
    #needs a lot of testing to ensure we always want num*percent as biggest clade. what if huge inclusion?

#artcutoff  =  (str tips in clade) / (total str tips) for determining monophyly with artifacts. eg (.7 for jumpy dataset, .9 for solid)?
#duocutoff = the percent tips-in-clade to total-str-tips a subtree must have to be considered likely true clade. eg(.1 to get many, .3 to get big ones)?ArithmeticError
def CheckArtifacts(string, stdict, monoinc, artcutoff, fewval = "no", verbose = False):
    dupcutoff = 1.0 - float(artcutoff)
    bestsubtree = monoinc[1]
    goodst = monoinc[2]
    goodtipsst = goodst[bestsubtree]
    totst = monoinc[3]
    totaltipsst = totst[bestsubtree]
    allstrtips = monoinc[4]
    valst = monoinc[5]
    stdict = monoinc[6]
    bestartval = float(goodtipsst)/float(len(allstrtips))
    if fewval == "few":
        cmono = CheckMonoSingleSubtree(bestsubtree, stdict[bestsubtree], allstrtips)
        if cmono == "Monophyly":
            return(string, "NoClade_few", "best:"+bestsubtree+"M", str(totst[bestsubtree]), str(len(allstrtips)-int(goodst[bestsubtree])), str(len(allstrtips)))
        else:
            return(string, "NoClade_few", "best:"+bestsubtree+"P", str(goodst[bestsubtree])+"/"+str(totst[bestsubtree]), str(len(allstrtips)-int(goodst[bestsubtree])), str(len(allstrtips)))
    elif bestartval > float(artcutoff):
        if goodtipsst == totaltipsst:
            if verbose == True:
                print("Mono with art")
            return(string, "Mono_Few_Exc", bestsubtree, str(goodtipsst), str(len(allstrtips)-goodtipsst), str(len(allstrtips)))
        else:
            if verbose == True:
                print("Para with art")
            return(string, "Para_Few_Exc", bestsubtree, str(goodtipsst)+"/"+str(totaltipsst), str(len(allstrtips)-goodtipsst), str(len(allstrtips)))
    else:
        if verbose == True:
            print("Not mono or para with Artifacts (cutoff="+str(artcutoff)+"), best tree ="+str(bestartval))
        numberclades = 1
        if verbose == True:
            print("dc="+str(dupcutoff))
        threshold = float(len(allstrtips))*float(dupcutoff)
        if goodst[bestsubtree] > threshold:
            pass
        else:
            if verbose == True:
                print("No clade meets min threshold of "+str(dupcutoff)+" of total group tips")
            if verbose == True:
                print("Best clade is: "+str(goodst[bestsubtree])+" but "+str(threshold)+" is needed.")
            cmono = CheckMonoSingleSubtree(bestsubtree, stdict[bestsubtree], allstrtips)
            if cmono == "Monophyly":
                return(string, "NoClade_thresh", "best:"+bestsubtree+"M", str(totst[bestsubtree]), str(len(allstrtips)-int(goodst[bestsubtree])), str(len(allstrtips)))
            else:
                return(string, "NoClade_thresh", "best:"+bestsubtree+"P", str(goodst[bestsubtree])+"/"+str(totst[bestsubtree]), str(len(allstrtips)-int(goodst[bestsubtree])), str(len(allstrtips)))
        cladeslist = [bestsubtree]
        itis = "on"
        notallowed = stdict[bestsubtree]
        if verbose == True:
            print("Best st is: "+str(bestsubtree)+" with tips: "+str(len(notallowed)))
        liststdict = []
        for i in stdict:
            exi = i
            liststdict.append(i)
        if verbose == True:
            print("There are "+str(len(liststdict))+" distinct subtrees to start")
        while itis == "on":
            #remove subtrees from consideration if they contain tips from already considered subtrees.
            if verbose == True:
                print("Removing intersecting subtrees based on "+str(len(notallowed))+" used tips")
            sttotal = 0
            stremoved = 0
            lis2 = []
            for item in liststdict:
                lis2.append(item)
            for subtreel in liststdict:
                sttotal +=1
                ok = "yes"
                for tip in stdict[subtreel]: 
                    for ntip in notallowed:
                        if tip == ntip: 
                            ok = "no"
##                            print(valst[subtreel])
##                            print("removing a subtree named"+str(subtreel)+"of type"+str(type(subtreel)))
                            if subtreel in lis2:
                                stremoved +=1
                                lis2.remove(subtreel)
                                del valst[subtreel]
                            break
                    if ok == "no":
                        break
            liststdict = lis2
            if verbose == True:
                print("Start: "+str(sttotal)+" Removed: "+str(stremoved)+" Considering "+str(len(liststdict))+" exclusive subtrees")
            #if all trees contain tip, and all trees overlap, valst will be entirely removed. skipp further verification, will return one cladee tree later. later
            if len(valst) == 0:
                itis = "off"
            else:
                #get topmost non-intersecting tree
                ke, val = max(valst.iteritems(), key=lambda x:x[1])
                #see if it beats threshold to be considered a clade in its own right
                #if yes, continue search. if no close out.
                if goodst[ke] > threshold:
                    numberclades+=1
                    if verbose == True:
                        print("Found a clade. Num "+str(ke)+"Val: "+str(valst[ke])+"Good: "+str(goodst[ke])+"Thr: "+str(threshold)+"total: "+str(numberclades))
                    cladeslist.append(ke)
                    if verbose == True:
                        print("Adding "+str(len(stdict[ke]))+" tips to notallowed from subtree "+str(ke))
    ##                print (stdict[ke])
                    for eachtip in stdict[ke]:
                        valst[ke] = 0
                        if eachtip in notallowed:
                            pass
                        else:
                            notallowed.append(eachtip)
    ##                print("not allowed:"+str(len(notallowed)))
                else:
                    itis = "off"
        #create info for table
        if verbose == True:
            print("total clades: "+str(numberclades))
        goodstl = ""
        totstl = ""
        trees = "["+str(numberclades)+"]: "
        gatip = 0
        gtot = ""
        for thethings in notallowed:
            if string in thethings:
                gatip+=1
        #return if one clade found, too messy for strict para/mono with artifacts.
        if str(numberclades) == "1":
            if verbose == True:
                print("Not polyphyletic: only one clade discovered(clade cutoff="+str(dupcutoff)+")")
            cmono = CheckMonoSingleSubtree(cladeslist[0], stdict[cladeslist[0]], allstrtips)
            if cmono == "Monophyly":
                
                return(string, "Monophyly_Many_Exc", cladeslist[0], str(totst[cladeslist[0]]), str(len(allstrtips)-goodst[cladeslist[0]]), str(len(allstrtips)))
            else:
                return(string, "Paraphyly_Many_Exc", cladeslist[0], str(goodst[cladeslist[0]])+"/"+str(totst[cladeslist[0]]), str(len(allstrtips)-goodst[cladeslist[0]]), str(len(allstrtips)))
        else:
            for item in cladeslist:
                cmono = CheckMonoSingleSubtree(item, stdict[item], allstrtips)
                if cmono == "Monophyly":
                    gtot = gtot+str(goodst[item])+", "
                    trees = trees+item+"M, "
                else:
                    gtot = gtot+str(goodst[item])+"/"+str(totst[item])+", "
                    trees = trees+item+"P, "
            trees = trees[:-2]
            gtot = gtot[:-2]   
            if verbose == True:
                print("Polyphyletic")
            return(string, "Polyphyly", trees, gtot, str(len(allstrtips)-gatip), str(len(allstrtips)))

def PerformScan(stringlist, treefile, artcutoff, fewnum, treesubfile, verbose = False):
    if stringlist[0] == "PerformGetRankFromSubtreesFile":
        ranknumstring = stringlist[1]
        
        ranknumlist = ranknumstring.split()
        perf = "yes"
    else:
        perf = "no"
    subtreefilename = treefile+"Subtrees.txt"
    #this bit will currently only work on my laptop (Artemis)
    #requires the following R script, which requires ape library installed:
###!/usr/bin/Rscript
###get arguments in a list, first nexus file then destination file
##myArgs <- commandArgs(trailingOnly = TRUE)
###print(myArgs)
##nexusfile<-myArgs[1]
##print("opening")
##print(nexusfile)
##subtreefile<-myArgs[2]
###create new empty destination file
##write("SUBTREES_BEGIN", file = subtreefile, append = FALSE)
###import necessary packages
##library(ape)
###open nexus file
##tree<-read.nexus(nexusfile)
##stree<-subtrees(tree)
##print("writing simple subtrees to file: ")
##print(subtreefile)
###write number of subtree followed by its tips, each on a newline.
##for(i in 1:length(stree)){
##  write(c("subtree#",i), file = subtreefile, append = TRUE)
##  write.tree(stree[[i]], file = subtreefile, append = TRUE)
##  for (n in 1:length(stree[[i]]$tip.label)){
##  write(stree[[i]]$tip.label[[n]], file = subtreefile, append = TRUE)
##  }}
##print("finished!")

    #with alias or symlink to keyword "subtreegen"
    if treesubfile == "NA":
        
        subtreefilename = treefile+"Subtrees.txt"
        subtreefilename = MakeSubtreesFile(treefile, subtreefilename, verbose)
    else:
        subtreefilename = treesubfile
    if perf == "yes":
        stringlist = []
##        print(ranknumlist)
        for rnum in ranknumlist:
##            print(rnum)
            stringlist.append("RANKLEVEL: "+str(rnum))
            ranknum = int(rnum)
            rstringlist = GetRankFromSubtreesFile(subtreefilename, ranknum, verbose)
            for rstring in rstringlist:
                stringlist.append(rstring)
    stdict, st_trees = MakeSubtreesDict(subtreefilename, verbose)
    from prettytable import PrettyTable
    tablelist = []
    rankinfolist = []
    start = "yes"
    table = PrettyTable(['Group', 'Topology', 'Subtree#', 'ST_Tips', 'Excluded', "Total"])
    print("Checking Mono/Para/Polyphyly for each given group")
    for string in stringlist:
        if "RANKLEVEL:" in string:
            table.sortby = "Topology"
            
            if start == "yes":
                start = "no"
            else:
                tablelist.append(table)
                table = PrettyTable(['Group', 'Topology', 'Subtree#', 'ST_Tips', 'Excluded', "Total"])
            tablelist.append(string)
            continue
        
        if verbose == True:
            print("Checking group: "+string)
        CM = CheckMonophyly(string, stdict, fewnum, verbose)
        tab = CM
        if CM == "NotMono":
            if verbose == True:
                print("Not monophyletic")
            CI = CheckInclusions(string, stdict, verbose)
            tab = CI
            if CI[0] == "NotParaAlone":
                if verbose == True:
                    print("Not paraphyletic")
                CA = CheckArtifacts(string, stdict, CI, artcutoff, verbose)
                tab = CA
        if CM == "few":
            CI = CheckInclusions(string, stdict, "few")
            CA = CheckArtifacts(string, stdict, CI, artcutoff, "few", verbose)
            tab = CA

        table.add_row(tab)
        #whichever is good, print(?) or just return I guess if we're gunna run though lots
    table.sortby = "Topology"
    tablelist.append(table)
   
##    for titem in tablelist:
##        print titem
    return tablelist

def GetRankFromSubtreesFile_SubSample(subtrees, ranknum, groupstring, verbose = False):
    print("Getting ranks from subtrees file seqIDS")
    ranklist = []
    lines = []
    new = "no"
    newtree = "no"
    with open(subtrees) as old:
        for line in old:
            if line == "SUBTREES_BEGIN\n":
                print("beginning...")
            elif new == "yes":
                new = "no"
                newtree = "yes"
            elif newtree == "yes":
                newtree = "no"
            elif line == "subtree#\n":
                new = "yes"
            else:
                if "\n" in line:
                    line = line[:-1]
                
                if line in lines:
                    pass
                else:
                    lines.append(line)
    for example in lines:
        if groupstring in example:
            pass
        else:
            continue
        import re     
        gitax = re.sub("(.*)(\|)(gi#?\|?)([0-9]*)(.*)", "\\4~\\1", example)
##  this bit only works if your seqID format is taxon|omic|info|gi#|whateverelse
##if switch to accession nums, replace gi(above) with acn or something.
        
        gi, tax = gitax.split("~")
        
        tn = tax[:-1]
        taxlist = tn.split("|")
        
        therankstr = taxlist[ranknum-1]
        therankstr = therankstr+"|"
        
        if therankstr in ranklist:
            pass
##        else:
##            if therankstr == "X|":
##                pass
##            else:
##                ranklist.append(therankstr)
        else:
            ranklist.append(therankstr)
##        
##        print(ranklist)
##    print("Added "+str(len(ranklist))+" groupnames (subset of "+groupstring+" to parsing queue")
    return ranklist

def PerformScan_SubSample(stringlist, treefile, artcutoff, fewnum, treesubfile, verbose = False):
    print("Beginning outer clade-determination")
    if stringlist[0] == "PerformGetRankFromSubtreesFile":
        ranknumstring = stringlist[1]
        
        ranknumlist = ranknumstring.split()
        perf = "yes"
    else:
        perf = "no"
    subtreefilename = treefile+"Subtrees.txt"
    #this bit will currently only work on my laptop (Artemis)
    if treesubfile == "NA":
        subtreefilename = treefile+"Subtrees.txt"
        subtreefilename = MakeSubtreesFile(treefile, subtreefilename, verbose)
    else:
        subtreefilename = treesubfile
    if perf == "yes":
        stringlist = []
##        print(ranknumlist)
        for rnum in ranknumlist:
##            print(rnum)
            stringlist.append("RANKLEVEL: "+str(rnum))
            ranknum = int(rnum)
            rstringlist = GetRankFromSubtreesFile(subtreefilename, ranknum, verbose)
            for rstring in rstringlist:
                stringlist.append(rstring)
    stdict, st_trees = MakeSubtreesDict(subtreefilename, verbose)
    from prettytable import PrettyTable
    tablelist = []
    rankinfolist = []
    start = "yes"
    table = PrettyTable(['Group', 'Topology', 'Subtree#', 'ST_Tips', 'Excluded', "Total"])
    subtrees_new_dict = {}
##    print("Checking Mono/Para/Polyphyly for each given group")
    for string in stringlist:
        if "RANKLEVEL:" in string:
            table.sortby = "Topology"
            
            if start == "yes":
                start = "no"
            else:
                tablelist.append(table)
                table = PrettyTable(['Group', 'Topology', 'Subtree#', 'ST_Tips', 'Excluded', "Total"])
            tablelist.append(string)
            continue
        
##        print("Checking group: "+string)
        CM = CheckMonophyly(string, stdict, fewnum, verbose)
        tab = CM
        if CM == "NotMono":
##            print("Not monophyletic")
            CI = CheckInclusions(string, stdict, verbose)
            tab = CI
            if CI[0] == "NotParaAlone":
##                print("Not paraphyletic")
                CA = CheckArtifacts(string, stdict, CI, artcutoff, verbose)
                tab = CA
        if CM == "few":
            CI = CheckInclusions(string, stdict, "few")
            CA = CheckArtifacts(string, stdict, CI, artcutoff, "few", verbose)
            tab = CA

        table.add_row(tab)
        subtrees_line = tab[2]
      
##        print(subtrees_line)
        subtrees_line_list = subtrees_line.split(",")
##        print(subtrees_line_list)
        for item in subtrees_line_list:
            if "[" in item:
                item = item[3:]
            nitem = re.sub("[A-Za-z \[\]:]*", "", item)
            if len(nitem) >4:
                print (nitem+"   "+item)
##            print(nitem)
            if nitem in subtrees_new_dict:
                 subtrees_new_dict[nitem] =  subtrees_new_dict[nitem]+" "+string
            subtrees_new_dict[nitem] = string
            #subtrees new dict contains 1,2,3 (subtree numbers) : actinobacteria (or "dogs cat rats")?
    
        #whichever is good, print(?) or just return I guess if we're gunna run though lots
    table.sortby = "Topology"
    tablelist.append(table)
   
##    for titem in tablelist:
##        print titem
##    print("SubNewDict")
####    print(subtrees_new_dict)
    return tablelist, subtrees_new_dict, st_trees


def PerformScan_SubSample_Inner(stringlist, treefile, artcutoff, fewnum, treesubfile, groupstring, tree_in_text, verbose = False):
    
    subtrees_new_dict = {}
    if stringlist[0] == "PerformGetRankFromSubtreesFile":
        ranknumstring = stringlist[1]
        ranknumlist = ranknumstring.split()
        newranknumlist = []
        for item in ranknumlist:
            item = int(item) + 1
            newranknumlist.append(item)
        ranknumlist = newranknumlist
        perf = "yes"
    else:
        perf = "no"
    subtreefilename = treefile+"Subtrees.txt"
    #this bit will currently only work on my laptop (Artemis)
    if treesubfile == "NA":
        subtreefilename = treefile+"Subtrees.txt"
        subtreefilename = MakeSubtreesFile(treefile, subtreefilename, verbose)
    else:
        subtreefilename = treesubfile
    if perf == "yes":
        stringlist = []
##        print(ranknumlist)
        for rnum in ranknumlist:
##            print(rnum)
            stringlist.append("RANKLEVEL: "+str(rnum))
            ranknum = int(rnum)
            rstringlist = GetRankFromSubtreesFile_SubSample(subtreefilename, ranknum, groupstring, verbose)
            for rstring in rstringlist:
                stringlist.append(rstring)
    print(groupstring)
    #stdict is ??
    #st_trees is number:newickstr
    stdict, st_trees = MakeSubtreesDict(subtreefilename, verbose)
    from prettytable import PrettyTable
    tablelist = []
    rankinfolist = []
    start = "yes"
    table = PrettyTable(['Group', 'Topology', 'Subtree#', 'ST_Tips', 'Excluded', "Total"])
    for string in stringlist:
        print("checking :"+string)
        if "RANKLEVEL:" in string:
            table.sortby = "Topology"
            if start == "yes":
                start = "no"
            else:
                tablelist.append(table)
                table = PrettyTable(['Group', 'Topology', 'Subtree#', 'ST_Tips', 'Excluded', "Total"])
            tablelist.append(string)
            continue
##        print("Checking group: "+string)
        CM = CheckMonophyly(string, stdict, fewnum, verbose)
        tab = CM
        if CM == "NotMono":
##            print("Not monophyletic")
            CI = CheckInclusions(string, stdict, verbose)
            tab = CI
            if CI[0] == "NotParaAlone":
##                print("Not paraphyletic")
                CA = CheckArtifacts(string, stdict, CI, artcutoff, verbose)
                tab = CA
        if CM == "few":
            CI = CheckInclusions(string, stdict, "few")
            CA = CheckArtifacts(string, stdict, CI, artcutoff, "few")
            tab = CA

        table.add_row(tab)
        subtrees_line = tab[2]
    
        subtrees_line_list = subtrees_line.split(",")
        for item in subtrees_line_list:
            if "[" in item:
                item = item[3:]
            else:     
                nitem = re.sub("[A-Za-z :\[\]]*", "", item)
                if nitem in subtrees_new_dict:
                   
                    subtrees_new_dict[nitem] =  subtrees_new_dict[nitem]+" "+string
                subtrees_new_dict[nitem] = string
        print("Inner subtree determined to be:"+nitem)
        
        #whichever is good, print(?) or just return I guess if we're gunna run though lots
    table.sortby = "Topology"
    tablelist.append(table)
   
##    for titem in tablelist:
##        print titem
    #if only one order found in class-group-search (eg, all Halos are Natrialbales) we should return two halo sequences, one basal and one random.
    #else we only return basal seq of all orders represented within class-clade.
    if len(subtrees_new_dict) == 1:
        subtrees_new_dict["alone"] = string
        st_trees["alone"] = tree_in_text
    return tablelist, subtrees_new_dict, st_trees

def LeastIndentTip(newicktree, alone, groupstring, verbose=False):
    indent = 0
    newtip = "yes"
    tipend = "no"
    ongoingtip = "no"
    tipdict = {}
    for item in newicktree:
        if item == "(":
            indent+=1
        elif item == ")":
            indent-=1
        elif item == ",":
            pass
        elif item == ":":
            ongoingtip = "no"
            if groupstring in currenttip:
                tipdict[currenttip] = indent
            else:
                pass
        elif ongoingtip == "no":
            if item == "A":
                newtip = "yes"
                currenttip = item
                ongoingtip = "yes"
            elif item == "B":
                newtip = "yes"
                currenttip = item
                ongoingtip = "yes"
            elif item == "E":
                newtip = "yes"
                currenttip = item
                ongoingtip = "yes"
            elif item == "X":
                newtip = "yes"
                currenttip = item
                ongoingtip = "yes"
            else:
                pass
        elif ongoingtip == "yes":
            currenttip = currenttip+item
        else:
            print("Error in indentcounting")
    if len(tipdict) == 0:
        print("Zero!"+groupstring)
        print(newicktree)
        return("none")
    if alone == "yes":
        import random
##        print (newicktree)
##        print (tipdict)
        shallowtip = random.choice(tipdict.keys())
    else:  
        shallowtip, value = min(tipdict.iteritems(), key=lambda x:x[1])
    if groupstring in shallowtip:
            return shallowtip
    else:
        print("groupstring not in chosen tip")
        raise SystemExit
                
                
        
##________________
##EASYMODE: subsampling
##0. PreformScan_SubSample -> subtrees_new_list (of all subtrees found at rank (class)
##1. for each item in subtrees_new_list, get subtrees_tree and run easymode:

def SubSampling_Master(stringlist, treefile, artcutoff, fewnum, treesubfile, inputfasta, verbose):
    tablelist, subtrees_new_dict, st_trees = PerformScan_SubSample(stringlist, treefile, artcutoff, fewnum, treesubfile, verbose)
    chosen_seqs = []
##    print(subtrees_new_dict)
    import random
    print("Clade - dictionaries were successfully created. Moving on to inner SS... ")
    with open(treefile+"SubSamp_INFO.txt", "w") as info:
        for item in subtrees_new_dict:
            groupstring = subtrees_new_dict[item]
            a = st_trees[item]
            with open ("Temp"+item+".txt", "w") as new:
                 new.write(a)
            dire = os.getcwd()
            os.system("ConvertPhylo "+dire+" Temp"+item+".txt newick nexus")
            nexus = "Temp"+item+".nexus"
            #run the thing on the nexus file with given string(s)
            inner_treefile=nexus
            inner_treesubfile="NA"
            inner_stringlist=stringlist
            
            info.write("\nfor group "+groupstring+" subtree "+item+" the following sequences have been subsampled:")
            inner_tablelist, inner_subtrees_new_dict, inner_st_trees = PerformScan_SubSample_Inner(inner_stringlist, inner_treefile, artcutoff, fewnum, inner_treesubfile, groupstring, a, verbose)
##            if verbose == True:
            print("Created inner clade - dictionaries for group: "+groupstring+" consisting of "+str(len(inner_subtrees_new_dict))+" subtrees")
            for tree in inner_subtrees_new_dict:
                inner_string = inner_subtrees_new_dict[tree]
                if tree == "alone":
                    info.write("\nOnly one subgroup found... keeping an additional random sequence.")
                    chosen_sequence = LeastIndentTip(inner_st_trees[tree], "yes", groupstring, verbose)
                else:
                    chosen_sequence = LeastIndentTip(inner_st_trees[tree], "no", groupstring, verbose)
                if chosen_sequence in chosen_seqs:
                    
                    info.write("\n"+inner_string+": "+chosen_sequence+" is a Duplicate... avoiding")
                elif chosen_sequence == "none":
                    pass
                else:
                    chosen_seqs.append(chosen_sequence)
                    info.write("\n"+inner_string+": "+chosen_sequence)
            #remove generated files
            os.system("rm "+nexus)
            os.system("rm Temp"+item+".txt")        
            subtreefilename = nexus+"Subtrees.txt"
            os.system("rm "+subtreefilename)
        info.write("\nFinished chosing extracted sequences will be at: "+inputfasta+"Extracted.fasta\nInfo is available at: "+treefile+"SubSamp_INFO.txt")

    print("Finished choosing sequences to keep.")
    print("Extracting given sequences.... ")
    keepseqs = ""
    for seqid in chosen_seqs:
        keepseqs = keepseqs+" "+seqid
    keepseqs=keepseqs[:-1]
    print(inputfasta)
    print(keepseqs)
    os.system("feast . "+inputfasta+" -ex \""+keepseqs+"\"")
    print("Finished, extracted sequences SHOULD be at: "+inputfasta+"Extracted.fasta")
    print("Info available at: "+treefile+"SubSampInfo.txt")
   





############september###########
class Subtree:
    def __init__(self, name):
    #all ids should be stripped and have ">" removed for reasons.
    #for now, sequences do not have any stripping applied
        self.number_name = name
        self.string_name = ""
        self.category = ""
        self.tips = []
        self.fasta = ""
        self.fasttree_st = ""
        self.raxml = ""
        self.species_tree = ""
    def set_species_tree(self,speciest):
        self.species_tree = ""
    def set_raxml(self, rax):
        self.raxml = rax
    def set_fasttree(self, ft):
        self.fasttree_st=ft
    def set_string(self, stringn):
        self.string_name = stringn
    def set_category(self, cat):
        self.category = cat
    def set_fasta(self, fasta):
        self.fasta = fasta
    def set_tips(self, tips):
        if type(tips) == str:
            self.tips.append(str)
        else:
            for item in tips:
                self.tips.append(item)
    def ret_string(self):
        return self.string_name
    def ret_cat(self):
        return self.category
    def ret_tips(self):
        return self.tips
    def ret_fasta(self):
        return self.fasta
    def ret_fasttree(self):
        return self.fasttree_st
    def ret_raxml(self):
        return self.raxml
    def ret_species_tree(self):
        return self.species_tree



def PerformScan_Complex_SubSample(stringlist, treefile, artcutoff, fewnum, treesubfile, verbose = False):
    print("Beginning outer clade-determination")
    if stringlist[0] == "PerformGetRankFromSubtreesFile":
        ranknumstring = stringlist[1]
        
        ranknumlist = ranknumstring.split()
        perf = "yes"
    else:
        perf = "no"
    subtreefilename = treefile+"Subtrees.txt"
    print(treesubfile)
    print(ranknumlist)
    #this bit will currently only work on my laptop (Artemis)
    if treesubfile == "NA":
        subtreefilename = treefile+"Subtrees.txt"
        subtreefilename = MakeSubtreesFile(treefile, subtreefilename, verbose)
    else:
        subtreefilename = treesubfile
    if perf == "yes":
        stringlist = []
##        print(ranknumlist)
        for rnum in ranknumlist:
##            print(rnum)
            stringlist.append("RANKLEVEL: "+str(rnum))
            ranknum = int(rnum)
            rstringlist = GetRankFromSubtreesFile(subtreefilename, ranknum, verbose)
            for rstring in rstringlist:
                stringlist.append(rstring)
    stdict, st_trees = MakeSubtreesDict(subtreefilename, verbose)
    from prettytable import PrettyTable
    tablelist = []
    rankinfolist = []
    start = "yes"
    table = PrettyTable(['Group', 'Topology', 'Subtree#', 'ST_Tips', 'Excluded', "Total"])
    subtrees_new_dict = {}
##    print("Checking Mono/Para/Polyphyly for each given group")
    for string in stringlist:
        if "RANKLEVEL:" in string:
            table.sortby = "Topology"
            
            if start == "yes":
                start = "no"
            else:
                tablelist.append(table)
                table = PrettyTable(['Group', 'Topology', 'Subtree#', 'ST_Tips', 'Excluded', "Total"])
            tablelist.append(string)
            continue
        
##        print("Checking group: "+string)
        CM = CheckMonophyly(string, stdict, fewnum, verbose)
        tab = CM
        if CM == "NotMono":
##            print("Not monophyletic")
            CI = CheckInclusions(string, stdict, verbose)
            tab = CI
            if CI[0] == "NotParaAlone":
##                print("Not paraphyletic")
                CA = CheckArtifacts(string, stdict, CI, artcutoff, verbose)
                tab = CA
        if CM == "few":
            CI = CheckInclusions(string, stdict, "few")
            CA = CheckArtifacts(string, stdict, CI, artcutoff, "few", verbose)
            tab = CA

        table.add_row(tab)
        subtrees_line = tab[2]
      
##        print(subtrees_line)
        subtrees_line_list = subtrees_line.split(",")
##        print(subtrees_line_list)
        for item in subtrees_line_list:
            if "[" in item:
                item = item[3:]
            nitem = re.sub("[A-Za-z \[\]:]*", "", item)
            if len(nitem) >4:
                print (nitem+"   "+item)
##            print(nitem)
            if nitem in subtrees_new_dict:
                 subtrees_new_dict[nitem] =  subtrees_new_dict[nitem]+" "+string
            subtrees_new_dict[nitem] = string
            #subtrees new dict contains 1,2,3 (subtree numbers) : actinobacteria (or "dogs cat rats")?
    
        #whichever is good, print(?) or just return I guess if we're gunna run though lots
    table.sortby = "Topology"
    tablelist.append(table)
   
##    for titem in tablelist:
##        print titem
##    print("SubNewDict")
####    print(subtrees_new_dict)
    return tablelist, subtrees_new_dict, st_trees, stdict

    
def ComplexSubSamplingMaster(stringlist, treefile, artcutoff, fewnum, treesubfile, inputfasta, verbose):
    tablelist, subtrees_new_dict, st_trees, stdict = PerformScan_Complex_SubSample(stringlist, treefile, artcutoff, fewnum, treesubfile, verbose)
    print(stdict)
    print(tablelist)
    print(subtrees_new_dict)
    print(st_trees)
    chosen_seqs = []
#here i want to make a bunch of objects of class Subtree, one for each idntified clade of interest.




    import random
    print("Clade - dictionaries were successfully created. Moving on to inner SS... ")
    with open(treefile+"SubSamp_INFO.txt", "w") as info:
        for item in subtrees_new_dict:
            groupstring = subtrees_new_dict[item]
            print(groupstring)
            a = st_trees[item]
            #write the .fasttree subtree to file. item = groupname
            with open ("Fasttree_subtree_"+groupstring+".newick", "w") as new:
                 new.write(a)
            dire = os.getcwd()
            
            break
        
    print("Finished choosing sequences to keep.")
    print("Extracting given sequences.... ")
    keepseqs = ""
    for seqid in chosen_seqs:
        keepseqs = keepseqs+" "+seqid
    keepseqs=keepseqs[:-1]
    print(inputfasta)
    print(keepseqs)
    os.system("feast . "+inputfasta+" -ex \""+keepseqs+"\"")
    print("Finished, extracted sequences SHOULD be at: "+inputfasta+"Extracted.fasta")
    print("Info available at: "+treefile+"SubSampInfo.txt")
   
##this does subsampling by making raxml and species trees for each subtree of given rank, and then







    
##
##    (where input is one level down, and tree given is actually a subtree)                                            
##
##    2. run CladeParser at order-level within best subtree
##    3. for each found clade, keep least indented tip sequence (given that it is member of group)
##def EASYMODE(subtree_tree, subtree_tips):
##    
##
##    subtree_file_name = ???
##    remove treefile
    
#something to either 1. be given a specific string (or list of strings) and run
#or
#be given a depth to search eg >Superkingdom|phylum|class|species_name|gi|etc
#preform search on all classes in entire file
##
##parser will have
##1. treefile
##2. stringlist or depth-search
##
##stringlist = ["PerformGetRankFromSubtreesFile", ranknum]
## or  = [actual, list, of, strings]

                   
 #############################parser################333                   
if __name__ == "__main__":

    print("Running in terminal")  
    import sys
    import argparse
    import os
    import re
    parser = argparse.ArgumentParser(description="All")
    parser.add_argument("directory", nargs='?', default=os.getcwd(), type=str, help="type name of directory to run in (where .nex resides)")
    parser.add_argument("NEXUS_TREE", type=str, help="type the name of your .nex file")
    parser.add_argument("-st", "--string", action = "store", help="to parse only given string, or strings space-seperated within quotes")
    parser.add_argument("-r", "--rank", action = "store", help="to parse all strings of given depth in tree. Euk|Meta|Chor|Mam = 1|2|3|4. Type one number or several space-sep in quotes.")
    parser.add_argument("-p", "--polyclade", type = float, action = "store", default = 0.1, help = "grouptips(in-clade/in-tree) required in subtree to count as a distinct clade. default is 0.1")
    parser.add_argument("-m", "--minimum", action = "store", default = 10, help="number of tips below which a group will not be considered. default is 10")
    parser.add_argument("-t", "--treesfile", action = "store", default = "NA", help="if you already have a subtrees file and dont want to rerun it")
    parser.add_argument("-ss", "--subsampling", action = "store_true", help="toggles one rank lower subsampling within each clade")
    parser.add_argument("-c", "--complexss", action = "store_true", help="toggles complex subsampling within each clade")
    parser.add_argument("-f", "--fasta", action = "store", help="if using subsampling, provide a fasta file to pull sequences from")
    parser.add_argument("-v", "--verbose", action = "store_true", help="prints more information - for debugging mostly.")
    
    args = parser.parse_args()
    
#making input list
    try:
        os.chdir(args.directory)
    except:
        print ("didn't change dir")
    if ".newick" in args.NEXUS_TREE:
        try:
            os.system("ConvertPhylo . "+args.NEXUS_TREE+" newick nexus")
            name, newick = args.NEXUS_TREE.split(".")
            args.NEXUS_TREE = name+".nexus"
        except:
            pass
    if str(args.string) == "None":
        if str(args.rank) == "None":
            print("You need to either specify a string with -s, or a rank-depth with -r")
            raise SystemExit
        arank = args.rank.split()
        stringlist = ["PerformGetRankFromSubtreesFile", args.rank]
    elif str(args.rank) == "None":
        stringlist = args.string.split()
    else:
        print("You may not use both -s and -r at the same time (yet), sorry!")
    artcutoff = 1.0 - float(args.polyclade)

    #complex subsampling needs to exist. creates a .raxml of the subtree and also makes a raxml species tree for those taxa.
    if args.subsampling == True:
        print("Beginning Smart SubSampling!")
        SubSampling_Master(stringlist, args.NEXUS_TREE, artcutoff, int(args.minimum), args.treesfile, args.fasta, args.verbose)
        #call master subsampling
    if args.complexss == True:
        ComplexSubSamplingMaster(stringlist, args.NEXUS_TREE, artcutoff, int(args.minimum), args.treesfile, args.fasta, args.verbose)
    else:
        tabl = PerformScan(stringlist, args.NEXUS_TREE, artcutoff, int(args.minimum), args.treesfile, args.verbose)      
        information = args.NEXUS_TREE+"_INFO.txt"
        print("Analysis Finished!")
        print("************************************************************")
        for item in tabl:
            print(item)
        with open (information, "w") as info:
            for titem in tabl:
                if type(titem) == str:
                    info.write("\n"+titem+"\n")
                else:
                    titem_text = titem.get_string()
                    info.write(titem_text)
            info.write(str(args))
        print("finished")

##*****W_art should check for polyphyly - could have 1 @ 80, 1 @ 20.
##should check:
##    mono
##    para_strict - should have an inclusions limit
##    para_some_excl
##    mono_some_excl  
##    polyphyly or oneclade_mess
##      Poly: check each st for mono/para. report "Polyphyly | [4] 82M, 32P, 89M, 4P  | 10/10 5/10 8/8 5/7"
##if mono, reply only the number in the clade (will match)
##        Onemess_ check mono or para, then return Mono_many_excl
##    TooFewTips -> check mono on largest 
##    NoClade@Threshold -> check 
##    
####    
##"""
##Definitions:
##Monophyly        - The best clade contains only group-tips and no other-tips
##Paraphyly        - The best clade contains group-tips and also some other-tips. Currently no restriction on
##                   number of inclusions is implemented.
##--Strict         - Applied to a mono or para best clade. Indicates that ALL group-tips in the tree are within
##                   that specific clade. 
##--Some_Exc       - Applied to a mono or para best clade, indicates that some group-tips
##                   (up to 10%) are found outside of the best clade. Best-clade must
##                   contain at least 90% of group-tips to return this modifier.
##--Many_Exc       - Applied to Mono or Para, indicates that though > 10% of group-tips are found outside of the
##                   best clade, they are not grouped such that there is any other clade containing at least 10%
##                   of all group tips.
##Polyphyly        - If > 1 clade contains at least 10% of group-tips, polyphyly will be returned. Each subtree
##                 - will be checked for mono/para status, and an M or P printed next to it's number.
##NoClade_few      - There were not enough tips (>10) containing the given string to reach consensus. 
##NoClade_thresh   - No clade was found to contain at least 10% of all group-tips.
##Note: 10% clade cutoff can be changed with args.polyclade;10 tips for consensus changed with args.fewnum.
##"""
##        
##    
