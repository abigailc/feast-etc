#!/usr/bin/python
#this is a conglomerate of fasta-fixing scripts, now called FISH (FASTA ID SWAPPING HELPER) because lol i can acronym.
##
#last edit abigailc@Actaeon Sept 7 2016

#things this doesn't do: play super nice with accession numbers instead of GI numbers. probably easy to convert, (see that one script that one time), but meh
#do it later.

#todo
#import functions from FEAST - append taxonomy, extract, split, etc


#when you create a Fasta object, initialize it with sequence ID and Data by using either gen_original_lists or blast2fasta
class Fasta:
    def __init__(self, name):
        #all ids and seqs should be stripped of leading and trailing whitespace and have ">" removed for reasons.
        #this is the name of the fasta, it can be anything, i'm not really using it right now.
        self.name = name
        #this is the to-be-modified version of sequence IDs and sequence-Data
        # ALWAYS keep IDS and SEQS the same length. id[1] should ALWAYS correspond to seq[1].
        self.ids = []
        self.seqs = []
        # these are the original SEQids and Sequences. They should never be modified after generation in gen_original_lists or blast_to_fasta
        self.original_ids = []
        self.original_seqs = []
        self.species_names = []
        self.numbers = []
        self.taxid = []
        self.taxonomy = []
    def ret_name(self):
        return self.name
    def gen_original_lists(self, fastaname):
        with open(fastaname) as fastafile:
            for line in fastafile:
                if "\n" == line:
                    pass
                if ">" in line:
                    #write the previous AA seq
                    try:
                        AAseq=AAseq.strip()
                        self.seqs.append(AAseq)
                        self.original_seqs.append(AAseq)
                    except:
                        pass
                        #initialize a new AAseq
                    AAseq = ""
                    #format the seqID
                    newline = line.strip()
                    newline = line.strip(">")
                    #write the seqID
                    self.ids.append(newline.strip())
                    self.original_ids.append(newline.strip())
                else:
                    AAseq = AAseq+line
                    AAseq=AAseq.strip()
            #catch the last AAseq pass
            self.seqs.append(AAseq)
            self.original_seqs.append(AAseq)
        print("Initial sequence and ID lists created. Contains "+str(len(self.ids))+" sequences")
    def manual_shorten(self, shorts):
        #the list of shorts will be provided like "Bacteria,Bac Eukarya,Euk"
        changes = shorts.split()
        for item in self.ids:
            newline = item
            index = self.ids.index(item)
            for change in changes:
                old, new = change.split(",")
                newline = newline.replace(old,new)
            self.ids[index] = newline
            #done
        print("Manual shorten complete")
    def gen_numbers(self):
        for item in self.ids:
            number = re.sub("(.*)(\|)(.*)","\\3", item)
            self.numbers.append(number)
    def gen_species_lists(self):
        self.species_names = []
        speclist = []
        for item in self.ids:
                                # item will be "Nostoc_punctiforme_PCC_73102|gi#|186468349" or "Blah|Rank|Nostoc_punctiforme_PCC_73102|gi#|186468349"
                                # for now, ignores anything that isn't Genus_species.
                                # for example, ignores strain, ., things with an extra
                                # word, etc.
            taxon = re.sub("([^_]*)([A-Z][a-z]*_[a-z]*)(.*)", "\\2", item)
            if "#" in taxon:
                print ("TAXON error in gen_species_lists():" + taxon)
            speclist.append(taxon)
            self.species_names.append(taxon)
        return speclist
    def common_shorten(self, verbose = False):
        #TODO: allow input of manual shorten-pairs, possibly in new function
        #put your conversions of common strings to shorten here
        inte = 0
        for item in self.ids:
            newline = item
            index = self.ids.index(item)
            newline = re.sub("bacteria\|", "bac|", newline)
            newline = re.sub("bacteriales\|", "bacl|", newline)
            newline = re.sub("bacteriaceae\|", "bacc|", newline)
            newline = re.sub("Bacteria\|", "Bac|", newline)
            newline = re.sub("Archaea\|", "Arc|", newline)
            newline = re.sub("Eukaryota\|", "Euk|", newline)
            newline = re.sub("Fungi\|", "Fun|", newline)
            newline = re.sub("Viridiplantae\|", "Vir|", newline)
            newline = re.sub("Metazoa\|", "Met|", newline)
            newline = re.sub("mycetes\|", "myc|", newline)
            newline = re.sub("mycetales\|", "mycl|", newline)
            newline = re.sub("mycetaceae\|", "mycc|", newline)
            newline = re.sub("Methanomassiliicoccaceae\|", "Methmasscoc|", newline)
                     
            #newline = re.sub("bacteriales\|", "bacles|", newline)
            #newline = re.sub("bacteriales\|", "bacles|", newline)
            #newline = re.sub("[+=\.]", "", newline)
            newline = re.sub("_enterica_subsp_enterica_serovar", "", newline)
            if newline == item:
                pass
            else:
                if verbose is True:
                    print(item)
                    print(newline)
                    inte +=1
            self.ids[index] = newline
        print("Common shorten complete")
        print("Fixed "+str(inte)+" lines")
        #this should have successfully modified the self.ids list to contain shortened sequence ids.
    def length_check(self, length, verbose):
        #needs to pass in a number... charnum
        toolong = 0
        length = int(length)
        print("trying to shorten to length "+str(length))
        for item in self.ids:
            index = self.ids.index(item)
            linelength = len(item)
            newline= item
            if int(linelength) > int(length):
                toolong +=1
                         
                #change all 12 to 14 if include \n at end of seqids... for now, we are not considering them.
                gi = newline[-12:]
                rest = re.sub("([^#]*)(#)(.*)", "\\1", newline)
                nogi = rest[:-3]
                newl = length-13
                                #12 instead of 12 to leave space for adding a bar.
                newnogi = nogi[:newl]
                if newnogi[-1:] == "|":
                    pass
                else:
                    newnodi = newnogi[:-1]
                newline = newnogi+"|"+gi
                if verbose == True:
                    print ("LENGTHERROR: "+item[:length]+" || "+item[length:])
                    print("Tried to fix: "+newline)
                self.ids[index] = newline
        #end
        print("Length-check complete, "+str(toolong)+" sequences were fixed")

    def weird_AA_check(self, verbose = False):
        lerr = 0
        errletters = []
        for item in self.seqs:
            #if you want to not remove "-" just add it to list of letters.
            listofletters = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'Y']
            newseq = ""
            #open sequences list
            index = self.seqs.index(item)
            anerror = "no"
            for letter in item:
                if letter in listofletters:
                    pass
                elif letter == "\n":
                    pass
                else:
                    if verbose == True:
                        if letter == "-":
                            pass
                        else:
                            print("LETTERERROR: "+letter)
                            anerror = "yes"
                    errletters.append(letter)
                    letter = ""
                    lerr +=1
                newseq = newseq+letter
            if verb == True:
                if anerror == "yes":
                    print(item)
            self.seqs[index] = newseq
        if verbose == True:
            from collections import Counter
            counta = Counter(errletters).most_common()
            print("There were "+str(lerr)+" letter errors as so:")
            print(type(counta))
            for thing in counta:
                print(thing)
        #end
        print("weird aa check done")

    def weird_ID_check(self, verb = False):
        errors = 0
        for item in self.ids:
            index = self.ids.index(item)
            newitem = re.sub("[\[\]]", "", item)
            newitem = re.sub("[:;=,/\+'\.\(\)]", "_", newitem)
            newitem = re.sub(" ", "_", newitem)
            newitem = re.sub("__", "_", newitem)
            if item == newitem:
                pass
            else:
                errors += 1
                if verb == True:
                    print("Replacing:\n"+item+"\n with:\n"+newitem)
            self.ids[index] = newitem
        if verb == True:
            print("there were "+str(errors)+" weird_ID errors")
        print("weird id check done")
        
    def duplicates_check(self, verb = False):
        listoflines = []
        rep = 0
        num = 0
        for line in self.ids:
            index = self.ids.index(line)
            if line in listoflines:
                num+=1
                rep = line+"v"+str(num)
                self.ids[index] = rep
            listoflines.append(line)
        if verb == True:
            print ("there were "+str(num)+" duplicate sequences that were numbered")
        #done
        print("duplicate check done")
    def duplicates_remove(self, verb=False):
        listoflines = []
        to_remove = []
        rep = 0
        num = 0
        for line in self.ids:
            index = self.ids.index(line)
            #if it is a duplicate...
            if line in listoflines:
                num+=1
                to_remove.append(index)
            listoflines.append(line)
        for number in to_remove:
            self.ids.pop(number)
            self.seqs.pop(number)
        if verb == True:
            print ("there were "+str(num)+" duplicate sequences that were removed")
        #done
        print("duplicate ID removal done")
    def index_shorted(self, replace):
        #currently does NOT work w/ accession numbers
        #here replace is depth and/or gi num eg "2 3 gi"
        CTdict = {}
        for line in self.ids:
            if "|gi#" in line:
                taxgi = re.sub("([^#]*)(\|gi#\|?)([0-9]*)(.*)", "\\1~\\3", line)
                tax, gi = taxgi.split("~")
                taxlist = tax.split("|")
                if replace == "gi":
                    CTdict[line] = gi
                if type(replace) is int:
                    CTdict[line] = taxlist[replace-1]
                if type(replace) is str:
                    listreplace = replace.split()
                    newid = ""
                    for item in listreplace:
                        if item == "gi":
                            newid = newid+"|"+gi
                        else:
                            newid = str(newid)+"|"+str(taxlist[int(item)-1])
                    newid = newid
                    CTdict[line] = newid
                    print(newid)
            else:
                tax = re.sub("([^#]*)(\|gi#\|?)([0-9]*)(.*)", "\\1", line)
                taxlist = tax.split("|")
                if replace == "gi":
                    pass
                if type(replace) is int:
                    CTdict[line] = taxlist[replace-1]
                if type(replace) is str:
                    listreplace = replace.split()
                    newid = ""
                    f = 1
                    for item in listreplace:
                        f += 1
                        if item == "gi":
                            newid = newid+"|NA"
                        else:
                            newid = str(newid)+"|"+str(taxlist[int(item)-1])
                            # #SPECIFICALLY FOR CURRENT USE_CASE, REMOVE LATER
                            # if f == 2:
                            #    newid = str(newid)+"|"+str(taxlist[int(item)-1])
                            # if f == 3:
                            #    newid = str(newid)+"|"+str(taxlist[int(item)])
                    newid = newid
                    CTdict[line] = newid
                    print(newid)
        for line in self.ids:
            index = self.ids.index(line)
            newestid = CTdict[line]
            self.ids[index] = newestid
        print("index check done")

    def ten_char(self):
        #something
        #this one should be done in a seperate loop
        CTdict = {}
        iteration = 0
        for line in self.ids:
            iteration +=1
            line = line.strip()

##                  #i have something like
##                  >Methanococcoides_burtonii|gi|909890
##                  #i want
##                  MethBurt00
            GenusSpecies = re.sub("([A-Z][a-z]*)(_)([A-Z]*[a-z]*)(.*)", "\\1~\\3", line)
            try:
                Genus, Species = GenusSpecies.split("~")
                g4 = Genus[:4]
                try:
                    s4 = Species[:4]
                    s3 = Species[:3]
                except:
                    s4 = Species[:2]
                    s3 = Species[:2]
                if iteration < 10:
                    newid = g4+s4.capitalize()+"0"+str(iteration)
                elif iteration > 99:
                    newid = g4+s3.capitalize()+str(iteration)
                else:
                    newid = g4+s4.capitalize()+str(iteration)
            except:
##                      print(GenusSpecies)
                gs8 = GenusSpecies[1:9]
                if iteration < 10:
                    newid = gs8+"0"+str(iteration)
                elif iteration > 99:
                    newid = gs8[:-1]+str(iteration)
                else:
                    newid = gs8+str(iteration)
##                      print(newid)

            CTdict[line] = newid
        for line in self.ids:
            index = self.ids(line)
            newestid = CTdict[line]
            self.ids[index] = newestid
        print("ten char done")
    def mb_version(self):
        #shorten seqids to 94 if not already done.
        self.length_check(94)
        #deal with any duplicates that may have caused
        self.duplicates_check()
        #remove the # and | characters that MrBayes Hates
        for line in self.ids:
            if "#" in nline:
                nline = re.sub("[#]", "", nline)
            if "|" in nline:
                nline = re.sub("\|", "_", nline)
        #tell you what to do
        print("MB version ids created")
        print("You should print this too .fasta format, and then convert to nexus however you want")

    def load_info_swap(self, info_file_in):
        #reads a file of form
        #   originalID
        #   changedID
        #and generates self.ids from that file.
        kid = "no"
        vid = "no"
        CTdict = {}
        with open (info_file_in) as old:
            for line in old:
                #first pass: gets key (original ID)
                #second pass: gets value (new ID)
                #if we have no info, get key
                if kid == "no":
                    key = line.strip()
                    kid = "yes"
                    continue
                elif kid == "yes":
                    #if we have key and value, record.
                    if vid == "yes":
                        CTdict[key]=value
                        vid = "no"
                        kid = "no"
                        continue
                    #if we have key but no value, get value.
                    if vid == "no":
                        value = line.strip()
                        vid = "yes"
            #catch the final pass
            CTdict[key]=value

        if self.original_ids == []:
            for thing in CTdict:
                self.ids.append(thing)
                self.original_ids.append(CTdict[thing])
        else:
            for item in self.original_ids:
                index = self.original_ids.index(item)
                newid = CTdict[item]
                self.ids[index] = newid
        print("original ids:")
        print(self.original_ids)
        print("new ids:")
        print(self.ids)
        #done
        #troubleshooting: do not preform this operation after any that change self.ids. this op must be done first, or in a seperate command.

    def gen_new_fasta(self, new_fasta_name):
        #this should print the changed seqids and changed AA sequences to file.
        newfasta = new_fasta_name
        # print(len(self.original_ids))
        # print(len(self.ids))
        # print(len(self.original_seqs))
        # print(len(self.seqs))
        with open (newfasta, "w") as new:
            for i in range(len(self.original_ids)):
                new.write(">"+self.ids[i].strip()+"\n")
                # print(i)      #
                #unclear if this needs a "\n" after it... check.#TODO
                                #print(self.seqs)
                                #print(type(self.seqs[i]))
                new.write(self.seqs[i]+"\n")
        print("Finished, your new fasta file is located at "+newfasta)
        #done
    def extract(self, list_of_keeps):
        keep_ids = []
        keep_seq = []
        success = 0
        suc_num = len(list_of_keeps)
        for item in list_of_keeps:
            item = item.strip()
            for thing in self.original_ids:
                if thing == item:
                    keep_ids.append(thing)
                    index = self.original_ids.index(item)
                    seq = self.original_seqs[index]
                    keep_seq.append(seq)
                    success += 1
        if suc_num == success:
            print("100% complete extract")
        else:
            print(str(success)+"out of "+str(suc_num)+" sequences extracted")
        self.ids = keep_ids
        self.seqs = keep_seq
        
    def swap_in_newick(self, old_newick_name, new_file_name):
        #this replaces the tip names in a newick file. sometimes works on nexus files too, but I havent extensively tested it.
        newick = old_newick_name
        newnewick = new_file_name
        with open (newick) as old:
            with open (newnewick, "w") as new:
                for line in old:
                    for item in self.original_ids:
                        index = self.original_ids.index(item)
                        line = line.replace(item, self.ids[index])
                    new.write(line)
        print("finished, tip-replaced-newick file at: "+newnewick)
        #done

    def swap_in_nexus(self):
        print ("You didn't implement this yet. try using newick replace, it might work")
        pass
        #something
        #to-do, try nexus replace in the meantime, it should work
    def gen_info(self, info_file_name):
        #writes a file of form
        #   originalID
        #   changedID
        with open(info_file_name, "w") as inf:
            listlength = len(self.original_ids)
            if listlength != len(self.ids):
                print ("List lengths do not match! FATAL ERROR")
                print (self.original_ids)
                print (self.ids)
                raiseSystemExit
            for i in range(listlength):
                inf.write(self.original_ids[i])
                inf.write(self.ids[i]+"\n")
        print("Info file was generated. Named "+info_file_name)
        #done
                
    def write_one_seq_per_file(self):
        geneflist = []
        genenames = []
        for i in range(len(self.ids)):
            with open("Seq" + str(i), "w") as new:
                new.write(">" + self.ids[i]+"\n")
                new.write(self.seqs[i]+"\n")
                name = re.sub("([^\|]*)(\|)(.*)", "\\1", self.ids[i])
                geneflist.append("Seq" + str(i))
                genenames.append(name)
        return geneflist, genenames
        print("one per file generated")

    def number_of_sites(self):
        return len(self.original_seqs[0])

    def shorten(self):
        print("shortening ids...")
        unk = "no"
        normal = 0
        ucount = 0
        for line in self.ids:
            index = self.ids.index(line)

            # this removes words in brackets that aren't Species_name
            # and then changes NCBI's default naming scheme to be
            #>Species_name|#########
            # and makes a list of all gi nums and all
            # duplicates
            # AAH91460.1 Ribosomal protein L3 [Danio rerio]
            if "gi|" in line:
                number = re.sub("(gi)(\|)([0-9]*)(\|)([A-Za-z]*)(\|)(.*)(\[\'?[A-Z]?[a-z]* ?.*\])(.*)", "\\3", line)
                num = number.strip()
                edit1 = re.sub("(gi)(\|)([0-9]*)(\|)([A-Za-z]*)(\|)(.*)(\[\'?[A-Z]?[a-z]* ?.*\])(.*)", "\\8\\2\\1#|", line)
            #get acc number
            else:
                number = re.sub("([^ ]*)(.*)(\[\'?[A-Z]?[a-z]* ?.*\])(.*)", "\\1", line)
                num = number.strip()
                #get edit | AAH91460.1 Ribosomal protein L3 [Danio rerio]
                edit1 = re.sub("([^ ]*)(.*)(\[\'?[A-Z]?[a-z]* ?.*\])(.*)", "\\3|", line)
            if "[" in edit1:
                unk = "no"
                normal += 1
            else:
                unk = "yes"
            edit2 = re.sub("[\[\]]", "", edit1)
            #for now, leave periods in name due to their necessity in acc numbers (????)
            edit3 = re.sub("[:;\.=,/\+'\(\)]", "_", edit2)
            edit4 = re.sub(" ", "_", edit3)
            edit4 = re.sub("__", "_", edit4)
            edit4 = edit4+num
            if unk == "no":
                self.ids[index] = edit4
            else:
                print("Unknown Species in ID:" + line)
        print("shortened: "+str(normal)+" sequence")
        
    def blast2fasta(self, blastlist, ENTREZ=False, num=False):
        # entrez is used to ensure that sequence saved uses correct TAXON, esp. if sequence is a MULTISPECIES entry.
        # entrex should be somethin like "Mycobacterium triplex"
        #take from MakeSPeciesTree.py version if you want a new sequence for each multispecies thing(!)
        # num is how many sequences to write. for species trees, we almost certainly only want one.
        # for converting full downloaded .fastas, we will want all of them (default = False means to do all of them)
        # Converts blast outfmt "6 sseqid stitle sseq" to original lists if
        # entrez = false

        #... now converting outfmt "6 sallseqid salltitles sseq" to sh fasta with selection of proper gi/acc/taxon
        # this should take format " " blast names and replace them with the proper
        # fasta shit
        ernum = 0
        # we open each file in a unique call to blast2fasta. files should be
        # deleted afterwards.
        bf = open(blastlist, 'r')
        error = 0
        end = "no"
        for line in bf:
            if end == "yes":
                break
            # gi|738518257|ref|WP_036466735.1|;gi|620038207|emb|CDO87046.1|   50S
            # ribosomal protein L15 [Mycobacterium triplex]<>50S ribosomal protein L15
            # [Mycobacterium triplex]

            gis = re.sub("(.*)(\t)(.*])(\t)([A-Z-]*)", "\\1", line)
            names = re.sub("(.*)(\t)(.*])(\t)([A-Z-]*)", "\\3", line)
            seq = re.sub("(.*)(\t)(.*])(\t)([A-Z-]*)", "\\5", line)
            # this removes sequences with no Species_name given, so as to avoid errors
            # downstream
            if "\t" in gis:
                error += 1
                print("ERROR in blast parsing: " + line)
                continue
            else:
                gilist = gis.split(";")
                namelist = names.split("<>")
                if ENTREZ is False:
                    index = 0
                else:
                    ENTREZ = ENTREZ.strip("\"")
                    for item in namelist:
                        if ENTREZ in item:
                            index = namelist.index(item)
                try:
                    seqi = gilist[index].strip() + namelist[index].strip()
                    #end = "yes"
                except UnboundLocalError:
                    error += 1
                    print("Name error... might fix")
                    if error == 5:
                        print("Serious ENTREZ error:")
                        print(ENTREZ)
                        print(namelist)
                        print("This gene wasn't found in this taxon, skipping")
                        break
                    continue
                    # goes to next line, abandoning this one
                seqid = re.sub("[ ]", "_", seqi)
                # strips for .fasta format
                seqid = seqid.strip()
                seqid = seqid.strip(">")
                # add the new sequence id to the list.
                self.ids.append(seqid)
                self.original_ids.append(seqid)
                # the new sequence
                slist = []
                count = 0
                newseq = ""
                for letter in seq:
                    if count > 79:
                        count = 0
                        newseq = newseq + ("\n")
                    newseq = newseq + letter
                    count += 1
                self.seqs.append(newseq.strip())
                self.original_seqs.append(newseq.strip())

        print("Blasttofasta id/seq loading complete!")
    def SetTaxID(self):
        self.taxid = []
        for item in self.numbers:
            GItoTAXID = "xmllint --xpath '/GBSet/GBSeq/GBSeq_feature-table/GBFeature/GBFeature_quals/GBQualifier[GBQualifier_name=\"db_xref\"]/GBQualifier_value/text()' \"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id="+item+"&retmode=xml\""
            futuretaxid = subprocess.check_output(GItoTAXID, shell=True)
            taxid = re.sub("(taxon:)([0-9]*)(.*)", "\\2", futuretaxid)
            self.taxid.append(taxid)
            
    def GetTaxonomy(self):
        self.taxonomy = []
        if self.taxid = []:
            print("You need to generate taxids first.. lets try")
            self.SetTaxID()
        for item in self.taxid:
            taxid = number
            ranklist = "superkingdom kingdom phylum class order family"
            ranklist = ranklist.split()
            for r in ranklist:
                TAXIDtoRANKNAME = "xmllint --xpath '/TaxaSet/Taxon/LineageEx/Taxon[Rank=\"" + r + \
                                "\"]/ScientificName/text()'  \"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=" + taxid + "\""
                try:
                    rankname = subprocess.check_output(TAXIDtoRANKNAME, shell=True)
                except:
                    rankname = "NA"
                    rankname = re.sub(" ", "_", rankname)
                taxdict = {}
                taxdict[r]=rankname
                self.taxonomy.append(taxdict)
    def AppendTaxonomy(self):
        for item in self.ids:
            index = self.ids.index(item)
            rankdict = self.taxonomy[index]
            newitem = rankdict["superkingdom"]+"|"+rankdict["kingdom"]+"|"+rankdict["phylum"]+"|"+rankdict["class"]+"|"+rankdict["order"]+"|"+rankdict["family"]+"|"+item
            self.ids[index] = newitem
#TODO:
#add get taxonomy to parser..
#this hasn't been implemented in class fasta, so I am leaving it commented out.. subtrees file might be easily replaced using replace.newick but it might take literally ages... unclear.

# def replace2(replace_file, dict_old_new, charnum, verb):
#    print("Recognized subtrees file, using subtrees varient")
#    outputlist = []
#    rep = 0
#    replist = []
#    newfilename = replace_file.split(".")
#    newfilename = newfilename[0]+str(charnum)+"limit."+newfilename[1]
#    with open(replace_file) as old:
#        if verb == True:
#            print("Opening "+replace_file)
#        with open (newfilename, "w") as new:
#            for line in old:
#                line = line.strip()
#                for item in dict_old_new:
#                    if item[:127] in line:
#                        if item[:127] in replist:
#                            pass
#                        else:
#                            replist.append(item[:127])
#                        rep+=1
# ##                            print(line)
#                        oldline = line
#                        line = line.replace(item[:127], dict_old_new[item])
# ##                        if verb == True:
# ##                            if len(line) <200:
# ##                                print oldline
# ##                                print item
# ##                                print dict_old_new[item]
# ##                                print(line)
# ##                                print("\n")
# ##                            print("\n")
#                new.write(line+"\n")
#            print("finished with "+newfilename+"made "+str(rep)+" replacements of "+str(len(replist))+" differnt patterns")
# ##                print(replist)

#    return newfilename



    # def gen_original_lists(self, fastaname):

    # def load_info_swap(info_file_in):

    # def duplicates_check(verb = False):
    # def weird_ID_check(verb = False):
    # def weird_AA_check(verbose = False):
    # def length_check(length, verbose=False):

    # def manual_shorten():
    # def common_shorten():

    # def mb_version():
    # def index_shorted(replace):
    # def ten_char():

    #    #write stuff
    # def gen_new_fasta(new_fasta_name):
    # def swap_in_nexus():
    # def swap_in_newick(old_newick_name, new_file_name):
    # def gen_info(info_file_name):

# example $ python FISH_2.py /path/to/my/files -fas input.fasta -wf output.fasta -dr
#to remove duplicates from a thing
    
if __name__ == "__main__":

    print("Running in terminal")
    import sys
    import argparse
    import os
    import re
    parser = argparse.ArgumentParser(description="All")
    #necessary bits
    parser.add_argument("directory", nargs='?', default=os.getcwd(), type=str, help="type name of directory to run in where fasta resides, if not pwd")
    parser.add_argument("-fas", "--fasta", action = "store", default = False, help="type the name of your .fasta file")
    #options to load changes from another file
    parser.add_argument("-i", "--infofile", action = "store", default = False, help="Provide an Info File (as generated by this script previously) to pull original and new sequences from")

        
    #options#  to check,fix,edit,etc the seqs or seqids
    # -length
    # -duplicate
    # -weirdaa
    # -weirdID
    parser.add_argument("-sh", "--shorten", action = "store_true", default=False, help="shortens blast (from online) seqIDs")
    parser.add_argument("-b2f", "--blast2fasta", action = "store_true", default=False, help="Blast+ output -> fasta download format BUGGY")
    parser.add_argument("-l", "--length", action = "store", default=False, help="Provide a max length for your sequenceIDs")
    parser.add_argument("-dn", "--duplicates_number", action = "store_true", help="Flag causes identical seqIDs to be numbered 1 2 3 etc to prevent program confusion")
    parser.add_argument("-dr", "--duplicates_remove", action = "store_true", help="Flag causes identical seqIDs to be removed")
    parser.add_argument("-fid", "--fixID", action = "store_true", help="Flag scans SeqIDs and removes weird characters like += etc")
    parser.add_argument("-faa", "--fixAA", action = "store_true", help="Flag scans Sequences and removes non-standard AA characters like X B &")
    #options to shorten specific words
    # -m manual_shorten
    # -c common_shorten

    parser.add_argument("-c", "--common", action = "store_true", help="Flag causes seqIDs to be shortened in a predefined manner, eg bacteriales->bacles ")
    parser.add_argument("-m", "--manual", default = False, action = "store", help="Provide a list of \"original,new\" things to shorten. eg \"Bacteria,Bac Eukaryota,Euk\"")
    #special shortening methods

    parser.add_argument("-t", "--tenchars", action = "store_true", help="Flag turns sequence IDs into ten character strings")
    parser.add_argument("-ba", "--bayes", action = "store_true", help="Flag turns sequences into form that will work as MrBayes input")
    parser.add_argument("-p", "--piece", default = False, action = "store", help="Provide taxonomy-depth, gi, or combo for shortening eg \"1 3 gi\"")

    #writing methods
    parser.add_argument("-wf", "--writefasta", action = "store", default=False, help="Provide name for new fasta file")
    parser.add_argument("-wn", "--writenewick", action = "store",default=False, help="Provide name of newick, name of newfile eg \"example.newick replaced.newick\"")
    parser.add_argument("-wi", "--writeinformation", action = "store", default=False, help="Provide name for this info_file")

    # -fasta
    # -newick replace
    # -info gen (should this always happen?)

    parser.add_argument("-v", "--verbose", action = "store_true", help="prints more information - for debugging mostly. might not be implemented yet")

    args = parser.parse_args()
#workflow: do all the things you want to do to change seqID/seq in one step, save the information and .fasta file.
#then, if desired, use that fasta as base to make ten-char shortened, MBversion, or depth-shortened files, also saving info file so they are reversable.


#actual work flow

#change dir if desiredprint(args.fasta)
    try:
        os.chdir(args.directory)
        if args.verbose == True:
            print("moved to dir: "+args.directory)
    except:
        print ("didn't change dir")
    if args.verbose:
        verb = True
    else:
        verb = False
#originate the fasta class instance

    MyFasta = Fasta("MyFastaName")

    if args.blast2fasta != False:
        MyFasta.blast2fasta(args.fasta)
    else:
        if args.fasta != False:
            MyFasta.gen_original_lists(args.fasta)


    #this should be done in conjunction w / write fasta or replace newick.
    if args.infofile != False:
        MyFasta.load_info_swap(args.infofile)

    #here are the error-fixing calls
    if args.duplicates_number == True:
        MyFasta.duplicates_check(verb)
    if args.duplicates_remove == True:
        MyFasta.duplicates_remove(verb)
    if args.fixID == True:
        MyFasta.weird_ID_check(verb)
    if args.fixAA == True:
        MyFasta.weird_AA_check(verb)
    #shortening calls
    if args.shorten == True:
        MyFasta.shorten()

    if args.common == True:
        MyFasta.common_shorten(verb)
    if args.manual != False:
        MyFasta.manual_shorten(args.manual)

    if args.piece != False:
        MyFasta.index_shorted(args.piece)

    if args.length != False:
        MyFasta.length_check(args.length, verb)

    #these should only be done on their own, not combined w the above. for mrbayes, anything that requires 10 characters.
    if args.bayes == True:
        MyFasta.mb_version()
    if args.tenchars == True:
        MyFasta.ten_char()

    #write stuff
    if args.writefasta != False:
        MyFasta.gen_new_fasta(args.writefasta)
    if args.writenewick != False:
        old, new = args.writenewick.split()
        MyFasta.swap_in_newick(old, new)
    if args.writeinformation != False:
        MyFasta.gen_info(args.writeinformation)
    print("All things finished, exiting...")


#TODO
#detailed information on how to use
#test everything
#????
# FISH FASTA ID SWAPPING HELPER

    #### this is becoming a dedicated tip-name-editing package.

    #requires that tips are in format given by FEAST's shorten or shorten-keep-info
    #things it can do:
    #   1. shorten too long seqids using common shortening phrases, or by removing info from the species-name (usually catches strain info)
    #   2. remove weird characters from seqIDS
    #   3. remove weird characters from AA sequences
    #   do this on 1. fasta files 2. nexus files (maybe? unclear) 3. newick files (maybe? unclear)
