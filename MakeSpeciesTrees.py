# #!/usr/bin/python

# last edit abigailc@Actaeon on october 21 2016


#major speed increase... first run should take 2-4 hours, subsequent runs from same directory ~1 hour.
#usage
#assume you want to use species that are found in the single_gene_fasta file example.fasta, and genes as found in Euk_Ribo.fasta
#always use -r flag for now, it is faster.
#if you want to use species from a list file (one species per line) use -s SPECIESFILE instead of -f FASTAFILE
# $python MakesSpeciesTrees.py -f example.fasta -g Euk_Ribo.fasta -p MyTree1 -r





# #to use this script, you need to do several things.
# 1. set the ssh path so you will be able to run things on the cluster.
# 2. set the cluster_head path so files will be transferred to/from your partition (not mine!)
# 3. run from the command line
# a. cd to directory you are keeping this script.
# b. identify the genes you want your species tree to be built from, and create a .fasta file containing them. For ease of understanding some produced files, I am using the sequence id convention
# >GeneName|other_sequence_information_blah_blah
# You can use the files I have here (Archaeal_Ribo.fasta, Bacterial_Ribo.fasta, Eukaryal_Ribo.fasta) or create your own query file. Please provide the full path to the query file, or keep it in the same folder as the directory you indicate when calling MakeSpeciesTree.
# indicate the gene-queries to use with the -genes tag.
# c.  identify the species that should be included in your species tree. They should either be in a shorten-formatted .fasta file, or in a plaintext file with one species name on each line of the file eg.
# Cat
# Dog
# Rat
# indicate the species file with -species or -fasta if you are giving a shortened fasta instead.
# d. if you want your output to be named something reasonable (like SOD_cyanos_blah) give it a prefix using -p
# e. run something like $python MakesSpeciesTrees.py -s species_file.txt -g bacterial_ribo_gene_seqs.fasta -p MyProject_SpeciesTree
# 4. wait. it might take a long time to run all of the blast searches. then the program will send them to the cluster to align (per gene), wait until they are done checking every... 5 minutes, concatenate the genes for each species, and then send the concatenated file to the cluster to run raxml, and download when it is complete. This might take some time, do not close your terminal or turn your computer off while it is running. subsequent runs with the same gene/will not need to re-do the blast search, so will proceed more quickly.


##this is also usable as a module, when calling from within another program it is probably best to pass a list of trees to make like so [ [spec1, gene1, prefix1], [spec2, gene2, prefix2] ] to get_multiple_concat_alignments()
#& then also run the run_raxml_on_cluster thingy with the output muscle_aligned_list_of_files.

# idea to save time; implement before final run through of giant project
# make sure to save each hit in a current-project file, but ALSO append it to an OVERALL file, and
# when we read a new taxa/gene pair, consider searching the OVERALL file
# before running a new blast DONE

# for expediency in problem-solving when sharing script, just keeping bits of Fasta module that are needed contained here.


######### PERSONAL_SETTINGS #########
ssh_inst = "ssh -l abigailc -i ~/.ssh/id_rsa eofe5.mit.edu"
clus_head = "abigailc@eofe5.mit.edu:/home/abigailc/"
Path_Blast = "/Users/abigailc/blast/"

#for local blast use:
#i have nr database downloaded (via blast+ tool '$update_blastdb.pl nr' 'ls *.gz|xargs -n1 tar -xzvf' 'rm *.gz.*') and the download location added to my .profile (export BLASTDB=Users/abigailc/blast)

##########

import sys
import argparse
import os
import re
import time

#imports
def check_directory_existance(prefix, ssh_inst):
    import os
    print(type(ssh_inst))
    print(type(prefix))
    print(prefix)
    os.system(ssh_inst+" \' mkdir Species_Trees;cd Species_Trees;mkdir "+prefix+"\'")
    
def remove_slurm_files(ssh_inst, prefix, pattern):
     os.system(ssh_inst+" \' cd Species_Trees;cd "+prefix+"; rm "+pattern+"\'")

    
class Fasta:
    def __init__(self, name="whatever"):
        # all ids should be stripped and have ">" removed for reasons.
        # for now, sequences do not have any stripping applied
        self.name = name
        self.ids = []
        self.original_ids = []
        self.original_seqs = []
        self.seqs = []

    def gen_original_lists(self, fastaname):
        try:
            with open(fastaname) as fastafile:
                for line in fastafile:
                    if "\n" == line:
                        pass
                    if ">" in line:
                        # write the previous AA seq
                        try:
                            AAseq = AAseq.strip()
                            self.seqs.append(AAseq)
                            self.original_seqs.append(AAseq)
                        except:
                            pass
                            # initialize a new AAseq
                        AAseq = ""
                        # format the seqID
                        newline = line.strip()
                        newline = newline.strip(">")
                        # write the seqID
                        self.ids.append(newline)
                        self.original_ids.append(newline)
                    else:
                        AAseq = AAseq + line
                AAseq=AAseq.strip()
                # catch the last AAseq pass
                self.seqs.append(AAseq)
                self.original_seqs.append(AAseq)
            print("Initial sequence and ID lists created. Contains " +
                  str(len(self.ids)) + " sequences")
        except UnboundLocalError:
            print("probably this file :" + fastaname +
                  " has nothing in it. skipping.")
            pass
        except IOError:
            print("no file named: " + fastaname +
                  " exists... creating a blank file")
            with open(fastaname, "w") as new:
                pass
            print("hopefully you intended that!")

    def gen_species_blast(self):
        self.species_names = []
        #add each thing - the bit between start and | in each id.
        for item in self.ids:
            nam = re.sub("([^\|]*)(\|)(.*)", "\\1", item)
            nam.strip()
            self.species_names.append(nam)
            
    def gen_raw_blast_lists(self, fastaname):
        bf = open(fastaname, 'r')
        self.species_names = []
        for line in bf:
            gis = re.sub("(.*)(\t)(.*)(\t)([A-Z-]*)", "\\1", line)
            names = re.sub("(.*)(\t)(.*)(\t)([A-Z-]*)", "\\3", line)
            seq = re.sub("(.*)(\t)(.*)(\t)([A-Z-]*)", "\\5", line)
            if "\t" in gis:
                print("ERROR in blast parsing: " + line)
                continue
            else:
                gilist = gis.split(";")
                namelist = names.split("<>")
            for name in namelist:
                index = namelist.index(name)
                gi = gilist[index]
                gi = re.sub("(gi\|)([0-9]*)(.*)", "\\2", gi)
                seqid = re.sub("[ ]", "_", name)
                seqid = seqid.strip()
                #test that the identified species is the same as you want it to be... if not, skip and try again.
                #because some dummy likes to name "gene like in [arabidopsis thaliana] [mus musculus]"
                species_as_found_in_gsl = re.sub("([^\[]*)(.*)", "\\2", seqid)
                species_as_found_in_gsl = re.sub("[\[\]]", "",  species_as_found_in_gsl)
                species_name = re.sub("[_]", " ",  species_as_found_in_gsl)
                if species_name == "":
                    continue
                #this will return a species name of form "Danio rerio"
                #lists
                species_name.strip()
                if species_name in self.species_names:
                    pass
                else:
                    newID = species_as_found_in_gsl+"|gi#|"+gi
                    self.ids.append(newID)
                    self.original_ids.append(newID)
                    self.seqs.append(seq)
                    self.original_seqs.append(seq)
                    self.species_names.append(species_name)
    
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

    def number_of_sites(self, num = 0):
        testseq = self.original_seqs[num]
        testseq = re.sub("\n", "", testseq)
        #print(testseq)
        #print (len(self.original_seqs[num]))
        #print (len(testseq))
        return len(testseq)

    def shorten(self):
        unk = "no"
        normal = 0
        ucount = 0
        for line in self.ids:
            index = self.ids.index(line)

            # this removes words in brackets that aren't Species_name
            # and then changes NCBI's default naming scheme to be
            #>Species_name|gi#|#########
            # and makes a list of all gi nums and all
            # duplicates
            number = re.sub(
                "(gi)(\|)([0-9]*)(\|)([A-Za-z]*)(\|)(.*)(\[\'?[A-Z]?[a-z]* ?.*\])(.*)", "\\3", line)
            num = number.strip()
            edit1 = re.sub(
                "(gi)(\|)([0-9]*)(\|)([A-Za-z]*)(\|)(.*)(\[\'?[A-Z]?[a-z]* ?.*\])(.*)", "\\8\\2\\1#|\\3", line)
            if "[" in edit1:
                unk = "no"
                normal += 1
            edit2 = re.sub("[\[\]]", "", edit1)
            edit3 = re.sub("[:;=,/\+'\.\(\)]", "_", edit2)
            edit4 = re.sub(" ", "_", edit3)
            edit4 = re.sub("__", "_", edit4)
            if unk == "no":
                self.ids[index] = edit4
            else:
                print("Unknown Species in ID:" + line)

    def blast2fasta(self, blastlist, ENTREZ=False, num=False):
        #returns "bad" if bad. else returns True
        # entrez is used to ensure that sequence saved uses correct TAXON, esp. if sequence is a MULTISPECIES entry.
        # entrex should be somethin like "Mycobacterium triplex"
        # num is how many sequences to write. for species trees, we almost certainly only want one.
        # for converting full downloaded .fastas, we will want all of them (default = False means to do all of them)
        # Converts blast outfmt "6 sseqid stitle sseq" to original lists if
        # entrez = false

        #... now converting outfmt "6 sallseqid salltitles sseq" to sh fasta with selection of proper gi/acc/taxon
        # this should take format " " blast names and replace them with the proper
        # fasta shit
    
        # we open each file in a unique call to blast2fasta. files should be
        # deleted afterwards.
        #print(blastlist)
        #print(os.getcwd())
        bf = open(blastlist, 'r')
        bf2 = open(blastlist, 'r')
        #print(bf)
        error = 0
        end = "no"
        bad = "no"
        length_bf = 0
        for line in bf2:
            length_bf+=1
        errnum = length_bf
        #print("Maxerror: "+str(length_bf))
        run = 0
        #this should be 5 if called from current setup, but might be less if less hits are found.
        for line in bf:
            run +=1
            if end == "yes":
                break
            # gi|738518257|ref|WP_036466735.1|;gi|620038207|emb|CDO87046.1|   50S
            # ribosomal protein L15 [Mycobacterium triplex]<>50S ribosomal protein L15
            # [Mycobacterium triplex]
            #i just removed a ] from group3/... used to be (.*]) but occasionally species name isn't given at end of thing causes errors.
            gis = re.sub("(.*)(\t)(.*)(\t)([A-Z-]*)", "\\1", line)
            names = re.sub("(.*)(\t)(.*)(\t)([A-Z-]*)", "\\3", line)
            seq = re.sub("(.*)(\t)(.*)(\t)([A-Z-]*)", "\\5", line)
            # this removes sequences with no Species_name given, so as to avoid errors
            # downstream
            if "\t" in gis:
                error += 1
                #print(error)
                print("ERROR in blast parsing: " + line)
                continue
            #check if entrez is in the thing
            else:
                gilist = gis.split(";")
                namelist = names.split("<>")
                if ENTREZ is False:
                    index = 0
                else:
                    ENTREZ = ENTREZ.strip("\"")
                    found = "no"
                    for item in namelist:
                        if ENTREZ+"]" in item:
                            index = namelist.index(item)
                            found = "yes"
                    if found == "no":
                        #entrez not found in this thing
                        error += 1
                        print("Name error... might fix")
                        print(namelist)
                        if error == errnum:
                            print("Serious ENTREZ error:")
                            print(ENTREZ)
                            print("This gene wasn't found in this taxon, skipping")
                            bad = "yes"
                            break
                        continue
                        
                        
                try:
                    seqi = gilist[index].strip() + namelist[index].strip()
                except UnboundLocalError:
                    #print(error)
                    print("Name error... might fix")
                    if error == errnum:
                        print("Serious ENTREZ error:")
                        print(ENTREZ)
                        print(namelist)
                        print("This gene wasn't found in this taxon, skipping")
                        bad = "yes"
                        break
                    continue
                except:
                    print(index)
                    error += 1
                    print("unknown index error")
                    if error == errnum:
                        break
                    continue
                #print("passed unbound error")
                    # goes to next line, abandoning this one
                
                # strips for .fasta format
                seqid = re.sub("[ ]", "_", seqi)
                seqid = seqid.strip()
                seqid = seqid.strip(">")
                #test that the identified species is the same as you want it to be... if not, skip and try again.
                #because some dummy likes to name "gene like in [arabidopsis thaliana] [mus musculus]"
                species_as_found_in_gsl = re.sub("([^\[]*)(.*)", "\\2", seqid)
                species_as_found_in_gsl = re.sub("[\[\]]", "",  species_as_found_in_gsl)
                species_as_found_in_gsl = re.sub("[_]", " ",  species_as_found_in_gsl)
                if ENTREZ is False:
                    pass
                else:
                    if ENTREZ == species_as_found_in_gsl:
                        pass
                    else:
                        #print(ENTREZ)
                        #print(species_as_found_in_gsl)
                        #print(seqid)
                        print("Odd, species identified did not match Entrez. Trying next seq, might fix.")
                        error += 1
                        #print(error)
                        end = "no"
                        if error == errnum:
                            print("Serious ENTREZ error:")
                            print(ENTREZ)
                            print(namelist)
                            print("This gene not found in this taxon, skipping")
                            print("Trying to label as non-existant to prevent re-searching later... ")
                            bad = "yes"
                            break
                        continue
               
                #print("passed entrez error")
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
                end = "yes"
                #print(newseq)
                #print(seqid)
                self.seqs.append(newseq.strip())
                self.original_seqs.append(newseq.strip())
        #print(error)
        #print(run)
        #how do you manage to get to here, without either being set "Bad" or set good...
        if bad == "yes":
            #print("Has been defined as bad")
            entrez = ENTREZ.replace(" ", "_")
            self.ids.append(entrez+"|gi#|000000000")
            self.original_ids.append(entrez+"|gi#|000000000")
            self.original_seqs.append("----")
            self.seqs.append("----")
            return "bad"
        else:
            #print("returning now")
            return True

    def load_info_swap(self, info_file_in):
        # reads a file of form
        #   originalID
        #   changedID
        # and generates self.ids from that file.
        kid = "no"
        vid = "no"
        CTdict = {}
        with open(infofile) as old:
            for line in old:
                # first pass: gets key (original ID)
                # second pass: gets value (new ID)
                # if we have no info, get key
                if kid == "no":
                    key = line.strip()
                    kid = "yes"
                    continue
                elif kid == "yes":
                    # if we have key and value, record.
                    if vid == "yes":
                        CTdict[key] = value
                        vid = "no"
                        kid = "no"
                        continue
                    # if we have key but no value, get value.
                    if vid == "no":
                        value = line.strip()
                        vid = "yes"
            # catch the final pass
            CTdict[key] = value
        for item in self.original_ids:
            index = self.original_ids.index(item)
            newid = CTdict[item]
            self.ids[index] = newestid
        # done
        # troubleshooting: do not preform this operation after any that change
        # self.ids. this op must be done first, or in a seperate command.

    def gen_new_fasta(self, new_fasta_name):
        # this should print the changed seqids and changed AA sequences to
        # file.
        newfasta = new_fasta_name
        # print(len(self.original_ids))
        # print(len(self.ids))
        # print(len(self.original_seqs))
        # print(len(self.seqs))
        with open(newfasta, "w") as new:
            for i in range(len(self.original_ids)):
                new.write(">" + self.ids[i].strip() + "\n")
                # print(i)  #
                # unclear if this needs a "\n" after it... check.#TODO
                new.write(self.seqs[i].strip()+"\n")
        print("Finished, your new fasta file is located at :" + newfasta)
        # done
        return newfasta

    def swap_in_newick(self, old_newick_name, new_file_name):
        # this replaces the tip names in a newick file. sometimes works on nexus
        # files too, but I havent extensively tested it.
        newick = old_newick_name
        newnewick = new_file_name
        with open(newick) as old:
            with open(newnewick, "w") as new:
                for line in old:
                    for item in self.original_ids:
                        index = self.original_ids.index(item)
                        line = line.replace(item, self.ids[index])
                        new.write(line)
        print("finished, tip-replaced-newick file at: " + newnewick)
        # done

    def swap_in_nexus(self):
        print (
            "You didn't implement this yet. try using newick replace, it might work")
        pass
        # something
        # to-do, try nexus replace in the meantime, it should work

    def gen_info(self, info_file_name):
        # writes a file of form
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
                inf.write(self.ids[i] + "\n")
        print("Info file was generated. Named " + info_file_name)
        # done
# this ought to be deleted - but needs to be removed from the main  code first.

    def pick_name(self, QUERY, n):
        pass
        #

    def gen_species_lists(self):
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
        return speclist

#########
########
# end of fasta functions
#########


####begin ConCat Functions########

# this will correlated IDs based on Genus_species across datasets.
def correlate_ids(list_of_id_lists):
    # so now we have a list of species from each Fasta, and we need to correlate them. should be it's own matching function?
    # go though first list, save name (index list1 ,index list2, indexlist3), and add name to "used" list
    # then go through second list, save same thing but skip any that are already in "used" list.
    # continue etc.
    numlists = len(list_of_id_lists)
    output_list = []
    used = []
    for item in list_of_id_lists:
        # list1
        for ids in item:

            # id1 of list1
            name = ids
            if name in used:
                pass
            else:
                used.append(name)
                id_index_list = []
                id_index_list.append(ids)
                # check the index of that id in each list and append to "id_index_list"
                # if the id is not in that list, should append "NA"
                for eachlist in list_of_id_lists:
                    try:
                        index = eachlist.index(ids)
                    except ValueError:
                        index = "NA"
                    id_index_list.append(index)
                # add the result of scanning that id to overall output list.
                output_list.append(id_index_list)
    # output list looks like:
    # outputlist = [ ["Cat",1,2,3,4] , ["Dog", 2,1,13,14] ]
    print("Concat")
    print(output_list)
    return output_list


# this should work but hasn't been tested. requires stuff above.
def Concatenate(listoffastafiles, new_cc_fasta_name, prefix, strain_limit=False):
        # for each thing in listoffastas, create a fasta object in a list.
    fasta_class_list = []
    list_species_lists = []
    # this will create an instance of Fasta from each fasta file
    for i in listoffastafiles:
        fasta_class_list.append(Fasta(i))
    # this will create the original lists (seqID and sequence) for each Fasta object
    # this will create a list of the "Genus_species" from each seq in each
    # Fasta object (pass true/false of strain control if need be)
    for f in fasta_class_list:
        f.gen_original_lists(f.name)
        species_list = f.gen_species_lists()
        list_species_lists.append(species_list)
    # so now we have a list of species from each Fasta, and we need to correlate them. should be it's own matching function?
    # go though first list, save name (index list1 ,index list2, indexlist3), and add name to "used" list
    # then go through second list, save same thing but skip any that are already in "used" list.
    # continue etc.
    # check that things i already did for greg... swapping placement in file
    # was similar.. where did i save that.
    indexed_ids_list = correlate_ids(list_species_lists)
    # indexed_ids_list in form [ ["Cat", 1,2,3],["dog",2,1,2] ]
    # this part will create a new, concatenated .fasta file
    # requires indexed_ids_list, fasta_class_list, function that returns
    # number of sites. self.number_of_sites
    #debugging:
    lenlist_final = []
     #ensure that for each fas in fasta class list, all sequences are of the same length.
    for fas in fasta_class_list:
        ndash1 = fas.number_of_sites
    with open(new_cc_fasta_name, "w") as new:
        for item in indexed_ids_list:
            lenlist = []
            # do i want to print species_name to file, or full
            # taxonomy? for now, just species_name is gunna happen,
            # but easy to swap i think
            new.write(">" + item[0].strip()+"\n")
            fas_num = 0
            allseq = ""
            #ensure that for each fas in fasta class list, all sequences are of the same length.
            for fas in fasta_class_list:
                fas_num += 1
                # fas_num keeps track of what number fasta we are on, which
                # correlates to the index of the index in indexed_ids_list
                search_index = item[fas_num]
                # search_index will be something like "22"

                # if search_index is NA, generate a str of "-" that is n
                # characters long, where n is the return from
                # fas.number_of_sites
                if search_index == "NA":
                    ndash = fas.number_of_sites()
                    retreived_seq = ""
                    for i in range(int(ndash)):
                        retreived_seq = retreived_seq + ("-")
                else:
                    retreived_seq = fas.seqs[search_index]
                    # retreived_seq wil be something like "the 22nd sequence in
                    # object Fas's sequence list... " or "BLAHSEQUENCEDATA"
                retreived_seq = re.sub("\n", "", retreived_seq)
                lenlist.append(len(retreived_seq))
                count = 0
                allseq = allseq + retreived_seq
            #print(lenlist)
            #print(len(allseq))
            lenlist_final.append(len(allseq))
            newseq = ""
            for letter in allseq:
                if count > 79:
                    count = 0
                    newseq = newseq + ("\n")
                newseq = newseq + letter
                count += 1
                
            new.write(newseq.strip()+"\n")
                # there might be too many line-breaks in this fashion. if so,
                # concatenated all retreived_seqs and then only call write()
                # once.
        for length in lenlist_final:
            if length == lenlist_final[0]:
                pass
            else:
                print("ERROR your concat sequences are not all of the same length something has gone horribly wrong! aborting.")
                raise SystemExit
        print("Should be finished generating new concatenated fasta at: " + new_cc_fasta_name)
    print("done w cc gen!")
    return new_cc_fasta_name

#########################MUSCLE STUFF####################


def muscle_align_on_cluster(end_file_list, prefix):
    #this creates dir you will use on the cluster.
    aligned_list = []
    for item in end_file_list:
        aligned_list.append(item+"_Muscle.fasta")
    check_directory_existance(prefix, ssh_inst)
    clus_path = "/Species_Trees"
    a = gen_muscle_script(prefix+"_Sc.sh", "~"+clus_path+"/"+prefix+"/"+prefix+"_Corr.txt", str(len(end_file_list)), prefix+"job")
    b = gen_correlate_file(end_file_list, prefix+"_Corr.txt")
    end_file_list.append(a)
    end_file_list.append(b)
    direct = os.getcwd()


    move_to_cluster(end_file_list, clus_path+"/"+prefix)
    print("everything should be generated and on the cluster")
    os.system(ssh_inst+" 'cd ~/Species_Trees/"+prefix+";echo $PWD;sbatch "+a+"'")
    finished = "start"
    #to see if the run is complete, see if each new file has been generated. check every 5 minutes for muscle.
    movehome = []
    for i in aligned_list:
        movehome.append(i)
    while finished is not True:
        for filename in movehome:
            os.system("scp "+clus_head[:-1]+clus_path+"/"+filename[2:]+" "+direct+"/"+prefix)
        for item in aligned_list:
            #see if it got moved home.
            exists = os.path.isfile(item)
            if exists is True:
                if item in movehome:
                    movehome.remove(item)
                finished = "yes"
            else:
                finished = False
                print("waiting 120 seconds and trying again")
                print("didn't find :"+item)
                break
        if finished == "yes":
            print("Should be done!")
            finished = True
        else:
            #wait five minutes and then try again.
            time.sleep(120)
            finished = "yes"
    #now, concatenate the output file. when doing large batches, might want to move each to new folder.... but maybe not for now.
    print("Beginning concatenation")
    c = Concatenate(aligned_list, prefix+"_CC.fasta", prefix)
    print("Your file is located at :"+c)
    return c

    #now... run_raxml_on_cluster.... take all the c's and run them.another layer of abstraction.
   # gen_raxml_script()
    
def gen_muscle_script(scriptfile, indexname, n, Jobname):
    #currently assuming you are running in the dir that files are in and should be returned to.
    # direct = os.getcwd()
    # host = socket.gethostname()
    # user = getpass.getuser()
    # addr = user+"@"+host+":"+direct
    #figure out how many need to be run. n = len(listoffilesinindex). n
    #figure out a name - scriptfile
    #figure out path to index file. indexname
    #thats it. just print the script. return its filename, which will need to be added to list of things to be moved to the cluster.
    addr = "PLACEHOLDER"
    

##example script
    a =  """#!/bin/bash                                                                                             
#SBATCH -p sched_mit_g4nier                                                                             
#SBATCH -t 0-00:50:00                                                                                   
#SBATCH -J Mus"""+Jobname+"""                                                                                         
##SBATCH -o Mus"""+Jobname+""".out
#SBATCH --array=1-"""+n+"""

. /etc/profile.d/modules.sh
module add engaging/muscle/3.8.31
##gets my array id, which is needed for use below. this will be, i think, a number like 1,2,3 etc
MY_ARRAY_ID=$SLURM_ARRAY_TASK_ID
echo $MY_ARRAY_ID

## given an index file formatted                                                                        
## <index> <filename>                                                                                   
## produce the filename for given index                                                                 
THE_INDEX="""+indexname+"""
THE_INPUT_FILE=$( cat $THE_INDEX | grep "^$MY_ARRAY_ID " | awk '{print $2}' )
echo $THE_INPUT_FILE
ENDING=_Muscle.fasta
echo $THE_INPUT_FILE$ENDING

muscle -in $THE_INPUT_FILE -out $THE_INPUT_FILE$ENDING

##scp $THE_INPUT_FILE$ENDING """+addr+"""

exit"""
    with open(scriptfile, "w") as script:
        script.write(a)
    return scriptfile




def gen_raxml_script(scriptfile, indexname, n, Jobname):
    #is jobname == prefix?
    #add second /prefix/ to indexname.
    #currently assuming you are running in the dir that files are in and should be returned to.
    # direct = os.getcwd()
    # host = socket.gethostname()
    # user = getpass.getuser()
    # addr = user+"@"+host+":"+direct
    #figure out how many need to be run. n = len(listoffilesinindex). n
    #figure out a name - scriptfile
    #figure out path to index file. indexname
    #thats it. just print the script. return its filename, which will need to be added to list of things to be moved to the cluster.
    

##example script
    a =  """#!/bin/bash                                                                                             
#SBATCH -p sched_mit_g4nier                                                                             
#SBATCH -t 0-10:00:00    
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
# #SBATCH --exclusive                                                                                   
#SBATCH -J RAX"""+Jobname+"""   
#SBATCH -o RAX"""+Jobname+""".out                                                                                         
#SBATCH --array=1-"""+n+"""

. /etc/profile.d/modules.sh
module add engaging/RAxML/8.2.9
##gets my array id, which is needed for use below. this will be, i think, a number like 1,2,3 etc
MY_ARRAY_ID=$SLURM_ARRAY_TASK_ID
echo $MY_ARRAY_ID

## given an index file formatted                                                                        
## <index> <filename>                                                                                   
## produce the filename for given index                                                                 
THE_INDEX="""+indexname+"""
THE_INPUT_FILE=$( cat $THE_INDEX | grep "^$MY_ARRAY_ID " | awk '{print $2}' )
echo $THE_INPUT_FILE
ENDING=_Muscle.fasta
echo $THE_INPUT_FILE$ENDING

NEW=${THE_INPUT_FILE%%.*}
echo $NEW
  
raxmlHPC-PTHREADS-AVX -T 20 -f a -m PROTGAMMALGF -p 12345 -x 12345 -#100 -n $NEW -s $THE_INPUT_FILE         

exit"""
    with open(scriptfile, "w") as script:
        script.write(a)
    return scriptfile



def gen_correlate_file(list_of_input_files, corr_file):
    #this should be in form
    #1 name
    #2 name
    #3 name
    #requires 1: list of files 2. name for corr_file.
    i = 0
    with open(corr_file, "w") as corr:
        for item in list_of_input_files:
            if "/" in item:
                itemlist = item.split("/")
                item = itemlist[-1]
            i += 1
            #make sure the \n spacing works correctly.
            corr.write(str(i)+" "+item+"\n")
    return corr_file


def move_to_cluster(list_of_files, clus_path):
    #requires os.system
    #requires scp(?)
    #do the thing
    for item in list_of_files:
        os.system("scp "+item+" "+clus_head+clus_path)
    print("Finished moving files to cluster in place:"+clus_path)

#
#################################
# this is main fxn, should be at bottom
###############################



def MakeSpeciesTreeCC(species_list, gene_sequences_file, local, prefix=False):


    # species_names_file should look like this:
    #    "homo sapiens"
    #    "boops boops"
    #    "scopthalamus rhombus"

    # OR LIKE
    # 12334
    # 22344
    # 45645
    # where the numbers are TAXIDs

    # gene_sequences_file should look like this:
    # >Gene_MAT|Species_Cat
    # SUCHAMINOACIDSWOWVERYSEQUENCE
    # >Gene_FOO|Species_Cat
    # SUCHAMINOACIDSWOWVERYSEQUENCE
    # >Gene_BOX|Species_Dog
    # SUCHAMINOACIDSWOWVERYSEQUENCE

    # notably, your output files will be labeled the thing between ">" and "|"
    # so in this case, they'll be called "Gene_MAT.fasta Gene_FOO.fasta Gene_BOX.fasta" etc
    # if prefix specified, "Cyano_Gene_MAT.fasta" etc.
    
    #make directory
    os.system("mkdir Stored_Blasts; mkdir "+prefix)
    #put species list in dir for future reference
    write_species_list(species_list, prefix)
    gene_fasta_list = []
    end_file_list = []
    # initialize fasta instance for genes
    GeneFasta = Fasta()
    if os.path.isfile(gene_sequences_file) is False:
        print("Error, cannot find specified gene_sequences_file:"+gene_sequences_file)
        raise SystemExit
    GeneFasta.gen_original_lists(gene_sequences_file)
    # write one fasta file per sequence
    gene_fasta_list, gene_name_list = GeneFasta.write_one_seq_per_file()
    # gene_fasta_list contins the names of each auto-generated .fasta file, each of which contains exactly one sequence... each representing a gene you want in your gene tree.
    # overallfastas... to cut down on search time if repeating for similar genes and species.
    # Make sure this doesn't cause errors by giving each gene a distinct name (EukS4)
    list_overall_gene_fastas = []
    list_raw_gene_fastas = []
    gene_overall_fasta_names = []
    for b in range(len(gene_name_list)):
        gene_name = gene_name_list[b].strip()
        gene_overall_fasta = "./Stored_Blasts/"+gene_name + "_All.fasta"
        gene_overall_fasta_names.append(gene_overall_fasta)
        list_overall_gene_fastas.append(Fasta(gene_overall_fasta))
    for c in list_overall_gene_fastas:
        c.gen_original_lists(c.name)
    print("beginning blast process")
    n = 0
    i = 0
    output_list = []
    ##########if local is true
    # first, initialize (or check existance of) files containing all genes we care about.
    #similar to listoverallgenefastas but instead of containing things we've used, they contain raw blast data.
    #1. check existance of blast-result-files
    #2. if no, do the blast
    #3. now, check if each query/species pair are in the overall files. yes - grab them. no - trawl through the raw blast data.
    if local is True:
        print("Beginning local blast!")
        #make dir
        os.system("mkdir Raw_Blasts")
        #set up lists
        raw_fasta_list = []
        gene_raw_fasta_names = []
        #get names of raw_blasts (all)
        for d in range(len(gene_name_list)):
            gene_name = gene_name_list[d].strip()
            gene_raw_fasta = "./Raw_Blasts/"+gene_name + "_raw.fasta"
            gene_raw_fasta_names.append(gene_raw_fasta)
            list_raw_gene_fastas.append(Fasta(gene_raw_fasta))
        #determine which raw_blasts still need to be run:
        run_blast = []
        query_runblast = []
        convert_blast=[]
        timer = 0
        print(gene_raw_fasta_names)
        for f in gene_raw_fasta_names:
            #check if the blast has already been run
            g = os.path.isfile(f+"_blast.txt")
            #if not, add it to the list of blasts to run
            if g == False:
               run_blast.append(f)
               query_runblast.append(gene_fasta_list[timer])
            #if blast exists, check that conversion has happened
            else:
                j = os.path.isfile(f+"_raw_blast")
                #if not, add to list of needs to be converted
                if j == False:
                    #add needs to be converted, not run_blast
                    convert_blast.append(f)
            timer+=1
                #we need to send which query to run blast as well as which filename? else blast_)query will fail
                #timer will set this correctly (probably)
        h = []
        timer2 = 0
       
        #run blast is a list of 30 29's right now. why.
        print("Need to run blast on "+str(len(run_blast))+" items")
        print(run_blast)
        time1 = time.clock()
        everything = ""
        #we will run 4 at a time for efficiency's sake, tracked by variable iteration
        iteration = 0
        for item in run_blast:
            QUERY = query_runblast[timer2]
            timer2+=1
            h.append(Fasta(item+"_raw_blast"))
            this_blast_query = "blastp -query " + QUERY + " -remote -db nr -out " + item+"_blast.txt -max_target_seqs 20000 -evalue 1e-4 -outfmt \"6 sallseqid salltitles sseq\""
            this_blast_query = this_blast_query+" & sleep 5 ; "
            iteration += 1
            #if we have 4 loaded: go
            if iteration > 3:
            #reset iteration
                iteration = 0
                everything = everything+this_blast_query
                everything = everything+"wait"
                print("Beginning 4 blast searches")
                print(everything)
                c = os.system(everything)
                everything = ""
            else:
                everything = everything+this_blast_query
        #catch the final blast(s)
        everything = everything+" wait"
        print("Beginning final blast searches")
        print(everything)
        c = os.system(everything)
        print("Done with all blast searches")
        time2 = time.clock()
        print("That took time = ")
        print(time2-time1)
            #run the blast
##ADD THE ASSOCIATE QUERY FILE TO THE FASTA??? or otherwise how do track which query matches
            ###
            #5 for test, should be ~20,000

            #to run locally:
        #     blast_query = "blastp -query " + QUERY + " -db "+Path_Blast+"nr -out " + item+"_blast.txt -max_target_seqs 5 -evalue 1e-4 -outfmt \"6 sallseqid salltitles sseq\""

            #to run remotely (faster, unless local db is on a supercomputer maybe)
            # blast_query = "blastp -query " + QUERY + " -remote -db nr -out " + item+"_blast.txt -max_target_seqs 20000 -evalue 1e-4 -outfmt \"6 sallseqid salltitles sseq\""
            # print(blast_query)
            # os.system(blast_query)
            # #gen these Fasta objects
            # h.append(Fasta(item+"_raw_blast"))
            
        for newthing in convert_blast:
            h.append(Fasta(newthing+"_raw_blast"))
        print("Need to do conversion on "+str(len(h))+" items")
        print(h)
        for i in h:            
            #gen the lists including species list from the raw blast data
            print("this")
            i.gen_raw_blast_lists(i.name[:-10]+"_blast.txt")
            #write the parsed blast to fasta. this fasta should have each species_name listed once. 
            i.gen_new_fasta(i.name)
        for e in list_raw_gene_fastas:
            e.gen_original_lists(e.name+"_raw_blast")
            e.gen_species_blast()
            
        #now we have a list of Fasta objects representing the raw fasta files that may exist (or not)
        #and each has 1. ids 2. species names 3. sequences & all indexed exactly the same.
        
        #now i want to 1. make a new .fasta output for each query
        #2. do the thing for each species in query 1, print it.

        
        i = 0
        for QUERY in gene_fasta_list:
            #for each gene you want in the overall list....
            current_overall_fasta = list_overall_gene_fastas[i]
            print("Current overall")
            print(current_overall_fasta)
            print(current_overall_fasta.ids)
            current_raw_fasta = list_raw_gene_fastas[i]
            NewFasta = Fasta()
            #set up the new fasta to write to (NewFasta)
            if prefix is not False:
                NewFastaName = "./"+prefix+"/"+prefix + gene_name_list[i].strip() + ".fasta"
            else:
                NewFastaName = "./"+prefix+"/"+gene_name_list[i].strip() + ".fasta"
            i += 1
            donum = str(i)
            totalnum = str(len(gene_fasta_list))

            print("Beginning on: " + NewFastaName + "which is number: " +donum + " of a total :" + totalnum)
            #for each species....
            for ENTREZ in species_list:
                in_overall = "no"
                no_under = re.sub(" ", "_", ENTREZ.strip("\""))
                #see if its already listed in current_overall_fasta
                for item in current_overall_fasta.ids:
                    #print(item)
                    if no_under in item:
                        #print(item)
                        in_overall = "yes"
                        index = current_overall_fasta.ids.index(item)
                        nseq = current_overall_fasta.seqs[index]
                        item = item.strip()
                        nseq = nseq.lstrip()
                        if nseq == "----":
                        #was searched before and got no valid hits
                            break
                        else:
                            NewFasta.original_seqs.append(nseq)
                            NewFasta.original_ids.append(item)
                            NewFasta.seqs.append(nseq)
                            NewFasta.ids.append(item)
                            break
                    if in_overall == "yes":
                        continue
                #if not in current overall, look in raw_blast
                else:
                    #print("Searching")
                    #do the searc
                    match = False
                    #see if its in current_raw_fasta
                    for spec in current_raw_fasta.species_names:
                        if spec == no_under:
                            match = True
                            index = current_raw_fasta.species_names.index(no_under)
                            nseq = current_raw_fasta.seqs[index]
                            nid = current_raw_fasta.ids[index]
                            nseq = nseq.strip()
                            nid = nid.strip()
                            NewFasta.original_seqs.append(nseq)
                            NewFasta.original_ids.append(no_under)
                            NewFasta.seqs.append(nseq)
                            NewFasta.ids.append(no_under)
                            current_overall_fasta.original_seqs.append(nseq)
                            current_overall_fasta.original_ids.append(no_under)
                            current_overall_fasta.seqs.append(nseq)
                            current_overall_fasta.ids.append(no_under)
                    if match == False:
                        #cry. 
                        print("Could not find: "+no_under)
            current_overall_fasta.shorten()
            NewFasta.shorten()
            end_file_list.append(NewFasta.gen_new_fasta(NewFastaName))
        print("Should be finished creating one fasta file for each gene in your input! They are called: ")
        print(end_file_list)
        for removethis in gene_fasta_list:
            os.system("rm " + removethis)
        # after everything else...
        for d in list_overall_gene_fastas:
            nam = d.name
            d.gen_new_fasta(nam)
        a = muscle_align_on_cluster(end_file_list, prefix)
        print (a)
        return a
        
    #########back to REMOTE SEARCH ######
    for QUERY in gene_fasta_list:
       
        final_e = "no"
        totalnum = str(len(gene_fasta_list))
        current_overall_fasta = list_overall_gene_fastas[i]
        NewFasta = Fasta()
        if prefix is not False:
            NewFastaName = "./"+prefix+"/"+prefix + gene_name_list[i].strip() + ".fasta"
        else:
            NewFastaName = "./"+prefix+"/"+gene_name_list[i].strip() + ".fasta"
        i += 1
        donum = str(i)
        print("Beginning on: " + NewFastaName + "which is number: " +
              donum + " of a total :" + totalnum)
        e = 0
        ten_blasts = ""
        post_b = []
        iteration_blast = 0
        local_spec_list = []
        for item in species_list:
            local_spec_list.append(item)
        for ENTREZ in species_list:
            e += 1
            n += 1
            elen = len(species_list)
            if e == elen:
                final_e = "yes"
            in_overall = "no"
            no_under = re.sub(" ", "_", ENTREZ.strip("\""))
            for item in current_overall_fasta.ids:
                if no_under in item:
                    in_overall = "yes"
                    local_spec_list.remove(ENTREZ)
                    index = current_overall_fasta.ids.index(item)
                    nseq = current_overall_fasta.seqs[index]
                    item = item.strip()
                    nseq = nseq.lstrip()
                    if nseq == "----":
                        #was searched before and got no valid hits
                        break
                    else:
                        NewFasta.original_seqs.append(nseq)
                        NewFasta.original_ids.append(item)
                        NewFasta.seqs.append(nseq)
                        NewFasta.ids.append(item)
                        break
            if in_overall == "yes":
                continue
            else:
                #
                print("Searching")
                OUTPUT = "blastp" + str(n) + ".fasta"
        # run the blast
                # sequentially WORKS rn


                ############JUST FOR NOW######
                ##continue
                ###############


                #do the thing remotely
                blast_query = "blastp -query " + QUERY + " -remote -db nr -out " + OUTPUT +" -max_target_seqs 5 -entrez_query " + ENTREZ + " -evalue 1e-4 -outfmt \"6 sallseqid salltitles sseq\""
                print(blast_query)
                os.system(blast_query)
                print("finished, moving on")
                a = current_overall_fasta.blast2fasta(OUTPUT, ENTREZ, 1)
                if a == "bad":
                    print("Bad was returned, should be saving")
                    pass
                else:
                    NewFasta.blast2fasta(OUTPUT, ENTREZ, 1)
                output_list.append(OUTPUT)
                #######end
        # # nr or refseq_protein database?
                # delete all the newly-generated files
                # after each search:
        # things to do at end of each gene.
        current_overall_fasta.shorten()
        NewFasta.shorten()
        end_file_list.append(NewFasta.gen_new_fasta(NewFastaName))
    print("Should be finished creating one fasta file for each gene in your input! They are called: ")
    print(end_file_list)
    for outp in output_list:
        os.system("rm "+outp)
    for removethis in gene_fasta_list:
        os.system("rm " + removethis)
        # after everything else...
    for d in list_overall_gene_fastas:
        nam = d.name
        d.gen_new_fasta(nam)
    
    a = muscle_align_on_cluster(end_file_list, prefix)
    
    #optimizing:  see how much cores muscle uses
        
    #should gen muscle script
    #should gen correlation file
    #should move the files (end_file_list) to the cluster and wait for "finished" prompt, then move the output files back.
    print (a)
    return a

    # to do this in bulk for eg. my project
    # i will need to save the critters in each species tree in a critter file titled "prefix"
    # but i won't actually have to pass the subtree here.
    # so, subtree gen needs to save 1. a subtree .newick 2. a list of species
    # to include in the species tree 3. a list of all their names 4. a map of
    # which ribo genes to use with each.

    # make another script that catches the list/correlation and runs the above
    # once for each instance.

    #TODO if Species_Trees doesn't exist, make it. if prefix doesn't exist, make it.




def get_multiple_concat_alignments(list_of_whatever, local):
    #pass in multiple sets of [species_list, genes_list, prefixname] and will run blast, align on cluster, concatenate.
    #outputs a listof concatenated species tree files.
    list_of_ccs = []
    #local will be True if should run locally, false otherwise.
    for item in list_of_whatever:
        species_list = item[0]
        gene_sequences_file = item[1]
        if len(item) == 3:
            prefix = item[2]
            a = MakeSpeciesTreeCC(species_list, gene_sequences_file, local, prefix)
        else:
            print("No prefix given!")
            a = MakeSpeciesTreeCC(species_list, gene_sequences_file, local)
        list_of_ccs.append(a)
       
    return list_of_ccs
    

def raxml_run_on_cluster(cc_file_list, prefix):
    #this creates dir you will use on the cluster.
    tree_list = []
    for item in cc_file_list:
        #this is going to be something different... like raxml_bipartitions.blah
        alist = item.split(".")
        athing = alist[0]
        tree_list.append("RAxML_bipartitions."+athing)
    print(tree_list)
    check_directory_existance(prefix, ssh_inst)
    clus_path = "/Species_Trees"
    a = gen_raxml_script(prefix+"_Rax_Sc.sh", "~"+clus_path+"/"+prefix+"/"+prefix+"_Rax_Corr.txt", str(len(cc_file_list)), prefix+"job")
    b = gen_correlate_file(cc_file_list, prefix+"_Rax_Corr.txt")
    cc_file_list.append(a)
    cc_file_list.append(b)
    direct = os.getcwd()
    move_to_cluster(cc_file_list, clus_path+"/"+prefix)
    print("everything should be generated and on the cluster. starting raxml.")
    os.system(ssh_inst+" 'cd ~/Species_Trees/"+prefix+";echo $PWD;sbatch "+a+"'")
    finished = "start"
    #to see if the run is complete, see if each new file has been generated. check every 5 minutes for muscle.
    movehome = []
    ret_list = []
    for i in tree_list:
        movehome.append(i)
    while finished is not True:
        for filename in movehome:
            os.system("scp "+clus_head[:-1]+clus_path+"/"+prefix+"/"+filename+" "+direct+"/"+prefix)
        for item in tree_list:
            #see if it got moved home.
            exists = os.path.isfile("./"+prefix+"/"+item)
            if exists is True:
                if item in movehome:
                    movehome.remove(item)
                finished = "yes"
                #lets move this returned file (c) from /prefix/ to ~ and rename it from RAXML_bipart.blah to blah_RAxML_bipart.newick
                rax, ext = item.split(".")
                os.system("cp ./"+prefix+"/"+item+" .")
                os.system("mv "+item+" "+ext+"_Rax_Bipart.newick")
                ret_list.append(ext+"_Rax_Bipart.newick")
            else:
                finished = False
                print("Rax not done yet. could not locate :./"+prefix+"/"+item+"checking again in 5 minutes")
                break
        if finished == "yes":
            print("Should be done!")
            finished = True
        else:
            #wait ten minutes and then try again.
            time.sleep(300)
            finished = "yes"
    c = ret_list
    print(c)
    print("RAXML finished, tree(s) should exist locally")
    print("removing .fasta, muscle, and slurm.out files from cluster...")
    remove_slurm_files(ssh_inst, prefix, "slurm*")
    remove_slurm_files(ssh_inst, prefix, "*.fasta*")
    return c


def write_species_list(species_list, prefix):
    with open("./"+prefix+"/"+prefix+"_Species_List.txt", "w") as new:
        for item in species_list:
            item = item.strip()
            item = item.strip("\"")
            new.write(item+"\n")



# run the shit



###############
#PARSER #TODO##
# ###############
# import os
# #import socket
# #import getpass
# import re
# import time


# os.chdir("/Users/abigailc/Dropbox/ClusterFriend/Cara/")




# prefix = "Initial_Ribo"
# spec_file = "Species_File_Test.txt"
# gene_file = "Bacterial_Ribo.fasta"

# cc_files = get_multiple_concat_alignments([[spec_file, gene_file, prefix]])
# raxml_run_on_cluster(cc_files, prefix)

###print("super done")

if __name__ == "__main__":

    print("Running in terminal")  
    import sys
    import argparse
    import os
    import re
    import time
    parser = argparse.ArgumentParser(description="All")
    parser.add_argument("directory", nargs='?', default=os.getcwd(), type=str, help="type name of directory to run in (where .nex resides)")
    parser.add_argument("-s", "--species", action = "store", default = False, help="give a species file")
    parser.add_argument("-f", "--fasta", action = "store", default = False, help="give a .fasta file (to generate species list from)")
    parser.add_argument("-g", "--genes", action = "store", default = False, help = "give a .fasta file containing the target genes to include in the species tree (probably ribosomal genes from species similar to those in your fasta)")
    parser.add_argument("-p", "--prefix", action = "store", default = "SpTr", help="Prefix for this species tree")
    parser.add_argument("-l", "--lists", action = "store", default = False, help="use this if you have a list of species trees to make. format: [[spec1, gene1, prefix1], spec2, gene2, prefix2]")
    parser.add_argument("-r", "--runlocal", action = "store_true", default = False, help="run this locally. you must have a database in your path. currently set to nr.")

#useage
#assume you want to use species that are found in the fasta file example.fasta, and genes as found in Euk_Ribo.fasta
#always use -r flag for now, it is faster.
# $python MakesSpeciesTrees.py -f example.fasta -g Euk_Ribo.fasta -p MyTree1 -r


    
    args = parser.parse_args()
    
#making input list
    try:
        os.chdir(args.directory)
    except:
        print ("didn't change dir")
    #do the thing if given a list.

    if args.lists != False:
        cc_files = get_multiple_concat_alignments(args.lists, args.runlocal)
        raxout = raxml_run_on_cluster(cc_files, args.prefix)
        print("Finished")
        raise SystemExit
    #do the thing if not given a list
    if args.species != False:
        species_file = args.species
        species_list = []
        with open(species_file) as species:
            for item in species:
                if item == "":
                    pass
                if item == "\n":
                    pass
                item = item.strip()
                item = item.replace("_", " ")
                item = "\"" + item + "\""
                species_list.append(item)
    else:
        if args.fasta != False:
            species_list = []
            MyFasta = Fasta(args.fasta)
            MyFasta.gen_original_lists(MyFasta.name)
            species_f = MyFasta.gen_species_lists()
            for item in species_f:
                item = item.strip()
                item = item.replace("_", " ")
                item = "\"" + item + "\""
                species_list.append(item)
    if args.genes == False:
        print("Please specify a list of query genes to build the species trees from")
        raise SystemExit
    mylist = []
    while '""' in species_list:
        species_list.remove('""')
    mylist.append(species_list)
    mylist.append(args.genes)
    mylist.append(args.prefix)
    print (mylist)
    #since we are running multiple_cc_align with only one species tree to be made, the list (mylist) needs to be again nested in a list. usually that would be a list of several info for multiple species trees to be made, but since we are just making one, its a list of length 1. don't worry about it.

    cc_files = get_multiple_concat_alignments([mylist], args.runlocal)
    raxout = raxml_run_on_cluster(cc_files, args.prefix)
    

    
##check that loading
#1. fasta
#2. species file
#will work
#list of genes, species, prefix works

    
