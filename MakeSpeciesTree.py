# #!/usr/bin/python

# last edit abigailc@Actaeon on september 14 2016
# lets make an alternative that is usable with EITHER
# 1 downloaded full genomes that you turn into a database
# or
# 2 run against remote NCBI database (currently doing this, because i can't do taxon-restriction against full databse and
# my project needs automated creation of holy wow a shitton of genes.
#
# how this will work
# 1. i need a list of species names. right now I don't care how they are obtained,
# but ill put that on the todo list.
# 2. i need a list of protein sequences we care about, in fasta format.
#
#
# idea to save time; implement before final run through of giant project
# make sure to save each hit in a current-project file, but ALSO append it to an OVERALL file, and
# when we read a new taxa/gene pair, consider searching the OVERALL file
# before running a new blast

# for expediency in problem-solving when sharing script, just keeping bits of Fasta module that are needed contained here.
# try:
#     import FastaClass.py
# except:

######this needs to be set by each person right now. this is how I ssh from terminal.

ssh_inst = "ssh -l abigailc -i ~/.ssh/id_rsa eofe5.mit.edu"
clus_head = "abigailc@eofe5.mit.edu:/home/abigailc/"

##########

def check_directory_existance(prefix, ssh_inst):
    os.system(ssh_inst+" \' mkdir Species_Trees;cd Species_Trees;mkdir "+prefix+"\'")
    

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

    def write_one_seq_per_file(self):
        geneflist = []
        genenames = []
        for i in range(len(self.ids)):
            with open("Seq" + str(i), "w") as new:
                new.write(">" + self.ids[i])
                new.write(self.seqs[i])
                name = re.sub("([^\|]*)(\|)(.*)", "\\1", self.ids[i])
                geneflist.append("Seq" + str(i))
                genenames.append(name)
        return geneflist, genenames
        print("one per file generated")

    def number_of_sites(self):
        return len(self.original_seqs[0])

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
        # entrez is used to ensure that sequence saved uses correct TAXON, esp. if sequence is a MULTISPECIES entry.
        # entrex should be somethin like "Mycobacterium triplex"
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
                    end = "yes"
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
                self.seqs.append(newseq)
                self.original_seqs.append(newseq)

        print("Blasttofasta id/seq loading complete!")

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
        print("Finished, your new fasta file is located at " + newfasta)
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


def reorder(fasta1, fasta2):
    FirstFasta = Fasta("firstfasta")
    FirstFasta.gen_original_lists(fasta1)
    SecondFasta = Fasta("second")
    SecondFasta.gen_original_lists(fasta2)
    for item in FirstFasta.ids:
        index_original = FirstFasta.original_ids.index(item)
        index_two = SecondFasta.original_ids.index(item)
        if index_original == index_two:
            pass
        else:
            SecondFasta.ids[
                index_original] = SecondFasta.original_ids[index_two]
            SecondFasta.seqs[
                index_original] = SecondFasta.original_seqs[index_two]
    newfastaname = fasta2[:-6] + "Reordered.fasta"
    SecondFasta.gen_new_fasta(newfastaname)
    print("should be complete")


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
                id_index_list = []
                id_index_list.append(ids)
                # check the index of that id in each list and append to "id_index_list"
                # if the id is not in that list, should append "NA"
                for eachlist in list_of_ids_lists:
                    try:
                        index = eachlist.index(ids)
                    except ValueError:
                        index = "NA"
                    id_index_list.append(index)
                # add the result of scanning that id to overall output list.
                output_list.append(id_index_list)
    # output list looks like:
    # outputlist = [ ["Cat",1,2,3,4] , ["Dog", 2,1,13,14] ]
    return output_list


# this should work but hasn't been tested. requires stuff above.
def Concatenate(listoffastafiles, new_cc_fasta_name, strain_limit=False):
        # for each thing in listoffastas, create a fasta object in a list.
    fasta_class_list = []
    list_species_lists = []
    # this will create an instance of Fasta from each fasta file
    for i in listoffastafiles:
        fasta_class_list.append(Fasta(i))
        # i am not sure, then, how to call each's inner function, other than
        # one-by-one in a list (but that is fine i think?)

    # this will create the original lists (seqID and sequence) for each Fasta object
    # this will create a list of the "Genus_species" from each seq in each
    # Fasta object (pass true/false of strain control if need be)
    for f in fasta_class_list:
        f.gen_original_lists(f.name)
        species_list = f.gen_species_lists(strain_limit)
        # todo:make this function
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
    with open(new_cc_fasta_name, "w") as new:
        for item in indexed_ids_list:
            # do i want to print species_name to file, or full
            # taxonomy? for now, just species_name is gunna happen,
            # but easy to swap i think
            new.write(">" + item[0])
            fas_num = 0
            for fas in fasta_class_list:
                fas_num += 1
                # fas_num keeps track of what number fasta we are on, which
                # correlates to the index of the index in indexed_ids_list
                search_index = item[fas_num]
                # search_index will be something like "22"

                # if search_index is NA, generate a str of "-" that is n
                # characters long, where n is the return from
                # fas.number_of_sites
                if seach_index == "NA":
                    ndash = fas.number_of_sites
                    retreived_seq = ""
                    for i in range(ndash):
                        retreived_seq = retreived_seq + ("-")
                else:
                    retreived_seq = fas.seqs[search_index]
                    # retreived_seq wil be something like "the 22nd sequence in
                    # object Fas's sequence list... " or "BLAHSEQUENCEDATA"
                new.write(retreived_seq)
                # there might be too many line-breaks in this fashion. if so,
                # concatenated all retreived_seqs and then only call write()
                # once.
        print(
            "Should be finished generating new concatedated fasta at: " + new_cc_fasta_name)
    print("done w cc gen!")
    return new_cc_fasta_name



def muscle_align_on_cluster(end_file_list, prefix):
    #this creates dir you will use on the cluster.
    aligned_list = []
    for item in end_file_list:
        aligned_list.append(item+"_Muscle.fasta")
    check_directory_existance(prefix, ssh_inst)
    clus_path = "/Species_Trees/"+prefix
    a = gen_muscle_script(prefix+"_Sc.sh", "~"+clus_path+"/"+prefix+"_Corr.txt", str(len(end_file_list)), prefix+"job")
    b = gen_correlate_file(end_file_list, prefix+"_Corr.txt")
    end_file_list.append(a)
    end_file_list.append(b)
    move_to_cluster(end_file_list, clus_path)
    print("everything should be generated and on the cluster")
    os.system("ssh -l abigailc -i ~/.ssh/id_rsa eofe5.mit.edu 'cd ~/Species_Trees/"+prefix+";echo $PWD;sbatch "+a+"'")
    finished = "start"
    #to see if the run is complete, see if each new file has been generated. check every 5 minutes for muscle.
    while finished is not True:
        for item in aligned_list:
            exists = os.path.isfile(item)
            if exists is True:
                finished = "yes"
            else:
                finished = False
                break
        if finished == "yes":
            print("Should be done!")
            finished = True
        else:
            time.sleep(300)
    #now, concatenate the output file. when doing large batches, might want to move each to new folder.... but maybe not for now.
    print("Beginning concatenation")
    c = Concatenate(aligned_list, prefix+"_CC.fasta")
    print("Your file is located at :"+c)
    return c

    #now... run_raxml_on_cluster.... take all the c's and run them.another layer of abstraction.
   # gen_raxml_script()
    
def gen_muscle_script(scriptfile, indexname, n, Jobname):
    #currently assuming you are running in the dir that files are in and should be returned to.
    direct = os.getcwd()
    host = socket.gethostname()
    user = getpass.getuser()
    addr = user+"@"+host+":"+direct
    #figure out how many need to be run. n = len(listoffilesinindex). n
    #figure out a name - scriptfile
    #figure out path to index file. indexname
    #thats it. just print the script. return its filename, which will need to be added to list of things to be moved to the cluster.
    

##example script
    a =  """#!/bin/bash                                                                                             
#SBATCH -p sched_mit_g4nier                                                                             
#SBATCH -t 0-00:50:00                                                                                   
#SBATCH -J """+Jobname+"""                                                                                         
##SBATCH -o Jobname.out
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

scp $THE_INPUT_FILE$ENDING """+addr+"""

exit"""
    with open(scriptfile, "w") as script:
        script.write(a)
    return scriptfile




def gen_raxml_script():
    pass


def gen_correlate_file(list_of_input_files, corr_file):
    #this should be in form
    #1 name
    #2 name
    #3 name
    #requires 1: list of files 2. name for corr_file.
    i = 0
    with open(corr_file, "w") as corr:
        for item in list_of_input_files:
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

#obsolete for now, scp moving from cluster should happen automatically at the end of each job.
def move_from_cluster():
    pass

# align with muscle, probably on the cluster in an array script
# to automate:
# get an example script that [greg] made... i think all that is unique is
# that i need to generate a file of type

# i want this program to a) make that file and b) make a unique script
# 1 first_input_name
# 2 second_input_name
# etc and then another post-concat.

# concatenate them all. maybe also on the cluster?
# raxml them all. probably on the cluster in an array script.

# add the new sequence(s) to the list... or a list for verification? #todo implement verification
# delete all the # generated files
#   #save the new .fasta file.

# DONE 1. split the .fasta into a bunch of temporary .fasta files
# DONE 2. for each new .fasta file, run a blast containing each query, and save the top hit (only)
# WORKINGON 2.5(check that top hit species == query name) in another temp fasta file.
# 3. run blast2fasta(self, blastlist) to add your new thing to the persistant class fasta
# 3. print that list to a new .fasta file.
# you will have a bunch of .fasta files like
# ribo1
# ribo2
# ribo3
# etc
# -> run each through a reformatting program like blasttofasta DONE and shorten TODO
# then
# -> align each
# # -> concatenate each#

# how to gen several class instances at once.
# class Enemy:
#     def __init__(self):
#         pass

#     def update(self):
#         pass

#     def draw(self, display_surface):
#         pass
# To create multiple enemies you just need to do create a list of enemies
# where you will store all your reference to your currently active enemy
# objects.

# enemies = []
# To add a single enemy, you just append a call of the constructor of
# enemy to your list. You can do this as many times as you want.

# enemies.append(Enemy())
# As @monkey also mentioned in the comments, you can easily use list
# comprehension to set your list of active enemies to multiple object
# instances at once.

# enemies = [Enemy() for x in range(10)]
# for enemy in enemies: # loop through all your active enemies
#     enemy.update() # Update the state of every enemy.

# for enemy in enemies: # once again loop through all your active enemies
# enemy.draw(display_surface) # Draw the image of every enemy to the
# display surface.





#
#################################
# this is main fxn, should be at bottom
###############################



def MakeSpeciesTree(species_names_file, gene_sequences_file, prefix=False):


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
    # initialize lists:
    species_list = []
    gene_fasta_list = []
    end_file_list = []
    # initialize fasta instance for genes
    GeneFasta = Fasta()
    GeneFasta.gen_original_lists(gene_sequences_file)
    # write one fasta file per sequence
    gene_fasta_list, gene_name_list = GeneFasta.write_one_seq_per_file()
    # gene_fasta_list contins the names of each auto-generated .fasta file, each of which contains exactly one sequence... each representing a gene you want in your gene tree.
# TODO - a toggle between Archaea/Bacteria/Eukarya because there will be
# different ribo proteins in each.

    # split up the species into a list
    with open(species_names_file) as species:
        for item in species:
            item = item.strip()
            item = "\"" + item + "\""
            species_list.append(item)
    # overallfastas... gene_name_list. as "a"
    list_overall_gene_fastas = []
    gene_overall_fasta_names = []
    for b in range(len(gene_name_list)):
        gene_name = gene_name_list[b].strip()
        gene_overall_fasta = gene_name + "_All.fasta"
        gene_overall_fasta_names.append(gene_overall_fasta)
        list_overall_gene_fastas.append(Fasta(gene_overall_fasta))
    for c in list_overall_gene_fastas:
        c.gen_original_lists(c.name)

    n = 0
    i = 0
    for QUERY in gene_fasta_list:
        final_e = "no"
        totalnum = str(len(gene_fasta_list))
        current_overall_fasta = list_overall_gene_fastas[i]
        NewFasta = Fasta()
        if prefix is not False:
            NewFastaName = prefix + gene_name_list[i].strip() + ".fasta"
        else:
            NewFastaName = gene_name_list[i].strip() + ".fasta"
        i += 1
        donum = str(i)
        print("Beginning on: " + NewFastaName + "which is number: " +
              donum + " of a total :" + totalnum)
        e = 0
        ten_blasts = ""
        post_b = []
        iteration_blast = 0
        for ENTREZ in species_list:
            e += 1
            elen = len(species_list)
            if e == elen:
                final_e = "yes"
            in_overall = "no"
            no_under = re.sub(" ", "_", ENTREZ.strip("\""))
            for item in current_overall_fasta.ids:
                if no_under in item:
                    in_overall = "yes"
                    
                    index = current_overall_fasta.ids.index(item)
                    nseq = current_overall_fasta.seqs[index]
                    item = item.strip()
                    nseq = nseq.lstrip()
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
        # maybe several at once? unclear preformance increase.
        # set to nr database for now
        # set to one target sequence... will only pull the best hit.

        # maybe get all of these as lists, and run multiple at a time managed by
        # another function...?

        # wait until i have ten query/output/entrez pairs, and then execute all of
        # them...?
# trying to speed this up is not working. converting to iteratively.
                # iteration_blast += 1
                # n += 1
                # single = "blastp -query " + QUERY + " -remote -db nr -out " + OUTPUT + \
                #     " -max_target_seqs 5 -entrez_query " + ENTREZ + \
                #     " -evalue 1e-4 -outfmt \"6 sallseqid salltitles sseq\" &"
                # #print (ten_blasts)
                # ten_blasts = ten_blasts + single
                # post_b.append([OUTPUT, ENTREZ])
                # if int(iteration_blast) > int(3):
                #     ten_blasts = ten_blasts + " wait"
                #     print(ten_blasts)
                #     os.system(ten_blasts)
                #     ten_blasts = ""
                #     post_b = []
                #     iteration_blast = 0
                #     for pair in post_b:
                #         current_overall_fasta.blast2fasta(pair[0], pair[1], 1)
                #         NewFasta.blast2fasta(pair[0], pair[1], 1)
                #         os.system("rm " + OUTPUT)
                # if final_e == "yes":
                #     # run anyways
                #     ten_blasts = ten_blasts + " wait"
                #     os.system(ten_blasts)
                #     ten_blasts = ""
                #     iteration_blast = 0
                #     for pair in post_b:
                #         current_overall_fasta.blast2fasta(pair[0], pair[1], 1)
                #         NewFasta.blast2fasta(pair[0], pair[1], 1)
                #         os.system("rm " + pair[0])
                
                # sequentially WORKS rn


                ############JUST FOR NOW######
                continue
            ###############
                os.system("blastp -query " + QUERY + " -remote -db nr -out " + OUTPUT +
                          " -max_target_seqs 5 -entrez_query " + ENTREZ + " -evalue 1e-4 -outfmt \"6 sallseqid salltitles sseq\"")
                current_overall_fasta.blast2fasta(OUTPUT, ENTREZ, 1)
                NewFasta.blast2fasta(OUTPUT, ENTREZ, 1)
                os.system("rm " + OUTPUT)

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
    for removethis in gene_fasta_list:
        os.system("rm " + removethis)
        # after everything else...
    for d in list_overall_gene_fastas:
        nam = d.name
        # print(d.seqs[0])
        # print("seq")
        # print(d.ids[0])
        # print("id")
        # print("was space?")
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









# run the shit
import os
import socket
import getpass
import re
import time


os.chdir("/Users/abigailc/Dropbox/ClusterFriend/Cara/")
MakeSpeciesTree("Species_File_Test.txt", "Bacterial_Ribo.fasta", "Example")
print("super done")

##errr w call gen muscle script
###err w gen fastas if in pre-screening
