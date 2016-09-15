#!/usr/bin/python
#this is a conglomerate of fasta-fixing scripts, now called FISH (FASTA ID SWAPPING HELPER) because lol i can acronym.
##
#last edit abigailc@Actaeon Sept 7 2016

#things this doesn't do: play super nice with accession numbers instead of GI numbers. probably easy to convert, (see that one script that one time), but meh
#do it later.

class Fasta:
	def __init__(self, name):
		#all ids should be stripped and have ">" removed for reasons.
		#for now, sequences do not have any stripping applied
		self.name = name
		self.ids = []
		self.original_ids = []
		self.original_seqs = []
		self.seqs = []
	def gen_original_lists(self, fastaname):
		with open(fastaname) as fastafile:
			for line in fastafile:
				if ">" in line:
					#write the previous AA seq
					try:
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
					self.ids.append(newline)
					self.original_ids.append(newline)
				else:
					AAseq = AAseq+line
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
	def common_shorten(self):
		#TODO: allow input of manual shorten-pairs, possibly in new function
		#put your conversions of common strings to shorten here
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
			newline = re.sub("bacteriales\|", "bacles|", newline)
			newline = re.sub("bacteriales\|", "bacles|", newline)
			#newline = re.sub("[+=\.]", "", newline)
			newline = re.sub("_enterica_subsp_enterica_serovar", "_ent", newline)
			self.ids[index] = newline
		print("Common shorten complete")
		#this should have successfully modified the self.ids list to contain shortened sequence ids.
	def length_check(self, length, verbose=False):
		#needs to pass in a number... charnum
		toolong = 0
		for item in self.ids:
			index = self.ids.index(item)
			linelength = len(item)
			if linelength > length:
				toolong +=1
				#change all 12 to 14 if include \n at end of seqids... for now, we are not considering them.
				gi = newline[-12:]
				rest = re.sub("([^#/])(.*)", "\\1", newline)
				nogi = rest[:-3]
				newl = length-12
				newnogi = nogi[:newl]
				newline = newnogi+gi
				if verbose == True:
					print ("LENGTHERROR: "+line[:length]+" || "+line[length:])
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
				listoflines.append(line)
				self.ids[index] = rep
		if verb == True:
			print ("there were "+str(num)+" duplicate sequences that were numbered")
		#done
		print("duplicate check done")

	def index_shorted(self, replace):
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
							#	 newid = str(newid)+"|"+str(taxlist[int(item)-1])
							# if f == 3:
							#	 newid = str(newid)+"|"+str(taxlist[int(item)])
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

##					#i have something like
##					>Methanococcoides_burtonii|gi|909890
##					#i want
##					MethBurt00
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
##						print(GenusSpecies)
				gs8 = GenusSpecies[1:9]
				if iteration < 10:
					newid = gs8+"0"+str(iteration)
				elif iteration > 99:
					newid = gs8[:-1]+str(iteration)
				else:
					newid = gs8+str(iteration)
##						print(newid)

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
		with open (infofile) as old:
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
		for item in self.original_ids:
			index = self.original_ids.index(item)
			newid = CTdict[item]
			self.ids[index] = newestid
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
				new.write(">"+self.ids[i])
				# print(i)		#
				#unclear if this needs a "\n" after it... check.#TODO
				new.write(self.seqs[i])
		print("Finished, your new fasta file is located at "+newfasta)
		#done

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

#this hasn't been implemented in class fasta, so I am leaving it commented out.. subtrees file might be easily replaced using replace.newick but it might take literally ages... unclear.

# def replace2(replace_file, dict_old_new, charnum, verb):
#	 print("Recognized subtrees file, using subtrees varient")
#	 outputlist = []
#	 rep = 0
#	 replist = []
#	 newfilename = replace_file.split(".")
#	 newfilename = newfilename[0]+str(charnum)+"limit."+newfilename[1]
#	 with open(replace_file) as old:
#		 if verb == True:
#			 print("Opening "+replace_file)
#		 with open (newfilename, "w") as new:
#			 for line in old:
#				 line = line.strip()
#				 for item in dict_old_new:
#					 if item[:127] in line:
#						 if item[:127] in replist:
#							 pass
#						 else:
#							 replist.append(item[:127])
#						 rep+=1
# ##							print(line)
#						 oldline = line
#						 line = line.replace(item[:127], dict_old_new[item])
# ##						if verb == True:
# ##							if len(line) <200:
# ##								print oldline
# ##								print item
# ##								print dict_old_new[item]
# ##								print(line)
# ##								print("\n")
# ##							print("\n")
#				 new.write(line+"\n")
#			 print("finished with "+newfilename+"made "+str(rep)+" replacements of "+str(len(replist))+" differnt patterns")
# ##				print(replist)

#	 return newfilename



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

	#	 #write stuff
	# def gen_new_fasta(new_fasta_name):
	# def swap_in_nexus():
	# def swap_in_newick(old_newick_name, new_file_name):
	# def gen_info(info_file_name):


if __name__ == "__main__":

	print("Running in terminal")
	import sys
	import argparse
	import os
	import re
	parser = argparse.ArgumentParser(description="All")
	#necessary bits
	parser.add_argument("directory", nargs='?', default=os.getcwd(), type=str, help="type name of directory to run in where fasta resides, if not pwd")
	parser.add_argument("fasta", type=str, help="type the name of your .fasta file")
	#options to load changes from another file
	parser.add_argument("-i", "--infofile", action = "store", default = False, help="Provide an Info File (as generated by this script previously) to pull original and new sequences from")

	#options#  to check,fix,edit,etc the seqs or seqids
	# -length
	# -duplicate
	# -weirdaa
	# -weirdID

	parser.add_argument("-l", "--length", action = "store", default=False, help="Provide a max length for your sequenceIDs")
	parser.add_argument("-d", "--duplicates", action = "store_true", help="Flag causes identical seqIDs to be numbered 1 2 3 etc to prevent program confusion")
	parser.add_argument("-fid", "--fixID", action = "store_true", help="Flag scans SeqIDs and removes weird characters like += etc")
	parser.add_argument("-faa", "--fixAA", action = "store_true", help="Flag scans Sequences and removes non-standard AA characters like X B &")
	#options to shorten specific words
	# -m manual_shorten
	# -c common_shorten

	parser.add_argument("-c", "--common", action = "store_true", help="Flag causes seqIDs to be shortened in a predefined manner, eg bacteriales->bacles ")
	parser.add_argument("-m", "--manual", default = False, action = "store", help="Provide a list of \"original,new\" things to shorten. eg \"Bacteria,Bac Eukaryota,Euk\"")
	#special shortening methods

	parser.add_argument("-t", "--tenchars", action = "store_true", help="Flag turns sequence IDs into ten character strings")
	parser.add_argument("-b", "--bayes", action = "store_true", help="Flag turns sequences into form that will work as MrBayes input")
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

#change dir if desired
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
	MyFasta.gen_original_lists(args.fasta)


	#this should be done in conjunction w / write fasta or replace newick.
	if args.infofile != False:
		MyFasta.load_info_swap(args.infofile)

	#here are the error-fixing calls
	if args.duplicates == True:
		MyFasta.duplicates_check(verb)
	if args.fixID == True:
		MyFasta.weird_ID_check(verb)
	if args.fixAA == True:
		MyFasta.weird_AA_check(verb)
	#shortening calls
	if args.common == True:
		MyFasta.common_shorten()
	if args.manual != False:
		MyFasta.manual_shorten(args.manual)

	if args.piece != False:
		MyFasta.index_shorted(args.piece)

	if args.length != False:
		MyFasta.length_check(args.length, verbose=False)

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
