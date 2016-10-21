# #!/usr/bin/python

# last edit abigailc@Actaeon on october 19 2016

#this script will take an alignment, send it to muscle align and raxml on the cluster, and then bring it home.

#set these yourself
ssh_inst = "ssh -l abigailc -i ~/.ssh/id_rsa eofe5.mit.edu"
clus_head = "abigailc@eofe5.mit.edu:/home/abigailc/"

#imports
import sys
import argparse
import os
import re
import time

#makes dir if need be
def check_directory_existance(ssh_inst):
	import os
	print("checking dirs")
	os.system(ssh_inst+" \'mkdir MusRax\'")

#removes extra files
def remove_slurm_files(ssh_inst,pattern):
	print("removing")
	os.system(ssh_inst+" \'cd MusRax; rm "+pattern+"\'")

#does everything
def musrax_on_cluster(filename, scriptname):
	print("starting")
	#define muscle-out, rax-out names
	outname = filename+"_Muscle.fasta"
	raxl = outname.split(".")
	raxn = raxl[0]
	raxoutname = "RAxML_bipartitions."+raxn
	#make dir on cluster if need be
	check_directory_existance(ssh_inst)
	#artefact
	clus_path = "/MusRax"
	#make the script
	a = gen_musrax_script(scriptname, filename, outname, raxn)
	#current dir
	direct = os.getcwd()
	#move files to cluster
	move_to_cluster([filename,a], clus_path)
	#submit the .sh file on the cluster
	os.system(ssh_inst+" 'cd ~/MusRax/;echo $PWD;sbatch "+a+"'")
	finished = "start"
	#to see if the run is complete, see if each new file has been generated. check every 5 minutes for raxml output file name
	while finished is not True:
		#try and copy it home
		os.system("scp "+clus_head[:-1]+clus_path+"/"+raxoutname+" "+direct)
		#see if it exists at home
		exists = os.path.isfile(raxoutname)
		if exists is True:
			finished = "yes"
		#if not, wait and try again
		else:
			finished = False
			print("waiting 5 minutes and trying again")
			time.sleep(300)
		if finished == "yes":
			print("Should be done!")
			finished = True
	print("Your file should exist")
	return raxoutname


#this creates a single-use script to run muscle and raxml (gammaLG + 100 rapid bootstraps + consensus tree)
def gen_musrax_script(scriptfile, inputfile, outputfile, raxn):


##example script
	a =  """#!/bin/bash																							 
#SBATCH -p sched_mit_g4nier																			 
#SBATCH -t 0-01:00:00	
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH -J RAX"""+scriptfile+"""   
#SBATCH -o RAX"""+scriptfile+""".out																						 

. /etc/profile.d/modules.sh
module add engaging/RAxML/8.2.9
module add engaging/muscle/3.8.31
muscle -in """+inputfile+""" -out """+outputfile+"""
raxmlHPC-PTHREADS-AVX -T 20 -f a -m PROTGAMMALG -p 12345 -x 12345 -#100 -n """+raxn+""" -s """+outputfile+"""		 

exit"""
	
	with open(scriptfile+".sh", "w") as script:
		script.write(a)
	return scriptfile+".sh"

#moves your .fasta and script to cluster
def move_to_cluster(list_of_files, clus_path):
	#requires os.system
	#requires scp(?)
	#do the thing
	for item in list_of_files:
		os.system("scp "+item+" "+clus_head+clus_path)
	print("Finished moving files to cluster in place:"+clus_path)

# parser

if __name__ == "__main__":

	print("Running in terminal")  
	import sys
	import argparse
	import os
	import re
	import time
	parser = argparse.ArgumentParser(description="All")
	parser.add_argument("directory", nargs='?', default=os.getcwd(), type=str, help="type name of directory to run in (where .nex resides)")
	parser.add_argument("-s", "--script", action = "store", default = False, help="give a name for your script/jop")
	parser.add_argument("-f", "--fasta", action = "store", default = False, help="give a .fasta file")


	
	args = parser.parse_args()
	#change dir if given
	try:
		os.chdir(args.directory)
	except:
		print ("didn't change dir")
	#run the thing
	musrax_on_cluster(args.fasta, args.script)
	print("done")
