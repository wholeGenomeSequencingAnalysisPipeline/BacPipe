# coding=utf-8
# coding=utf-8
from appJar import gui
import time
import os
import subprocess
import sys
from sys import argv
from subprocess import call
import glob
import shutil
from argparse import (ArgumentParser, FileType)
import logging
import yaml
import re
import thread
import threading
from threading import Thread

global pipeline_path
pipeline_path=os.getcwd()

#Unix_mac
global Plasmidfinder_samples
Plasmidfinder_samples=""
global resfinder_samples
resfinder_samples=""
global virulence_samples
virulence_samples=""
global mlst_samples
mlst_samples=""
global emm_samples
emm_samples=""
global resfams_samples
resfams_samples=""
global virdb_samples
virdb_samples=""
global card_samples
card_samples=""
global OSystem
OSystem= sys.argv[1] 
global cutadapt
global cutadapt_path
global trimgalore_path
global spades_path
global mlst_ncbi_path
global plasmidfinder_ncbi_path
global resistancefinder_ncbi_path
global barrnap_path
global hmmer_path
global blastp_ncbi_path
global parSNP_path
global sample_name
global parSNP_tree
parSNP_tree=""

sample_name=""
if OSystem == "mac":
	cutadapt=os.path.join(pipeline_path+'/cutadapt-1.12/cutadapt/cutadapt')
	cutadapt_path=os.path.join(pipeline_path+'/cutadapt-1.12/')
	trimgalore_path=os.path.join(pipeline_path+'/trim_galore_v0.4.2/trim_galore')
	spades_path=os.path.join(pipeline_path+'/SPAdes-3.9.1-Darwin/bin/spades.py')
	mlst_ncbi_path=os.path.join(pipeline_path+'/mlst/blast-2.2.26/')
	plasmidfinder_ncbi_path=os.path.join(pipeline_path+'/mlst/blast-2.2.26/')
	resistancefinder_ncbi_path=os.path.join(pipeline_path+'/mlst/blast-2.2.26/')
	barrnap_path=str(pipeline_path+'/barrnap-master/bin/')
	hmmer_path=os.path.join(pipeline_path+'/hmmer-3.1b2-macosx-intel/binaries/hmmscan')
	blastp_ncbi_path=os.path.join(pipeline_path+'/mlst/blast-2.2.26/bin/blastp')
	parSNP_path=os.path.join(pipeline_path+'/parsnp_OSX64_v1_2/parsnp')
elif OSystem == "unix":
	cutadapt=str(pipeline_path+'/cutadapt-1.12/cutadapt/cutadapt')
	cutadapt_path=str(pipeline_path+'/cutadapt-1.12/')
	trimgalore_path=str(pipeline_path+'/trim_galore_v0.4.2/trim_galore')
	spades_path=os.path.join(pipeline_path+'/SPAdes-3.13.0-Linux/bin/spades.py')
	mlst_ncbi_path=os.path.join(pipeline_path+'/mlst/blast-2.2.26/')
	plasmidfinder_ncbi_path=os.path.join(pipeline_path+'/mlst/blast-2.2.26/')
	resistancefinder_ncbi_path=os.path.join(pipeline_path+'/mlst/blast-2.2.26/')
	barrnap_path=str(pipeline_path+'/barrnap_0_6/bin/')
	hmmer_path=os.path.join(pipeline_path+'/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmscan')
	blastp_ncbi_path=os.path.join(pipeline_path+'/mlst/blast-2.2.26/bin/blastp')
	parSNP_path=os.path.join(pipeline_path+'/parsnp_Linux64_v1_2/parsnp')
else:
	print "Please select unix or mac as first argument. E.g. python Pipeline.py mac"
	sys.exit()
#define variables
global counter
counter=0
global mlst_switch
global plasmids_finder_switch
global virulencefinder_switch
global resfinder_switch
global emmTyping_switch
global parSNP_switch
global resfams_switch
global cardSearch_switch
global VirDBSearch_switch

global Quit_Flag
Quit_Flag=0
global input_path
input_path=""
global output_path
output_path=""
global prokka_full_path
prokka_full_path="/BacPipe/prokka/bin/prokka"
global parSNP_reference_path
parSNP_reference_path=""
global parSNP_reference_fsa_path
parSNP_reference_fsa_path=""
global TG_switch
TG_switch=""
global TG_paired
TG_paired=""
global SPAdes_mode
SPAdes_mode=""
global spades_switch
spades_switch=""
global plasmids_finder_
plasmids_finder_=""
global prokka_output
prokka_output=""
global resistanceFinder_output
resistanceFinder_output=""
global mlst_output
mlst_output=""
global plasmidsFinder_output
plasmidsFinder_output=""
global virulenceFinder_output
virulenceFinder_output=""


# function called by pressing the buttons
def trim_galore_help(habal):
	app.infoBox("TrimGalore help", """Trim Galore: quality trimming\n
INPUT:\n  *_R1_001.fastq.gz  \n*_R2_001.fastq.gz\n
OUTPUT: \n  *_R1_001_val_1.fq.gz \n  *_R2_001_val_2.fq.gz""")
def spades_help(habal):
	app.infoBox("SPAdes help", """SPAdes: assembly and scaffolding\n
INPUT:\n  *_R1_001_val_1.fq.gz\n  *_R2_001_val_2.fq.gz\n
OUTPUT:\n  *.fasta""")
def MLST_help(habal):
	app.infoBox("MLST help", """MLST typing:\n
INPUT:\n  *.fasta (DNA sequences)""")
def emm_help(habal):
	app.infoBox("emm Typing help", """emm typing (for Streptococcus):\n
INPUT:\n  *.fasta (DNA sequences)""")
def Plasmindfinder_help(habal):
	app.infoBox("Plasmindfinder help", """Plasmids Finder:\n
INPUT:\n *.fasta (DNA sequences)""")
def resfinder_help(habal):
	app.infoBox("ResFinder help", """ResFinder (Plasmid mediated Resistance only):\n
INPUT:\n  *.fasta (DNA sequences)""")
def CARD_help(habal):
	app.infoBox("CARD help", """CARD help:\n
INPUT: *.fa (protein sequences)""")
def virulenceFinder_help(habal):
	app.infoBox("VirulenceFinder help", """Virulence Finder:\n
INPUT:\n  *.fasta (DNA sequences)""")
def VirDB_help(habal):
	app.infoBox("VirDB help", """VirDB search:\n
INPUT:\n  *.faa (poteins sequences)""")
def BarSNP_help(habal):
	app.infoBox("BarSNP help", """BarSNP: SNP-based phylogy assessment\n
INPUT:\n 2 or more fasta files """)
def Prokka_help(habal):
	app.infoBox("Prokka help", """Prokka (annotation):\n
INPUT:\n  *.fasta (DNA sequences)\n
  PATH to Prokka executable""")
def resFams_help(habal):
	app.infoBox("ResFams help", """ResFams search (extensive resistance search):\n
INPUT:\n  *.faa (poteins sequences)""")
def Summarize_help(habal):
	app.infoBox("Summary help", """Summarize output in excel file:\n
OUTPUT:\nExcel file for each sample (each selected tools shown in one sheet)""")
def Proc_help(habal):
	app.infoBox("Processors help", """Please specify number of processors\n""")

####################
def browse_input (habal):
	global input_path
	input_path=app.directoryBox()
def browse_output (habal):
	global output_path
	output_path=app.directoryBox()
def browse_prokka (habal):
	global prokka_full_path
	prokka_full_path=app.directoryBox()+"/bin/prokka"
def browse_parSNP (habal):
	global parSNP_reference_path
	#parSNP_reference_path=app.directoryBox()
	parSNP_reference_path=app.openBox(title="Please select Genbank file",fileTypes=[('Genbank', '*.gbk')], asFile=False)
def browse_parSNP_fsa (habal):
	global parSNP_reference_fsa_path
	#parSNP_reference_fsa_path=app.directoryBox()
	parSNP_reference_fsa_path=app.openBox(title="Please select fasta file",fileTypes=[('fasta', '*.fsa'),('fasta', '*.fasta'),('fasta', '*.fa')], asFile=False)




def show_results():
	app.setTabbedFrameDisabledTab("TabbedFrame", "Results", disabled=False)
	global output_path
	if(mlst_switch=="no"):
		try:
		 #global output_path
			global mlst_samples
			mlst_output_ = output_path+"/mlst_output.txt"
			try:
				os.remove(mlst_output_)
			except:
				print ""
			Output_organizer_path=os.path.join(pipeline_path+'/Output_organizer.pl')
			call(["perl",Output_organizer_path,mlst_samples,mlst_output_,"none"])

			#pipe = subprocess.Popen(["perl", Output_organizer_path,mlst_samples,mlst_output_], stdout=subprocess.PIPE)
			file=open (mlst_output_,'r')
			array = file.readlines()

			i=0
			try:
				app.addGridRow("Results grid", ["MLST", " "])
			except:
				   print ""
			while i < len(array)+1:
				tmp = array[i].rstrip().split("	")
				try:
					app.addGridRow("Results grid", tmp)
				except:
					print ""
				
				i=i+1
		except:
			print "Can't find MLST output"


	if(plasmids_finder_switch=="no"):
		try:
		 #global output_path
			global Plasmidfinder_samples
			plasmidsFinder_output_ = output_path+"/plasmids_finder_output.txt"
			try:
				os.remove(plasmidsFinder_output_)
			except:
				print ""
			Output_organizer_path=os.path.join(pipeline_path+'/Output_organizer.pl')
			call(["perl",Output_organizer_path,Plasmidfinder_samples,plasmidsFinder_output_,"none"])

			#pipe = subprocess.Popen(["perl", Output_organizer_path,Plasmidfinder_samples,plasmidsFinder_output_], stdout=subprocess.PIPE)

			file=open(plasmidsFinder_output_,'r')
			array = file.readlines()

			i=0
			try:
				app.addGridRow("Results grid", ["Plasmids", " "])
			except:
				   print "3abat"
			while i < len(array)+1:
				tmp = array[i].rstrip().split("	")
				try:
					app.addGridRow("Results grid", tmp)
				except:
					print ""
				
				i=i+1
		except:
			print "Can't find PlasmidFinder output"

	if(virulencefinder_switch=="no"):
		try:
		 #global output_path
			global virulence_samples
			virulencefinder_output_ = output_path+"/virulencefinder_output.txt"
			try:
				os.remove(virulencefinder_output_)
			except:
				print ""
			Output_organizer_path=os.path.join(pipeline_path+'/Output_organizer.pl')
			call(["perl",Output_organizer_path,virulence_samples,virulencefinder_output_,"none"])

			#pipe = subprocess.Popen(["perl", Output_organizer_path,virulence_samples,virulencefinder_output_], stdout=subprocess.PIPE)
			file=open (virulencefinder_output_,'r')
			array = file.readlines()

			i=0
			try:
				app.addGridRow("Results grid", ["Virulence", " "])
			except:
				   print ""
			while i < len(array)+1:
				tmp = array[i].rstrip().split("	")
				try:
					app.addGridRow("Results grid", tmp)
				except:
					print ""
				
				i=i+1
		except:
			print "Can't find VirulenceFinder output"

	if(resfinder_switch=="no"):
		try:
		 #global output_path
			global resfinder_samples
			resfinder_output_ = output_path+"/resfinder_output.txt"
			try:
				os.remove(resfinder_output_)
			except:
				print ""
			Output_organizer_path=os.path.join(pipeline_path+'/Output_organizer.pl')
			call(["perl",Output_organizer_path,resfinder_samples,resfinder_output_,"none"])

			#pipe = subprocess.Popen(["perl", Output_organizer_path,resfinder_samples,resfinder_output_], stdout=subprocess.PIPE)
			file=open (resfinder_output_,'r')
			array = file.readlines()

			i=0
			try:
				app.addGridRow("Results grid", ["Resistance", " "])
			except:
				   print ""
			while i < len(array)+1:
				tmp = array[i].rstrip().split("	")
				try:
					app.addGridRow("Results grid", tmp)
				except:
					print ""
				
				i=i+1
		except:
			print "Can't find ResFinder output"

	if(emmTyping_switch=="no"):
		try:
		 #global output_path
			global emm_samples
			emm_output_ = output_path+"/emm_output.txt"
			try:
				os.remove(emm_output_)
			except:
				print ""
			Output_organizer_path=os.path.join(pipeline_path+'/Output_organizer.pl')
			call(["perl",Output_organizer_path,emm_samples,emm_output_,"none"])

			#pipe = subprocess.Popen(["perl", Output_organizer_path,emm_samples,emm_output_], stdout=subprocess.PIPE)
			file=open (emm_output_,'r')
			array = file.readlines()

			i=0
			try:
				app.addGridRow("Results grid", ["EMM", " "])
			except:
				   print ""
			while i < len(array)+1:
				tmp = array[i].rstrip().split("	")
				try:
					app.addGridRow("Results grid", tmp)
				except:
					print ""
				
				i=i+1
		except:
			print "Can't find EMM output"

	if(parSNP_switch=="no"):
		global parSNP_tree
	 	tree_file=parSNP_tree+'.txt'
		try:
			os.remove(tree_file)
		except:
			print ""
		#Newick_visualization_path=os.path.join(pipeline_path+'/Newick_visualization.pl')
		#call(["perl",Newick_visualization_path,parSNP_tree])
		if os.path.isfile(tree_file):
			f = open(tree_file, 'r')
			read_data = f.read()
			f.close()

			app.clearTextArea("tree")
			app.setTextArea("tree", read_data)
		else:
			print "ParSNP tree not generated"


	if(resfams_switch=="no"):
		try:
		 #global output_path
			global resfams_samples
			resfams_output_ = output_path+"/resfams_output.txt"
			try:
				os.remove(resfams_output_)
			except:
				print ""
			Output_organizer_path=os.path.join(pipeline_path+'/Output_organizer.pl')
			call(["perl",Output_organizer_path,resfams_samples,resfams_output_,"resfams"])

			file=open (resfams_output_,'r')
			array = file.readlines()

			i=0
			try:
				app.addGridRow("Results grid", ["Resistance (ResFams)", " "])
			except:
				   print ""
			while i < len(array)+1:
				tmp = array[i].rstrip().split("	")
				try:
					app.addGridRow("Results grid", tmp)
				except:
					print ""
				
				i=i+1
		except:
			print "Can't find ResFams output"


	if(cardSearch_switch=="no"):
		try:
		 #global output_path
			global card_samples
			cardSearch_output_ = output_path+"/cardSearch_output.txt"
			try:
				os.remove(cardSearch_output_)
			except:
				print ""
			Output_organizer_path=os.path.join(pipeline_path+'/Output_organizer.pl')
			call(["perl",Output_organizer_path,card_samples,cardSearch_output_,"card"])

			file=open (cardSearch_output_,'r')
			array = file.readlines()

			i=0
			try:
				app.addGridRow("Results grid", ["Resistance (CARD)", " "])
			except:
				   print ""
			while i < len(array)+1:
				tmp = array[i].rstrip().split("	")
				try:
					app.addGridRow("Results grid", tmp)
				except:
					print ""
				
				i=i+1
		except:
			print "Can't find CARD output"

	if(VirDBSearch_switch=="no"):
		try:
		 #global output_path
			global virdb_samples
			VirDBSearch_output_ = output_path+"/VirDBSearch_output.txt"
			try:
				os.remove(VirDBSearch_output_)
			except:
				print ""
			Output_organizer_path=os.path.join(pipeline_path+'/Output_organizer.pl')
			call(["perl",Output_organizer_path,virdb_samples,VirDBSearch_output_,"vfdb"])

			file=open (VirDBSearch_output_,'r')
			array = file.readlines()

			i=0
			try:
				app.addGridRow("Results grid", ["Virulence (VFDB)", " "])
			except:
				   print ""
			while i < len(array)+1:
				tmp = array[i].rstrip().split("	")
				try:
					app.addGridRow("Results grid", tmp)
				except:
					print ""
				
				i=i+1
		except:
			print "Can't find CARD output"

#####
#####
#####

'''
def Resistance_detail(btn):
 	global counter
	global resistanceFinder_output
	Piechart_path=os.path.join(pipeline_path+'/PieChart.pl')
	PieChart_out=os.path.join(pipeline_path+'/PieChart_out.txt')
	counter_=str(counter)
	pipe = subprocess.Popen(["perl", Piechart_path,resistanceFinder_output,PieChart_out,counter_], stdout=subprocess.PIPE)
	label1="Resistance genes"+str(counter)
	try:
		app.destroySubWindow(label1)
	except:
		print ""
	app.startSubWindow(label1, modal=True)
	label="Resistance genes details"+str(counter)
	app.addLabel(label, "Resistance genes details")

	res_f = open(PieChart_out, 'r')
	resistance_results = res_f.read()
	res_f.close()
	exec(resistance_results) in globals(), locals()
	app.stopSubWindow()
	app.showSubWindow(label1)
	counter=counter+1
'''
def check_progress(short_log,out_log,counter_steps):
	app.setTabbedFrameSelectedTab("TabbedFrame","Progress")
	app.setTextArea("t1", short_log)
	global Quit_Flag
	while Quit_Flag==0:
		current_=0
		app.setTextArea("t1", "Analysis log")
		read_data=""
		if os.path.isfile(short_log):
			f = open(short_log, 'r')
			read_data = f.read()
			f.close()
		app.clearTextArea("t1")
		app.setTextArea("t1", read_data)
		current_sample_name=sample_name
		if 'parSNP finished '+sample_name in read_data:
		#app.setMeter("progress_", 100, text="Run finished")
			current_=current_+1
		if 'Summerization finished '+sample_name in read_data:
		#app.setMeter("progress_", 100, text="Run finished")
			current_=current_+1
		if 'VirDB search finished '+sample_name in read_data:
		#app.setMeter("progress_", 95,  text=None)
			current_=current_+1
		if 'CARD search finished '+sample_name in read_data:
		#app.setMeter("progress_", 90,  text=None)
			current_=current_+1
		if 'Resfams finished '+sample_name in read_data:
		#app.setMeter("progress_", 80, text=None)
			current_=current_+1
		if 'VirulenceFinder finished '+sample_name in read_data:
		#app.setMeter("progress_", 70, text=None)
			current_=current_+1
		if 'Annotation finished '+sample_name in read_data:
		#app.setMeter("progress_", 60, text=None)
			current_=current_+1
		if 'Resistance finished '+sample_name in read_data:
		#app.setMeter("progress_", 50, text=None)
			current_=current_+1
		if 'PlasmidFinder finished '+sample_name in read_data:
		#app.setMeter("progress_", 40, text=None)
			current_=current_+1
		if 'MLST finished '+sample_name in read_data:
		#app.setMeter("progress_", 30, text=None)
			current_=current_+1
		if 'SPAdes finished '+sample_name in read_data:
		#app.setMeter("progress_", 20, text=None)
			current_=current_+1
		if 'TrimGalore finished '+sample_name in read_data:
		#app.setMeter("progress_", 10, text=None)
			current_=current_+1
		if 'emmTyping search finished '+sample_name in read_data:
		#app.setMeter("progress_", 10, text=None)
			current_=current_+1
		
		else :
		 #app.setMeter("progress_", 0, text=None)
			print ""
		meter_=10.00
		meter_=current_*100
		meter_=meter_//counter_steps
		meter_=int(meter_)
		meter_Str="Analyzing: "+sample_name+" ("+str(meter_)+"%)"
		app.setMeter("progress_", meter_, text=meter_Str)
		time.sleep(5.0)

	app.setMeter("progress_", 100, text="Run finished")
	read_data=""
	if os.path.isfile(short_log):
		f = open(short_log, 'r')
		read_data = f.read()
		f.close()
	app.clearTextArea("t1")
	app.setTextArea("t1", read_data)
	read_data=""
	if os.path.isfile(out_log):
		f = open(out_log, 'r')
		read_data = f.read()
		f.close()
	info_="Detailed log file information"
	app.setTextArea("t1", info_+read_data)
	show_results()
	app.setTabbedFrameSelectedTab("TabbedFrame","Results")
	pass

def  pipeline(config_f, thread_):
	with open(config_f, 'r') as f:
		config = yaml.load(f)
		
	#Trim adapter and quality filtering
	def trim_galore(raw_reads_1, raw_reads_2):
		# Define variables from configuration file
		reads = str("--" + config['trim_galore']['reads_type'])
		quality = str(config["trim_galore"]["quality_threshold"])
		logging.info("Trimgalore is running...")
		#cutadapt=os.path.join(pipeline_path+'/cutadapt-1.12/cutadapt/cutadapt')
		#trimgalore_path=os.path.join(pipeline_path+"/trim_galore_v0.4.2/trim_galore")
		with open(file_out, 'a') as logf:
			if reads == '--paired':
				subprocess.call([trimgalore_path, "--paired", "-q", quality, "--suppress_warn", "-o", sample_directory, raw_reads_1, raw_reads_2], stdout=logf)
			else:
				subprocess.call([trimgalore_path, "-q", quality, "--suppress_warn", "-o", sample_directory, raw_reads_1], stdout=logf)
		logging.info("Reads trimming and quality filtering completed")
		os.chdir(sample_directory)
		#Correct filenames for spades input
		val_reads_1 = raw_reads_1.replace("_R1_001.fastq.gz", "_R1_001_val_1.fq.gz")
		val_reads_2 = raw_reads_2.replace("_R2_001.fastq.gz", "_R2_001_val_2.fq.gz")
		logging.info("TrimGalore finished "+sample_name)
		if spades_flag == 'False':
			spades(val_reads_1, val_reads_2)

	def call_script(args):
		subprocess.call(args)


	#Genome assembly
	def spades(spades_reads_pe1, spades_reads_pe2):
		# Define variables from configuration file
		kmer_size = str(config["spades"]["kmer"])
		os.chdir(sample_directory)
		spades_output = os.path.join(sample_directory, 'spades_assembly')
		if not os.path.exists(spades_output):
			os.makedirs(spades_output)
		logging.info("SPAdes genome assembler is running...")
		#spades_path=os.path.join(pipeline_path+'/SPAdes-3.9.1-Darwin/bin/spades.py')
		Sp_mode = config['spades']['Mode']
		with open(file_out, 'a') as logf:
			if Sp_mode == 'paired':
				subprocess.call([spades_path, "--careful", "-k", kmer_size, "-t", str(thread_), "-m", "24", "-1", spades_reads_pe1, "-2", spades_reads_pe2, "-o", spades_output], stdout=logf)
			elif Sp_mode == 'single':
				spades_reads_pe1 = spades_reads_pe1.replace("_R1_001_val_1.fq.gz","_R1_001_trimmed.fq.gz")
				subprocess.call([spades_path, "--careful", "-k", kmer_size, "-t", str(thread_), "-m", "24", "-s", spades_reads_pe1, "-o", spades_output], stdout=logf)
#			elif Sp_mode == 'pacbio':
#				spades_reads_pe1 = spades_reads_pe1.replace("_R1_001_val_1.fq.gz","_R1_001_trimmed.fq.gz")
#				subprocess.call([spades_path, "--careful", "-k", kmer_size, "-t", str(thread_), "-m", "24", "--pacbio", spades_reads_pe1, "-o", spades_output], stdout=logf)
#			elif Sp_mode == 'nanopore':
#				spades_reads_pe1 = spades_reads_pe1.replace("_R1_001_val_1.fq.gz","_R1_001_trimmed.fq.gz")
#				subprocess.call([spades_path, "--careful", "-k", kmer_size, "-t", str(thread_), "-m", "24", "--nanopore", spades_reads_pe1, "-o", spades_output], stdout=logf)
#			elif Sp_mode == 'nanopore':
#				spades_reads_pe1 = spades_reads_pe1.replace("_R1_001_val_1.fq.gz","_R1_001_trimmed.fq.gz")
#				subprocess.call([spades_path, "--careful", "-k", kmer_size, "-t", str(thread_), "-m", "24", "--sanger", spades_reads_pe1, "-o", spades_output], stdout=logf)'''
			elif Sp_mode == 'iontorrent':
				spades_reads_pe1 = spades_reads_pe1.replace("_R1_001_val_1.fq.gz","_R1_001_trimmed.fq.gz")
				subprocess.call([spades_path, "--careful", "-k", kmer_size, "-t", str(thread_), "-m", "24", "--iontorrent","--s", spades_reads_pe1, "-o", spades_output], stdout=logf)
			else:
				print "Error: unable to identify SPAdes mode"
			   
			
		os.chdir(spades_output)
		# Copy fasta assembly in the assemblies folder and rename it with sample id
		logging.info("SPAdes finished "+sample_name)
		if os.path.isfile('scaffolds.fasta'):
			assembly_fasta = 'scaffolds.fasta'
			shutil.copy(assembly_fasta, assemblies_directory)
			os.chdir(assemblies_directory)
			renamed_assembly = assembly_fasta.replace('scaffolds', sample_name)
			shutil.move(assembly_fasta, renamed_assembly)
			logging.info("Genome assembly completed "+sample_name)
			try:
				print ""
				parsnp_processors=int(thread_)
				thread1=threading.Thread(target=mlst_typing, args=(renamed_assembly,))
				thread2=threading.Thread(target=plasmids_finder, args=(renamed_assembly,))
				thread3=threading.Thread(target=resistance_finder, args=(renamed_assembly,))
				thread4=threading.Thread(target=virulence_finder, args=(renamed_assembly,))
				#thread5=threading.Thread(target=parSNP, args=(assemblies_directory,parsnp_processors,))
				thread6=threading.Thread(target=emmTyping, args=(renamed_assembly,))



				if mlst_typing_flag == 'False':
					#thread.start_new_thread(mlst_typing,(renamed_assembly,))
					thread1.start()
					parsnp_processors=parsnp_processors-1
				if plasmids_finder_flag == 'False':
					#thread.start_new_thread(plasmids_finder,(renamed_assembly,))
					thread2.start()
					parsnp_processors=parsnp_processors-1
				if resfinder_flag == 'False':
					#thread.start_new_thread(resistance_finder,(renamed_assembly,))
					thread3.start()
					parsnp_processors=parsnp_processors-1
				if virulencefinder_flag == 'False':
					#thread.start_new_thread(virulence_finder,(renamed_assembly,))
					thread4.start()
					parsnp_processors=parsnp_processors-1
				if emmTyping_flag == 'False':
					thread6.start()
					parsnp_processors=parsnp_processors-1
				#if parSNP_flag == 'False':
					#if parsnp_processors <1:
						#parsnp_processors=1
					#parsnp_processors=str(parsnp_processors)
					#thread5=threading.Thread(target=parSNP, args=(assemblies_directory,parsnp_processors,))
					#thread5.start()
				if mlst_typing_flag == 'False':
					thread1.join()
				if plasmids_finder_flag == 'False':
					thread2.join()
				if resfinder_flag == 'False':
					thread3.join()
				if virulencefinder_flag == 'False':
					thread4.join()
				#if parSNP_flag == 'False':
					#thread5.join()
				if emmTyping_flag == 'False':
					thread6.join()
 
			except:
				print "Error: unable to start thread after SPAdes"
			if prokka_flag == 'False':
				genome_annotation(renamed_assembly)
			elif Output_flag == 'False':
				Output()
			else:
				pass
		else:
			logging.info("**********************SPAdes assembly not found**********************")
			pass
###################
###################
###################
	def Post_assembled_tools(renamed_assembly):
		if os.path.isfile(renamed_assembly):
			try:
				print ""
				parsnp_processors=int(thread_)
				thread1=threading.Thread(target=mlst_typing, args=(renamed_assembly,))
				thread2=threading.Thread(target=plasmids_finder, args=(renamed_assembly,))
				thread3=threading.Thread(target=resistance_finder, args=(renamed_assembly,))
				thread4=threading.Thread(target=virulence_finder, args=(renamed_assembly,))
				#thread5=threading.Thread(target=parSNP, args=(assemblies_directory,parsnp_processors,))
				thread6=threading.Thread(target=emmTyping, args=(renamed_assembly,))



				if mlst_typing_flag == 'False':
					#thread.start_new_thread(mlst_typing,(renamed_assembly,))
					thread1.start()
					parsnp_processors=parsnp_processors-1
				if plasmids_finder_flag == 'False':
					#thread.start_new_thread(plasmids_finder,(renamed_assembly,))
					thread2.start()
					parsnp_processors=parsnp_processors-1
				if resfinder_flag == 'False':
					#thread.start_new_thread(resistance_finder,(renamed_assembly,))
					thread3.start()
					parsnp_processors=parsnp_processors-1
				if virulencefinder_flag == 'False':
					#thread.start_new_thread(virulence_finder,(renamed_assembly,))
					thread4.start()
					parsnp_processors=parsnp_processors-1
				if emmTyping_flag == 'False':
					thread6.start()
					parsnp_processors=parsnp_processors-1
				#if parSNP_flag == 'False':
					#if parsnp_processors <1:
						#parsnp_processors=1
					#parsnp_processors=str(parsnp_processors)
					#thread5=threading.Thread(target=parSNP, args=(assemblies_directory,parsnp_processors,))
					#thread5.start()
				if mlst_typing_flag == 'False':
					thread1.join()
				if plasmids_finder_flag == 'False':
					thread2.join()
				if resfinder_flag == 'False':
					thread3.join()
				if virulencefinder_flag == 'False':
					thread4.join()
				#if parSNP_flag == 'False':
					#thread5.join()
				if emmTyping_flag == 'False':
					thread6.join()
 
			except:
				print "Error: unable to start thread"
			if Output_flag == 'False':
			 	if Resfams_flag == 'False' or cardSearch_flag == 'False' or VirDBSearch_flag == 'False':
					print ""
				else:
					Output()
		else:
			print "File not found"
###################
###################
###################
	def mlst_typing(mlst_assembly):
		# Define variables from configuration file
		organism = str(config["mlst_typing"]["organism"])
		mlst_directory = os.path.join(sample_directory, 'mlst_typing')
		if not os.path.exists(mlst_directory):
			os.makedirs(mlst_directory)
		logging.info("MLST typing...")
		#mlst_path=os.path.join(pipeline_path+'/mlst/mlst.py')
		mlst_path=os.path.join(pipeline_path+'/mlst/mlst.pl')
		#mlst_DB_path=os.path.join(pipeline_path+'/mlst/mlst_db/')
		mlst_DB_path=os.path.join(pipeline_path+'/mlst/database/')
		#mlst_ncbi_path=os.path.join(pipeline_path+'/mlst/blast-2.2.26/')
		with open(file_out, 'a') as logf:
			subprocess.call([mlst_path, "-i", mlst_assembly, "-o", mlst_directory, "-d", mlst_DB_path, "-s", organism, "-b", mlst_ncbi_path], stdout=logf)
			#print "python3.5 %s -i %s -o %s -p  %s -s %s -x" % (mlst_path,mlst_assembly,mlst_directory,mlst_DB_path,organism)
			#subprocess.call(["python3.5",mlst_path,"-i",mlst_assembly,"-o", mlst_directory, "-p", mlst_DB_path,"-s", organism,"-x"], stdout=logf)

		logging.info("MLST finished "+sample_name)
		global mlst_output
		mlst_output=mlst_directory+'/results_tab.txt'
		#mlst_output=mlst_directory+'/results_tab.tsv'
		global mlst_samples
		mlst_samples=mlst_samples+mlst_output+","

	def plasmids_finder(plasmid_assembly):
		# Define variables from configuration file
		database = str(config["plasmids_finder"]["plasmids_database"])
		threshold = str(config["plasmids_finder"]["identity_threshold"])
		plasmids_directory = os.path.join(sample_directory, 'plasmids')
		if not os.path.exists(plasmids_directory):
			os.makedirs(plasmids_directory)
		logging.info("Finding plasmids...")
		plasmidfinder_path=os.path.join(pipeline_path+'/plasmidfinder/plasmidfinder.py')
		plasmidfinder_DB_path=os.path.join(pipeline_path+'/plasmidfinder/plasmidfinder_db')
		#plasmidfinder_ncbi_path=os.path.join(pipeline_path+'/mlst/blast-2.2.26/')
		with open(file_out, 'a') as logf:
		 #plasmidfinder.py -i ./test/test.fsa -o ./test/ -mp blastn -x -p ./database/ -q
			#print "python3.5 %s -i %s -o %s -p %s -mp blastn -x" % (plasmidfinder_path,plasmid_assembly,plasmids_directory,plasmidfinder_DB_path)
			subprocess.call(["python3.5",plasmidfinder_path, "-i", plasmid_assembly, "-o", plasmids_directory, "-p", plasmidfinder_DB_path, "-x","-mp", "blastn"], stdout=logf)
		logging.info("PlasmidFinder finished "+sample_name)
		global plasmidsFinder_output
		plasmidsFinder_output=plasmids_directory+"/results_tab.tsv"
		
		global Plasmidfinder_samples
		Plasmidfinder_samples=Plasmidfinder_samples+plasmidsFinder_output+","

	def resistance_finder(res_assembly):
		# Define variables from configuration file
		database = str(config["resfinder"]["resistance_database"])
		threshold = str(config["resfinder"]["identity_threshold"])
		minimum_overlap_length = str(config["resfinder"]["min_length"])
		res_directory = os.path.join(sample_directory, 'resistance_profile')
		if not os.path.exists(res_directory):
			os.makedirs(res_directory)
		logging.info("Finding antimicrobial resistance genes...")
		resistancefinder_path=os.path.join(pipeline_path+'/resfinder/resfinder.pl')
		resistancefinder_DB_path=os.path.join(pipeline_path+'/resfinder/resfinder_db')
		#resistancefinder_ncbi_path=os.path.join(pipeline_path+'/mlst/blast-2.2.26/')
		with open(file_out, 'a') as logf:
		 subprocess.call([resistancefinder_path, "-i", res_assembly, "-o", res_directory, "-d", resistancefinder_DB_path, "-k", threshold, "-l", minimum_overlap_length, "-a", database, "-b", resistancefinder_ncbi_path], stdout=logf)
		
		#subprocess.call([resistancefinder_path, "-i", res_assembly, "-o", res_directory, "-d", resistancefinder_DB_path, "-k", threshold, "-l", minimum_overlap_length, "-a", database, "-b", resistancefinder_ncbi_path], stdout=logf)
		logging.info("Resistance finished "+sample_name)
		global resistanceFinder_output
		resistanceFinder_output=res_directory+'/results_tab.txt'
		global resfinder_samples
		resfinder_samples=resfinder_samples+resistanceFinder_output+","

	def virulence_finder(vir_assembly):
		# Define variables from configuration file
		database = str(config["virulencefinder"]["virulence_database"])
		threshold = str(config["virulencefinder"]["identity_threshold"])#modify YAML
		vir_directory = os.path.join(sample_directory, 'virulence_profile')
		os.chdir(assemblies_directory)
		if not os.path.exists(vir_directory):
			os.makedirs(vir_directory)
		logging.info("Finding virulence...")
		#virulencefinder_path=os.path.join(pipeline_path+'/virulencefinder/virulencefinder.py')
		virulencefinder_path=os.path.join(pipeline_path+'/virulencefinder/virulencefinder.pl')
		#virulencefinder_DB_path=os.path.join(pipeline_path+'/virulencefinder/virulencefinder_db')
		virulencefinder_DB_path=os.path.join(pipeline_path+'/virulencefinder/database')
		virulencefinder_ncbi_path=os.path.join(pipeline_path+'/mlst/blast-2.2.26/')
		with open(file_out, 'a') as logf:
			#print "perl %s -i %s -o %s -d %s -k %s -s %s -b %s" % (virulencefinder_path,vir_assembly,vir_directory,virulencefinder_DB_path,threshold,database,virulencefinder_ncbi_path)
			subprocess.call(["perl",virulencefinder_path, "-i", vir_assembly, "-o", vir_directory, "-d", virulencefinder_DB_path, "-k", threshold, "-s", database, "-b", virulencefinder_ncbi_path], stdout=logf)
			#subprocess.call(["python3.5",virulencefinder_path, "-i", vir_assembly,"-o",vir_directory,"-mp", "blastn", "-x", "-p", virulencefinder_DB_path], stdout=logf)
			#print "python3.5 %s -i  %s -o  %s -mp blastn -x -p  %s" % (virulencefinder_path,vir_assembly,vir_directory,virulencefinder_DB_path)
		logging.info("VirulenceFinder finished "+sample_name)
		global virulenceFinder_output
		virulenceFinder_output=sample_directory+'/virulence_profile'+"/results_tab.txt"
		global virulence_samples
		virulence_samples=virulence_samples+virulenceFinder_output+","
	
	def emmTyping(emmTyping_assembly):
		# Define variables from configuration file
		emmTyping_directory = os.path.join(sample_directory, 'emmTyping_directory')
		#ann_directory = os.path.join(sample_directory, 'genome_annotation')
		#os.chdir(ann_directory)
		if not os.path.exists(emmTyping_directory):
			os.makedirs(emmTyping_directory)
		logging.info("emm Typing via blast...")
		results=os.path.join(emmTyping_directory + '/'+'results.txt')
		output=os.path.join(emmTyping_directory + '/'+'results_tab.txt')
		emmTyping_path=os.path.join(pipeline_path+'/emm_typing/Blast_emm.pl')
		blastall_ncbi_path=os.path.join(pipeline_path+'/mlst/blast-2.2.26/bin/blastall')
		mkblastdb_ncbi_path=os.path.join(pipeline_path+'/mlst/blast-2.2.26/bin/makeblastdb')
		emm_db_path=os.path.join(pipeline_path+'/emm_typing/emm_trimmed.fasta')
		with open(file_out, 'a') as logf:
			subprocess.call(["perl", emmTyping_path, emmTyping_assembly, output,emm_db_path,str(thread_), results,blastall_ncbi_path,mkblastdb_ncbi_path], stdout=logf)

		logging.info("emmTyping search finished "+sample_name)
		global emm_samples
		emm_samples=emm_samples+output+","

	def genome_annotation(annotation_assembly):
		ann_directory = os.path.join(sample_directory, 'genome_annotation')
		if not os.path.exists(ann_directory):
			os.makedirs(ann_directory)
		logging.info("Genome annotation with PROKKA...")
		prokka_path= str(config['prokka']['prokka_path'])
		prokka_path=prokka_path
		#barrnap_path=str(pipeline_path+'/barrnap-master/bin/')
		os.environ["PATH"] += os.pathsep + prokka_path
		os.environ["PATH"] += os.pathsep + barrnap_path
		with open(file_out, 'a') as logf:
			subprocess.call([prokka_path, "--quiet", "--cpus", str(thread_), "--kingdom", "Bacteria", "--addgenes", "--outdir", ann_directory, "--strain", sample_name, "--force", "--centre", "C", "--locustag", "L","--prefix","results",annotation_assembly], stdout=logf)
		os.chdir(ann_directory)
		# Rename gbk file with sample name
		gbk_file = glob.glob("*.gbk")
		gff_file = glob.glob("*.gff")
		faa_file = glob.glob("*.faa")
		renamed_gbk = gbk_file[0].replace(gbk_file[0], sample_name + '.gbk')
		renamed_gff = gff_file[0].replace(gff_file[0], sample_name + '.gff')

		shutil.move(gbk_file[0],renamed_gbk)
		shutil.move(gff_file[0],renamed_gff)
		Restfams_faa=os.path.join(ann_directory, '/',faa_file[0])
		logging.info("Annotation finished "+sample_name)
			##if Resfams_flag == 'False':
			##Resfams(faa_file[0])
			##if cardSearch_flag == 'False':
			##cardSearch(faa_file[0])
			##if VirDBSearch_flag == 'False':
			##VirDBSearch(faa_file[0])


		try:
			thread6=threading.Thread(target=Resfams, args=(faa_file[0],))
			thread7=threading.Thread(target=cardSearch, args=(faa_file[0],))
			thread8=threading.Thread(target=VirDBSearch, args=(faa_file[0],))

			if Resfams_flag == 'False':
				thread6.start()
				#thread.start_new_thread(Resfams,(faa_file[0],))
			if cardSearch_flag == 'False':
				thread7.start()
				#thread.start_new_thread(cardSearch,(faa_file[0],))
			if VirDBSearch_flag == 'False':
				thread8.start()
				#thread.start_new_thread(VirDBSearch,(faa_file[0],))
			if Resfams_flag == 'False':
				thread6.join()
			if cardSearch_flag == 'False':
				thread7.join()
			if VirDBSearch_flag == 'False':
				thread8.join()

		except:
			print "Error: unable to start thread after Prokka"

		if Output_flag == 'False':
			Output()

	def post_annotation_analysis(annotation_assembly_faa):
		if os.path.isfile(annotation_assembly_faa):
			try:
				thread6=threading.Thread(target=Resfams, args=(annotation_assembly_faa,))
				thread7=threading.Thread(target=cardSearch, args=(annotation_assembly_faa,))
				thread8=threading.Thread(target=VirDBSearch, args=(annotation_assembly_faa,))

				if Resfams_flag == 'False':
					thread6.start()
					#thread.start_new_thread(Resfams,(faa_file[0],))
				if cardSearch_flag == 'False':
					thread7.start()
					#thread.start_new_thread(cardSearch,(faa_file[0],))
				if VirDBSearch_flag == 'False':
					thread8.start()
					#thread.start_new_thread(VirDBSearch,(faa_file[0],))
				if Resfams_flag == 'False':
					thread6.join()
				if cardSearch_flag == 'False':
					thread7.join()
				if VirDBSearch_flag == 'False':
					thread8.join()

			except:
				print "Error: unable to start thread after Prokka"

		if Output_flag == 'False':
			Output()


	def Resfams(Resfams_faa):
		# Define variables from configuration file
		Resfams_directory = os.path.join(sample_directory, 'Resfams_directory')
		#ann_directory = os.path.join(sample_directory, 'genome_annotation')
		#os.chdir(ann_directory)
		if not os.path.exists(Resfams_directory):
			os.makedirs(Resfams_directory)
		logging.info("denovo annotation Resfams...")
		results=os.path.join(Resfams_directory + '/'+'results.txt')
		output=os.path.join(Resfams_directory + '/'+'output.txt')
		#hmmer_path=os.path.join(pipeline_path+'/hmmer-3.1b2-macosx-intel/binaries/hmmscan')
		resfams_path=os.path.join(pipeline_path+'/ResFam/Resfams.hmm')
		with open(file_out, 'a') as logf:
			subprocess.call([hmmer_path, "--cut_ga", "--tblout", results,"-o",output, resfams_path,Resfams_faa], stdout=logf)
		logging.info("Resfams finished "+sample_name)
		global resfams_output
		resfams_output=Resfams_directory+'/results.txt'
		global resfams_samples
		resfams_samples=resfams_samples+resfams_output+","

	def cardSearch(card_faa):
		# Define variables from configuration file
		cardSearch_directory = os.path.join(sample_directory, 'CARDsearch_directory')
		#ann_directory = os.path.join(sample_directory, 'genome_annotation')
		#os.chdir(ann_directory)
		if not os.path.exists(cardSearch_directory):
			os.makedirs(cardSearch_directory)
		logging.info("CARD search via blast...")
		results=os.path.join(cardSearch_directory + '/'+'results.txt')
		output=os.path.join(cardSearch_directory + '/'+'output.txt')
		CARDSearch_path=os.path.join(pipeline_path+'/cardSearch/Blast_card.pl')
		blastp_ncbi_path=os.path.join(pipeline_path+'/mlst/blast-2.2.26/bin/blastp')
		card_db_path=os.path.join(pipeline_path+'/cardSearch/card_protien.fasta')
		with open(file_out, 'a') as logf:
			subprocess.call(["perl", CARDSearch_path, card_faa, output,card_db_path,str(thread_), results,blastp_ncbi_path], stdout=logf)
		#pipe = subprocess.Popen(["perl", CARDSearch_path,card_faa,output,card_db_path,str(thread_),results], stdout=subprocess.PIPE)
		logging.info("CARD search finished "+sample_name)
		global card_output
		card_output=cardSearch_directory+'/results.txt'
		global card_samples
		card_samples=card_samples+card_output+","

	def VirDBSearch(VirDB_faa):
		# Define variables from configuration file
		VirDBSearch_directory = os.path.join(sample_directory, 'VirDBSearch_directory')
		#ann_directory = os.path.join(sample_directory, 'genome_annotation')
		#os.chdir(ann_directory)
		if not os.path.exists(VirDBSearch_directory):
			os.makedirs(VirDBSearch_directory)
		logging.info("VirDB search via blast...")
		results=os.path.join(VirDBSearch_directory + '/'+'results.txt')
		output=os.path.join(VirDBSearch_directory + '/'+'output.txt')
		VIRDBSearch_path=os.path.join(pipeline_path+'/VirDB/Blast_VirDB.pl')
		blastp_ncbi_path=os.path.join(pipeline_path+'/mlst/blast-2.2.26/bin/blastp')
		vir_db_path=os.path.join(pipeline_path+'/VirDB/VFDB_setA_pro.fas')
		pipe = subprocess.Popen(["perl", VIRDBSearch_path,VirDB_faa,output,vir_db_path,str(thread_),results,blastp_ncbi_path], stdout=subprocess.PIPE)
		logging.info("VirDB search finished "+sample_name)
		global virDB_output
		virDB_output=VirDBSearch_directory+'/results.txt'
		global virdb_samples
		virdb_samples=virdb_samples+virDB_output+","

	def parSNP(parSNP_fsa,parsnp_pros):
		# Define variables from configuration file
		parSNP_directory = os.path.join(output_directory, 'parSNP_directory')
		##parSNP_directory = os.path.join(sample_directory, 'parSNP_directory')
		#ann_directory = os.path.join(sample_directory, 'genome_annotation')
		#os.chdir(ann_directory)
		if not os.path.exists(parSNP_directory):
			os.makedirs(parSNP_directory)
		logging.info("SNP analysis via parSNP ...")
		results=os.path.join(parSNP_directory + '/'+'results.txt')
		output=os.path.join(parSNP_directory + '/'+'output.vcf')
		global parSNP_tree
		parSNP_tree=os.path.join(parSNP_directory + '/'+'parsnp.tree')
		#parSNP_path=os.path.join(pipeline_path+'/harvesttools_OSX64_v1_3/harvesttools')
		#parSNP_path=os.path.join(pipeline_path+'/parsnp_OSX64_v1_2/parsnp')
		parSNP_reference= str(config['parSNP']['parSNP_reference'])
		parSNP_reference_fsa= str(config['parSNP']['parSNP_reference_fsa'])

		with open(file_out, 'a') as logf:
		 
			if parSNP_reference == "None" :
				subprocess.call([parSNP_path, "-r", parSNP_reference_fsa, "-d", parSNP_fsa, "-o", parSNP_directory,"-p",parsnp_pros, "-c"], stdout=logf)
		
			else:
				subprocess.call([parSNP_path, "-g", parSNP_reference, "-d", parSNP_fsa, "-o", parSNP_directory,"-p",parsnp_pros, "-c"], stdout=logf)
			#print "%s -g %s -d %s -o %s -p %s" % (parSNP_path,parSNP_reference,parSNP_fsa,parSNP_directory,parsnp_pros)
		logging.info("parSNP finished "+sample_name)

	def Output():

		logging.info("Summarizing outputs...")
		tools=""
		if mlst_typing_flag == 'False':
			tools=tools+"0,"
		if plasmids_finder_flag == 'False':
			tools=tools+"1,"
		if resfinder_flag == 'False':
			tools=tools+"2,"
		if virulencefinder_flag == 'False':
			tools=tools+"3,"
		if prokka_flag == 'False':
			tools=tools+"4,"
		if Resfams_flag == 'False':
			tools=tools+"5,"
		if cardSearch_flag == 'False':
			tools=tools+"6,"
		if VirDBSearch_flag == 'False':
			tools=tools+"7,"
			#if emmTyping_flag == 'False':
			#tools=tools+"8,"
	
				 #			if parSNP_flag == 'False':
				#				tools=tools+"8,"
		if tools == "":
			print ""
		else:
			tools=tools[:-1]
			Output_batch_path=os.path.join(pipeline_path+'/Output_batch.pl')
			pipe = subprocess.Popen(["perl", Output_batch_path,output_directory,sample_name, tools], stdout=subprocess.PIPE)
			logging.info("summarization finished "+sample_name)

#os.chdir(input_directory)
		logging.info("-----Finish processing sample " + sample_name + "-----")


	## Create directories
	output_directory = config["directories"]["output"]
	if not os.path.exists(output_directory):
		os.makedirs(output_directory)

	input_directory = config["directories"]["reads"]

	assemblies_directory = os.path.join(output_directory, 'genome_assemblies')
	if not os.path.exists(assemblies_directory):
		os.makedirs(assemblies_directory)

	summary_directory = os.path.join(output_directory, 'Summary')
	if not os.path.exists(summary_directory):
		os.makedirs(summary_directory)

	os.chdir(output_directory)

	logfile ="pipeline.log"

	if logging.root:
		del logging.root.handlers[:]

	logging.basicConfig(
		filename=logfile,
		level=logging.DEBUG,
		filemode='w',
		format='%(asctime)s %(message)s',
		datefmt='%m/%d/%Y %H:%M:%S')

	class StreamToLogger(object):
	   """
	   Fake file-like stream object that redirects writes to a logger instance.
	   """
	   def __init__(self, logger, log_level=logging.INFO):
		  self.logger = logger
		  self.log_level = log_level
		  self.linebuf = ''
	 
	   def write(self, buf):
		  for line in buf.rstrip().splitlines():
			 self.logger.log(self.log_level, line.rstrip())

	stdout_logger = logging.getLogger('STDOUT')
	sl = StreamToLogger(stdout_logger, logging.INFO)
	sys.stdout = sl
	 
	stderr_logger = logging.getLogger('STDERR')
	sl = StreamToLogger(stderr_logger, logging.ERROR)
	sys.stderr = sl

	#STDOUT
	file_out=os.path.join(output_directory,'log.txt')

	#Store in a variable the text file with filenames of reads already processed
	files_processed = open('files_processed.txt', 'a+')
	if os.path.getsize('files_processed.txt') == 0:
		files_processed.write('----- List of processed files -----' + '\n')
	files_processed.close()

	os.chdir(input_directory)

	#deactivated tools
	if 'deactivate' in config['trim_galore']:
		trim_galore_flag= str(config['trim_galore']['deactivate'])
	else:
		trim_galore_flag= 'False'
	if 'deactivate' in config['spades']:
		spades_flag= str(config['spades']['deactivate'])
	else:
		spades_flag= 'False'
	if 'deactivate' in config['mlst_typing']:
		mlst_typing_flag= str(config['mlst_typing']['deactivate'])
	else:
		mlst_typing_flag= 'False'
	if 'deactivate' in config['plasmids_finder']:
		plasmids_finder_flag= str(config['plasmids_finder']['deactivate'])
	else:
		plasmids_finder_flag= 'False'
	if 'deactivate' in config['cardSearch']:
		cardSearch_flag= str(config['cardSearch']['deactivate'])
	else:
		cardSearch_flag= 'False'
	if 'deactivate' in config['VirDBSearch']:
		VirDBSearch_flag= str(config['VirDBSearch']['deactivate'])
	else:
		VirDBSearch_flag= 'False'
	if 'deactivate' in config['parSNP']:
		parSNP_flag= str(config['parSNP']['deactivate'])
	else:
		parSNP_flag= 'False'
	if 'deactivate' in config['emmTyping']:
		emmTyping_flag= str(config['emmTyping']['deactivate'])
	else:
		emmTyping_flag= 'False'
	if 'deactivate' in config['resfinder']:
		resfinder_flag= str(config['resfinder']['deactivate'])
	else:
		resfinder_flag= 'False'
	if 'deactivate' in config['prokka']:
		prokka_flag= str(config['prokka']['deactivate'])
	else:
		prokka_flag= 'False'
	if 'deactivate' in config['Resfams']:
		Resfams_flag= str(config['Resfams']['deactivate'])
	else:
		Resfams_flag= 'False'
	if 'deactivate' in config['virulencefinder']:
		virulencefinder_flag= str(config['virulencefinder']['deactivate'])
	else:
		virulencefinder_flag= 'False'
	if 'deactivate' in config['Output']:
		Output_flag= str(config['Output']['deactivate'])
	else:
		Output_flag= 'False'

	#printing tools deactivated
	if trim_galore_flag == 'True':
		logging.info("deactivating trim galore as specified")
	if spades_flag == 'True':
		logging.info("deactivating SPAdes as specified")
	if mlst_typing_flag == 'True':
		logging.info("deactivating MLST typing as specified")
	if plasmids_finder_flag == 'True':
		logging.info("deactivating PlasmidsFinder  as specified")
	if resfinder_flag == 'True':
		logging.info("deactivating ResFinder as specified")
	if cardSearch_flag == 'True':
		logging.info("deactivating CARD Search as specified")
	if VirDBSearch_flag == 'True':
		logging.info("deactivating VirDB Search as specified")
	if parSNP_flag == 'True':
		logging.info("deactivating parSNP as specified")
	if emmTyping_flag == 'True':
		logging.info("deactivating M-Typing as specified")
	if prokka_flag == 'True':
		logging.info("deactivating Prokka as specified")
	if Resfams_flag == 'True':
		logging.info("deactivating Resfams as specified")
	if virulencefinder_flag == 'True':
		logging.info("deactivating VirulenceFinder as specified")
	if Output_flag == 'True':
		logging.info("deactivating Output summarization as specified")

	#Check for paired end reads to process in the input directory
	global sample_name
	if trim_galore_flag == 'False':
		tracker=0
		list = [f for f in glob.glob("*_R1_*") if "fastq.gz" in f or "fq.gz" in f ]
		for illumina_reads_1 in list:
			logging.info("Reading files ending with _R1_001.fastq.gz for forward and _R1_001.fastq.gz for reverse within the input directory")
			illumina_reads_2 = illumina_reads_1.replace("_R1_", "_R2_")
			sample_name = illumina_reads_1[:-9]
			os.chdir(output_directory)
			files_processed = open('files_processed.txt', 'a+')
			row = [x.strip('\n') for x in files_processed.readlines()]
			if str(illumina_reads_1) in row:
				logging.info(sample_name + " already processed")
				continue
			else:
				sample_directory = os.path.join(output_directory, illumina_reads_1[:-9])
				if not os.path.exists(sample_directory):
					os.makedirs(sample_directory)
				os.chdir(input_directory)
				logging.info("-----Start processing sample " + sample_name + "-----")
				trim_galore(illumina_reads_1, illumina_reads_2)
				files_processed.write(illumina_reads_1 + '\n')
				files_processed.write(illumina_reads_2 + '\n')
				files_processed.close()
			tracker=1
		if parSNP_flag == 'False' and tracker == 1:
			logging.info("Reading assembled files ending .fasta within the genome_assemblies directory")
			logging.info("-----Start processing assembled genomes -----")
			Processors_=str(app.getEntry('Processors'))
			parSNP_input=os.path.join(output_directory, 'genome_assemblies')
			parSNP(parSNP_input,Processors_)
		elif tracker == 0:
			print "No input files for TrimGalore"
		else:
			print ""
 
	#Start from trim galore output
	elif spades_flag == 'False':
		list = [f for f in glob.glob("*_R1_*") if "fastq.gz" in f or "fq.gz" in f ]
		for illumina_reads_1 in list:
			logging.info("Reading files ending with _R1_001_val_1.fq.gz for forward and _R2_001_val_2.fq.gz for reverse within the input directory")
			illumina_reads_2 = illumina_reads_1.replace("_R1_", "_R2_")
			sample_name = illumina_reads_1[:-12]
			os.chdir(output_directory)
			files_processed = open('files_processed.txt', 'a+')
			row = [x.strip('\n') for x in files_processed.readlines()]
			if str(illumina_reads_1) in row:
				logging.info(sample_name + " already processed")
				continue
			else:
				sample_directory = os.path.join(output_directory, illumina_reads_1[:-12])
				if not os.path.exists(sample_directory):
					os.makedirs(sample_directory)
				os.chdir(input_directory)
				logging.info("-----Start processing sample " + sample_name + "-----")
				illumina_reads_1=os.path.join(input_directory,illumina_reads_1)
				illumina_reads_2=os.path.join(input_directory,illumina_reads_2)
				spades(illumina_reads_1, illumina_reads_2)
				files_processed.write(illumina_reads_1 + '\n')
				files_processed.write(illumina_reads_2 + '\n')
				files_processed.close()
		if parSNP_flag == 'False':
			logging.info("Reading assembled files ending .fasta within the genome_assemblies directory")
			logging.info("-----Start processing assembled genomes -----")
			Processors_=str(app.getEntry('Processors'))
			parSNP_input=os.path.join(output_directory, 'genome_assemblies')
			parSNP(parSNP_input,Processors_)

	else:


		#Start from Spades output 701_S8_L001_R1_001.fasta
		if mlst_typing_flag == 'False':
		 	os.chdir(input_directory)
			for illumina_reads_1 in glob.glob("*.fasta"):
				logging.info("Reading assembled files ending .fasta within the input directory")
				sample_name = illumina_reads_1[:-6]
				os.chdir(output_directory)
				files_processed = open('files_processed.txt', 'a+')
				row = [x.strip('\n') for x in files_processed.readlines()]
				if str(illumina_reads_1) in row:
					logging.info(sample_name + " already processed")
					continue
				else:
					sample_directory = os.path.join(output_directory, illumina_reads_1[:-6])
					if not os.path.exists(sample_directory):
						os.makedirs(sample_directory)
					os.chdir(input_directory)
					logging.info("-----Start processing sample " + sample_name + "-----")
					illumina_reads_1=os.path.join(input_directory,illumina_reads_1)
					#mlst_typing(illumina_reads_1)
					Post_assembled_tools(illumina_reads_1)
					files_processed.write(illumina_reads_1 + '\n')
					files_processed.close()

		#Start from Spades output 701_S8_L001_R1_001.fasta
		elif plasmids_finder_flag == 'False':
			os.chdir(input_directory)
			for illumina_reads_1 in glob.glob("*.fasta"):
				logging.info("Reading assembled files ending .fasta within the input directory")
				sample_name = illumina_reads_1[:-6]
				os.chdir(output_directory)
				files_processed = open('files_processed.txt', 'a+')
				row = [x.strip('\n') for x in files_processed.readlines()]
				if str(illumina_reads_1) in row:
					logging.info(sample_name + " already processed")
					continue
				else:
					sample_directory = os.path.join(output_directory, illumina_reads_1[:-6])
					if not os.path.exists(sample_directory):
						os.makedirs(sample_directory)
					os.chdir(input_directory)
					logging.info("-----Start processing sample " + sample_name + "-----")
					illumina_reads_1=os.path.join(input_directory,illumina_reads_1)
					#plasmids_finder(illumina_reads_1)
					Post_assembled_tools(illumina_reads_1)
					files_processed.write(illumina_reads_1 + '\n')
					files_processed.close()

		elif resfinder_flag == 'False':
		 	os.chdir(input_directory)
			for illumina_reads_1 in glob.glob("*.fasta"):
				logging.info("Reading assembled files ending .fasta within the input directory")
				sample_name = illumina_reads_1[:-6]
				os.chdir(output_directory)
				files_processed = open('files_processed.txt', 'a+')
				row = [x.strip('\n') for x in files_processed.readlines()]
				if str(illumina_reads_1) in row:
					logging.info(sample_name + " already processed")
					continue
				else:
					sample_directory = os.path.join(output_directory, illumina_reads_1[:-6])
					if not os.path.exists(sample_directory):
						os.makedirs(sample_directory)
					os.chdir(input_directory)
					logging.info("-----Start processing sample " + sample_name + "-----")
					illumina_reads_1=os.path.join(input_directory,illumina_reads_1)
					#resistance_finder(illumina_reads_1)
					Post_assembled_tools(illumina_reads_1)
					files_processed.write(illumina_reads_1 + '\n')
					files_processed.close()

		elif virulencefinder_flag == 'False':
		 	os.chdir(input_directory)
			for illumina_reads_1 in glob.glob("*.fasta"):
				logging.info("Reading assembled files ending .fasta within the input directory")
				sample_name = illumina_reads_1[:-6]
				os.chdir(output_directory)
				files_processed = open('files_processed.txt', 'a+')
				row = [x.strip('\n') for x in files_processed.readlines()]
				if str(illumina_reads_1) in row:
					logging.info(sample_name + " already processed")
					continue
				else:
					sample_directory = os.path.join(output_directory, illumina_reads_1[:-6])
					if not os.path.exists(sample_directory):
						os.makedirs(sample_directory)
					os.chdir(input_directory)
					logging.info("-----Start processing sample " + sample_name + "-----")
					illumina_reads_1=os.path.join(input_directory,illumina_reads_1)
					#virulence_finder(illumina_reads_1)
					Post_assembled_tools(illumina_reads_1)
					files_processed.write(illumina_reads_1 + '\n')
					files_processed.close()

		elif emmTyping_flag == 'False':
		 	os.chdir(input_directory)
			for illumina_reads_1 in glob.glob("*.fasta"):
				logging.info("Reading assembled files ending .fasta within the input directory")
				sample_name = illumina_reads_1[:-6]
				os.chdir(output_directory)
				files_processed = open('files_processed.txt', 'a+')
				row = [x.strip('\n') for x in files_processed.readlines()]
				if str(illumina_reads_1) in row:
					logging.info(sample_name + " already processed")
					continue
				else:
					sample_directory = os.path.join(output_directory, illumina_reads_1[:-6])
					if not os.path.exists(sample_directory):
						os.makedirs(sample_directory)
					os.chdir(input_directory)
					logging.info("-----Start processing sample " + sample_name + "-----")
					illumina_reads_1=os.path.join(input_directory,illumina_reads_1)
					#emmTyping(illumina_reads_1)
					Post_assembled_tools(illumina_reads_1)
					files_processed.write(illumina_reads_1 + '\n')
					files_processed.close()
		else:
			print ""

		if parSNP_flag == 'False':
			os.chdir(input_directory)
			if os.path.isdir(input_directory):
				sample_name = os.path.basename(input_directory)
			#for illumina_reads_1 in glob.glob("*.fasta"):
			logging.info("Reading assembled files ending .fasta within the input directory")
			#os.chdir(output_directory)
			#sample_directory = os.path.join(output_directory, sample_name)
			#if not os.path.exists(sample_directory):
			#	os.makedirs(sample_directory)
			os.chdir(input_directory)
			logging.info("-----Start processing sample " + sample_name + "-----")
			Processors_=str(app.getEntry('Processors'))
			parSNP(input_directory,Processors_)


		if prokka_flag == 'False':
		 	os.chdir(input_directory)
			for illumina_reads_1 in glob.glob("*.fasta"):
				logging.info("Reading assembled files ending .fasta within the input directory")
				sample_name = illumina_reads_1[:-6]
				os.chdir(output_directory)
				files_processed = open('files_processed.txt', 'a+')
				row = [x.strip('\n') for x in files_processed.readlines()]
				if str(illumina_reads_1) in row:
					logging.info(sample_name + " already processed")
					continue
				else:
					sample_directory = os.path.join(output_directory, illumina_reads_1[:-6])
					if not os.path.exists(sample_directory):
						os.makedirs(sample_directory)
					os.chdir(input_directory)
					logging.info("-----Start processing sample " + sample_name + "-----")
					illumina_reads_1=os.path.join(input_directory,illumina_reads_1)
					genome_annotation(illumina_reads_1)
					files_processed.write(illumina_reads_1 + '\n')
					files_processed.close()

		else:
			if Resfams_flag == 'False':
		 		os.chdir(input_directory)
				for illumina_reads_1 in glob.glob("*.faa"):
					logging.info("Reading protien files ending .faa within the input directory")
					sample_name = illumina_reads_1[:-4]
					os.chdir(output_directory)
					files_processed = open('files_processed.txt', 'a+')
					row = [x.strip('\n') for x in files_processed.readlines()]
					if str(illumina_reads_1) in row:
						logging.info(sample_name + " already processed")
						continue
					else:
						sample_directory = os.path.join(output_directory, illumina_reads_1[:-4])
						if not os.path.exists(sample_directory):
							os.makedirs(sample_directory)
						os.chdir(input_directory)
						logging.info("-----Start processing sample " + sample_name + "-----")
						illumina_reads_1=os.path.join(input_directory,illumina_reads_1)
						post_annotation_analysis(illumina_reads_1)
						#Resfams(illumina_reads_1)
						files_processed.write(illumina_reads_1 + '\n')
						files_processed.close()

			elif cardSearch_flag == 'False':
			 	os.chdir(input_directory)
				for illumina_reads_1 in glob.glob("*.faa"):
					logging.info("Reading protien files ending .faa within the input directory")
					sample_name = illumina_reads_1[:-4]
					os.chdir(output_directory)
					files_processed = open('files_processed.txt', 'a+')
					row = [x.strip('\n') for x in files_processed.readlines()]
					if str(illumina_reads_1) in row:
						logging.info(sample_name + " already processed")
						continue
					else:
						sample_directory = os.path.join(output_directory, illumina_reads_1[:-4])
						if not os.path.exists(sample_directory):
							os.makedirs(sample_directory)
						os.chdir(input_directory)
						logging.info("-----Start processing sample " + sample_name + "-----")
						illumina_reads_1=os.path.join(input_directory,illumina_reads_1)
						post_annotation_analysis(illumina_reads_1)
						#cardSearch(illumina_reads_1)
						files_processed.write(illumina_reads_1 + '\n')
						files_processed.close()

			elif VirDBSearch_flag == 'False':
			 	os.chdir(input_directory)
				for illumina_reads_1 in glob.glob("*.faa"):
					logging.info("Reading protien files ending .faa within the input directory")
					sample_name = illumina_reads_1[:-4]
					os.chdir(output_directory)
					files_processed = open('files_processed.txt', 'a+')
					row = [x.strip('\n') for x in files_processed.readlines()]
					if str(illumina_reads_1) in row:
						logging.info(sample_name + " already processed")
						continue
					else:
						sample_directory = os.path.join(output_directory, illumina_reads_1[:-4])
						if not os.path.exists(sample_directory):
							os.makedirs(sample_directory)
						os.chdir(input_directory)
						logging.info("-----Start processing sample " + sample_name + "-----")
						illumina_reads_1=os.path.join(input_directory,illumina_reads_1)
						#VirDBSearch(illumina_reads_1)
						post_annotation_analysis(illumina_reads_1)
						files_processed.write(illumina_reads_1 + '\n')
						files_processed.close()


			else:
 				print ""
				 


	if (Output_flag == 'True' and Resfams_flag == 'True' and prokka_flag == 'True' and virulencefinder_flag == 'True' and emmTyping_flag == 'True' and parSNP_flag == 'True' and VirDBSearch_flag == 'True' and cardSearch_flag == 'True' and resfinder_flag == 'True' and plasmids_finder_flag == 'True' and mlst_typing_flag == 'True' and spades_flag == 'True' and trim_galore_flag == 'True'):
		logging.info("Nothing to do here!!")
	logging.info("----------------------Script Completed---------------------")
	global Quit_Flag
	Quit_Flag=1

'''
		if Output_flag == 'False':
		 	os.chdir(input_directory)
			for illumina_reads_1 in glob.glob("*.fasta"):
				sample_name = illumina_reads_1[:-6]
				os.chdir(output_directory)
				files_processed = open('files_processed.txt', 'a+')
				row = [x.strip('\n') for x in files_processed.readlines()]
				if str(illumina_reads_1) in row:
					logging.info(sample_name + " already processed")
					continue
				else:
					sample_directory = os.path.join(output_directory, illumina_reads_1[:-6])
					if not os.path.exists(sample_directory):
						os.makedirs(sample_directory)
					os.chdir(input_directory)
					logging.info("-----Start processing sample " + sample_name + "-----")
					illumina_reads_1=os.path.join(input_directory,illumina_reads_1)
					Output()
					files_processed.write(illumina_reads_1 + '\n')
					files_processed.close()
'''


def press(btn):
	#Resistance_detail()
	global Plasmidfinder_samples
	Plasmidfinder_samples=""
	global resfinder_samples
	resfinder_samples=""
	global virulence_samples
	virulence_samples=""
	global mlst_samples
	mlst_samples=""
	global card_samples
	card_samples=""
	global virdb_samples
	virdb_samples=""
	global resfams_samples
	resfams_samples=""
	global output_path
	app.setMeter("progress_", 0, text=None)
	app.clearTextArea("t1")
	if btn=="Cancel":
		app.stop()
	else:
#		Cutadapt_flag=app.getCheckBox("Install Cutadapt (first time run)")
#		if Cutadapt_flag:
#			os.chdir(cutadapt_path)
#			print "python"
#			call(["python","setup.py","install","--user"])
		if input_path == "":
			print "Please select a valid path for your inputs\n"
			app.errorBox("Error", "Please select a valid path for your inputs")
			return
		global output_path
		if output_path == "":
			print "Please select a valid path for your outputs\n"
			app.errorBox("Error", "Please select a valid path for your outputs")
			return
		TG_switch=app.getRadioButton("Trim Galore")
		TG_paired=app.getOptionBox("Reads")
		TG_qual_thred=app.getEntry('TG_qual')
		if(TG_switch == "yes" and TG_qual_thred == ""):
			print "Please select a valid quality cutoff or you can  switch off Trim Galore\n"
			app.errorBox("Error", "Please select a valid quality cutoff or you can  switch off Trim Galore")
			return

		spades_switch=app.getRadioButton("spades")
		SPAdes_mode=app.getOptionBox("Mode")
		spades_kmer=app.getEntry('spades_kmer_')
		if(spades_switch == "yes" and spades_kmer == ""):
			print "Please select a valid cutoff or you can  switch off Spades\n"
			app.errorBox("Error", "Please select a valid cutoff or you can  switch off Spades")
			return
   
		global mlst_switch
		mlst_switch=app.getRadioButton('mlst_')
		mlst_organism=app.getOptionBox("Species")

		global plasmids_finder_switch
		plasmids_finder_switch=app.getRadioButton('plasmids_finder_')
		plasmids_database_type="NA"#app.getOptionBox('Database ')
		plasmids_finder_threshold=app.getEntry('plasmids_finder_thre')
		if(plasmids_finder_switch == "yes" and plasmids_finder_threshold == ""):
			 print "Please select a valid cutoff or you can  switch off PlasmidsFinder\n"
			 app.errorBox("Error", "Please select a valid cutoff or you can  switch off PlasmidsFinder")
			 return
		global resfinder_switch
		resfinder_switch=app.getRadioButton('resfinder_')
		resfinder_threshold=app.getEntry('resfinder_thre')

		global emmTyping_switch
		emmTyping_switch=app.getRadioButton('emm_')

		global parSNP_switch
		parSNP_switch=app.getRadioButton('parSNP_')

		global cardSearch_switch
		cardSearch_switch=app.getRadioButton('card_')

		global VirDBSearch_switch
		VirDBSearch_switch=app.getRadioButton('VirDB_')


		resistance_database_input=''
		if(app.getProperty("Antibiotics","aminoglycoside")):
			resistance_database_input=resistance_database_input+'aminoglycoside,'
		if(app.getProperty("Antibiotics","beta-lactam")):
			resistance_database_input=resistance_database_input+'beta-lactam,'
		if(app.getProperty("Antibiotics","colistin")):
			resistance_database_input=resistance_database_input+'colistin,'
		if(app.getProperty("Antibiotics","fosfomycin")):
			resistance_database_input=resistance_database_input+'fosfomycin,'
		if(app.getProperty("Antibiotics","fusidicacid")):
			resistance_database_input=resistance_database_input+'fusidicacid,'
		if(app.getProperty("Antibiotics","macrolide")):
			resistance_database_input=resistance_database_input+'macrolide,'
		if(app.getProperty("Antibiotics","nitroimidazole")):
			resistance_database_input=resistance_database_input+'nitroimidazole,'
		if(app.getProperty("Antibiotics","oxazolidinone")):
			resistance_database_input=resistance_database_input+'oxazolidinone,'
		if(app.getProperty("Antibiotics","phenicol")):
			resistance_database_input=resistance_database_input+'phenicol,'
		if(app.getProperty("Antibiotics","quinolone")):
			resistance_database_input=resistance_database_input+'quinolone,'
		if(app.getProperty("Antibiotics","rifampicin")):
			resistance_database_input=resistance_database_input+'rifampicin,'
		if(app.getProperty("Antibiotics","sulphonamide")):
			resistance_database_input=resistance_database_input+'sulphonamide,'
		if(app.getProperty("Antibiotics","tetracycline")):
			resistance_database_input=resistance_database_input+'tetracycline,'
		if(app.getProperty("Antibiotics","trimethoprim")):
			resistance_database_input=resistance_database_input+'trimethoprim,'
		if(app.getProperty("Antibiotics","glycopeptide")):
			resistance_database_input=resistance_database_input+'glycopeptide,'

		resistance_database_input = resistance_database_input[:-1]

		resfinder_min_length=str(0.60)#app.getEntry('resfinder_min_lngth')

		if(resfinder_switch == "yes" and (resistance_database_input =="" or resfinder_min_length=="" or resfinder_threshold == "")):
			print "Please select (at least) one database or insert a valid cutoff/min-length alternatively you can  switch off ResFinder\n"
			app.errorBox("Error", "Please select (at least) one database or insert a valid cutoff/min-length alternatively you can  switch off ResFinder")
			return
		 
		global virulencefinder_switch
		virulencefinder_switch=app.getRadioButton('virulencefinder_')
		virulencefinder_threshold=app.getEntry('virulencefinder_thre')
		virulencefinder_database=app.getOptionBox('DB: ')
		if(virulencefinder_switch == "yes" and virulencefinder_threshold == ""):
			print "Please select a valid cutoff or you can  switch off VirulenceFinder\n"
			app.errorBox("Error", "Please select a valid cutoff or you can  switch off VirulenceFinder")
			return
   
		prokka_switch=""
		prokka_switch=app.getRadioButton('prokka_')
		global prokka_full_path
		global parSNP_reference_path
		global parSNP_reference_fsa_path

		if prokka_full_path =="" and prokka_switch =="yes":
			print "Please select a valid prokka path or switch it off"
			app.errorBox("Error", "Please select a valid prokka path or switch it off")
			return
		if parSNP_switch =="yes":
		 	if parSNP_reference_fsa_path == "" and parSNP_reference_path =="":
				print "Please select a valid parSNP reference genbank or fasta file or switch it off"
				app.errorBox("Error", "Please select a valid parSNP reference genbank or fasta file or switch parSNP off")

				return

		output_switch="yes"#app.getRadioButton('output_')

		global resfams_switch
		resfams_switch=app.getRadioButton('resfams_')
		if resfams_switch== "yes" and prokka_switch == "no" and (TG_switch=="yes" or spades_switch == "yes"):
			print "Please switch \"on\" prokka and select its path or deactivate Resfams\n"
			app.errorBox("Error", "Please switch \"on\" prokka and select its path or deactivate Resfams")
			return
		if VirDBSearch_switch== "yes" and prokka_switch == "no" and (TG_switch=="yes" or spades_switch == "yes"):
			print "Please switch \"on\" prokka and select its path or deactivate VirDB\n"
			app.errorBox("Error", "Please switch \"on\" prokka and select its path or deactivate VirDB")
			return
		if cardSearch_switch== "yes" and prokka_switch == "no" and (TG_switch=="yes" or spades_switch == "yes"):
			print "Please switch \"on\" prokka and select its path or deactivate CARD\n"
			app.errorBox("Error", "Please switch \"on\" prokka and select its path or deactivate CARD")
			return
#changing the switch to the opposite to make sense
		steps_counter=0
		if(TG_switch=="no"):
			TG_switch="yes"
		else:
			TG_switch="no"
			steps_counter=steps_counter+1
		if(spades_switch=="no"):
			spades_switch="yes"
		else:
			spades_switch="no"
			steps_counter=steps_counter+1
		if(mlst_switch=="no"):
			mlst_switch="yes"
		else:
			mlst_switch="no"
			steps_counter=steps_counter+1
		if(plasmids_finder_switch=="no"):
			plasmids_finder_switch="yes"
		else:
			plasmids_finder_switch="no"
			steps_counter=steps_counter+1
		if(resfinder_switch=="no"):
			resfinder_switch="yes"
		else:
			resfinder_switch="no"
			steps_counter=steps_counter+1
		if(cardSearch_switch=="no"):
			cardSearch_switch="yes"
		else:
			cardSearch_switch="no"
			steps_counter=steps_counter+1
		if(emmTyping_switch == "no"):
			emmTyping_switch="yes"
		else:
			emmTyping_switch="no"
			steps_counter=steps_counter+1
		if(VirDBSearch_switch=="no"):
			VirDBSearch_switch="yes"
		else:
			VirDBSearch_switch="no"
			steps_counter=steps_counter+1
		if(parSNP_switch=="no"):
			parSNP_switch="yes"
		else:
			parSNP_switch="no"
			steps_counter=steps_counter+1
		if(virulencefinder_switch=="no"):
			virulencefinder_switch="yes"
		else:
			virulencefinder_switch="no"
			steps_counter=steps_counter+1
		if(prokka_switch=="no"):
			prokka_switch="yes"
		else:
			prokka_switch="no"
			steps_counter=steps_counter+1
		if(output_switch=="no"):
			output_switch="yes"
		else:
			output_switch="no"
			#steps_counter=steps_counter+1
		if(resfams_switch=="no"):
			resfams_switch="yes"
		else:
			resfams_switch="no"
			steps_counter=steps_counter+1
		#creating the yaml info
		os.chdir(pipeline_path)
		yaml_file="directories:\n  reads: \""+input_path+"\"\n  output: \""+output_path+"\"\n\ntrim_galore:\n  deactivate: "+TG_switch+"\n  reads_type: "+TG_paired+"\n  quality_threshold: "+TG_qual_thred+"\n\nspades:\n  deactivate: "+spades_switch+"\n  Mode: "+SPAdes_mode+"\n  kmer: \""+spades_kmer+"\"\n\nmlst_typing:\n  deactivate: "+mlst_switch+"\n  organism: "+mlst_organism+"\n\nplasmids_finder:\n  deactivate: "+plasmids_finder_switch+"\n  plasmids_database: "+plasmids_database_type+"\n  identity_threshold: "+plasmids_finder_threshold+"\n\nresfinder:\n  deactivate: "+resfinder_switch+"\n  identity_threshold: "+resfinder_threshold+"\n  resistance_database: "+resistance_database_input+"\n  min_length: "+resfinder_min_length+"\n\ncardSearch:\n  deactivate: "+cardSearch_switch+"\n\nemmTyping:\n  deactivate: "+emmTyping_switch+"\n\nvirulencefinder:\n  deactivate: "+virulencefinder_switch+"\n  identity_threshold: "+virulencefinder_threshold+"\n  virulence_database: "+virulencefinder_database+"\n\nparSNP:\n  deactivate: "+parSNP_switch+"\n  parSNP_reference: "+parSNP_reference_path+"\n  parSNP_reference_fsa: "+parSNP_reference_fsa_path+"\n\nVirDBSearch:\n  deactivate: "+VirDBSearch_switch+"\n\nprokka:\n  deactivate: "+prokka_switch+"\n  prokka_path: "+prokka_full_path+"\n\nOutput:\n  deactivate: "+output_switch+"\n\nResfams:\n  deactivate: "+resfams_switch+"\n"
		files_processed = open('configure.yaml', 'w')
		files_processed.write(yaml_file)
		files_processed.close()
		app.setTabbedFrameDisabledTab("TabbedFrame", "Results", disabled=True)
		log_short=output_path+'/pipeline.log'
		log=output_path+'log.txt'
		if os.path.isfile(log):
			files_processed = open(log, 'w')
			files_processed.write('')
			files_processed.close()
		if os.path.isfile(log_short):
			files_processed = open(log_short, 'w')
			files_processed.write('')
			files_processed.close()

############################
		Processors_=app.getEntry('Processors')
	 	global Quit_Flag
		Quit_Flag=0
		habl_=str(Processors_)

		try:
			thread.start_new_thread(pipeline,("configure.yaml",habl_,))
			thread.start_new_thread(check_progress,(log_short,log,steps_counter,))
			print ""
		except:
			print "Error: unable to start thread"

######################Main script GUI
######################Main script GUI
######################Main script GUI

app = gui("Bacterial Whole Genome Sequencing Analysis Pipeline",'780x800')#,'962x785')#,'569x625')#,'569x625'
#app.setSticky("w")
app.setExpand("both")
#app.setStretch("none")
app.setPadding([20,10]) # 20 pixels padding outside the widget [X, Y]
app.setInPadding([20,10])
app.setAllRadioButtonWidths(0)
#app.setAllRadioButtoAlign("left")

app.startTabbedFrame("TabbedFrame")
app.startTab("Settings")





text="Welcome to BacPipe \n Whole Genome Sequencing Pipeline \n"
#app.showSplash(text, fill='blue', stripe='black', fg='white', font=44)
raw = 0
#app.addLabel("title", "Welcome to BacPipe: Whole Genome Sequencing Pipeline for Clinical Diagnostics", raw, 0, 2)  # Row 0,Column 0,Span 2
raw=raw+1

app.addVerticalSeparator(3, 2, 0, 40 , colour="black")

#app.addButtons( ["login"], log_in,raw,1)
##app.addLabel("input", "Input:", raw, 2)
#app.setLabelAlign("input", "e")# Row 2,Column 0

# Row 1,Column 0
app.addButtons(["Choose input folder"], browse_input, raw, 3) # Row 3,Column 0,Span 2
#app.addEntry("input", 1, 1)                           # Row 1,Column 1
#raw=raw+1
#app.addLabel("output", "Output:", raw, 4)
#app.setLabelAlign("output", "e")# Row 2,Column 0
app.addButtons(["Choose output folder "], browse_output, raw, 4, 1) # Row 3,Column 0,Span 2
#app.addEntry("output", 2, 1)                     # Row 2,Column 1
#app.addEntry("TG_", 3, 1)                     # Row 2,Column 1





raw=raw+1
app.addHorizontalSeparator(raw,0,6, colour="black")
raw=raw+1
app.addLabel("Trim Galore", "Trim Galore", raw, 1,1,0)
app.setLabelAlign("Trim Galore", "e")
app.setLabelBg("Trim Galore", "white")

app.addNamedButton("?","TG_H", trim_galore_help, raw, 2)

app.addRadioButton("Trim Galore","yes",raw,4)
app.addRadioButton("Trim Galore","no",raw,3)
#app.addCheckBox("Install Cutadapt (first time run)",raw,5)
app.setRadioTick("Trim Galore", tick=True)
app.setRadioButtonAlign("Trim Galore", "e")
app.setRadioButtonAnchor("Trim Galore", "e")
#app.setRadioButtonFunction("Trim Galore", press_Trim_Galore)
raw=raw+1
#app.addLabelOptionBox("Reads", ["- Illumina -", "paired", "single", "- Pacbio -", "single"], raw,3)
app.addLabelOptionBox("Reads", ["paired", "single"], raw,3)
app.setOptionBox("Reads",0, value=True)

#raw=raw+1
app.addLabel("TG_qual", "Q.cutoff:", raw, 4,1)              # Row 2,Column 0
app.addEntry("TG_qual", raw, 5,0)                     # Row 2,Column 1
app.setEntry("TG_qual","25")
app.setLabelAlign("TG_qual", "e")
#app.setLabelAnchor("TG_qual", "w")
raw=raw+1
app.addHorizontalSeparator(raw,0,6, colour="black")
raw=raw+1
app.addLabel("spades", "SPAdes", raw, 1,1,0)
app.setLabelAlign("spades", "e")
app.setLabelBg("spades", "white")

app.addNamedButton("?", "SP_H",spades_help, raw, 2)

app.addRadioButton("spades","yes",raw,4)
app.addRadioButton("spades","no",raw,3)
app.setRadioTick("spades", tick=True)
raw=raw+1
#app.addLabelOptionBox("Mode", ["- Illumina -", "paired", "single", "- Pacbio -", "pacbio","- Iontorrent -", "iontorrent","- Nanopore -", "nanopore","- Sanger -", "sanger"], raw,3)
app.addLabelOptionBox("Mode", ["- Illumina -", "paired", "single", "- Pacbio -", "single","- Iontorrent -", "single","- Nanopore -", "single","- Sanger -", "single"], raw,3)
app.setOptionBox("Mode",1, value=True)

app.addLabel("spades_kmer_", "Kmer:", raw, 4)              # Row 2,Column 0
app.addEntry("spades_kmer_", raw, 5)                     # Row 2,Column 1
app.setEntry("spades_kmer_","77")#21,33,55,77
app.setLabelAlign("spades_kmer_", "w")

raw=raw+1
app.addHorizontalSeparator(raw,0,6, colour="black")
raw=raw+1
app.addLabel("mlst_", "MLST typing", raw, 1,1,0)
app.setLabelAlign("mlst_", "e")
app.setLabelBg("mlst_", "white")
app.addRadioButton("mlst_","yes",raw,4)
app.addRadioButton("mlst_","no",raw,3)
app.setRadioTick("mlst_", tick=True)
raw=raw+1
app.addLabelOptionBox("Species", ["ecoli","ecoli_2","kpneumoniae","paeruginosa","pfluorescens","abaumannii","abaumannii_2","senterica","saureus","spneumoniae","spyogenes","efaecalis","efaecium","----------------------------","achromobacter","aeromonas","afumigatus","aphagocytophilum","arcobacter","bcc","bcereus","bhampsonii","bhenselae","bhyodysenteriae","bintermedia","blicheniformis","bordetella","borrelia","bpilosicoli","bpseudomallei","brachyspira","bsubtilis","calbicans","campylobacter","cbotulinum","cconcisus","cdifficile","cdiphtheriae","cfetus","cfreundii","cglabrata","chelveticus","chlamydiales","chyointestinalis","cinsulaenigrae","ckrusei","clanienae","clari","cmaltaromaticum","cronobacter","csepticum","csinensis","csputorum","ctropicalis","cupsaliensis","ecloacae","fpsychrophilum","hcinaedi","hinfluenzae","hparasuis","hpylori","hsuis","kkingae","koxytoca","kseptempunctata","leptospira","leptospira_2","leptospira_3","llactis","lmonocytogenes","lsalivarius","mabscessus","magalactiae","mbovis","mcatarrhalis","mhaemolytica","mhyopneumoniae","mhyorhinis","mmassiliense","mplutonius","mpneumoniae","neisseria","orhinotracheale","otsutsugamushi","pacnes","pgingivalis","plarvae","pmultocida_multihost","pmultocida_rirdc","ppentosaceus","ranatipestifer","sagalactiae","sbsec","scanis","sdysgalactiae","sepidermidis","sgallolyticus","shaemolyticus","shominis","sinorhizobium","slugdunensis","smaltophilia","soralis","spseudintermedius","ssuis","sthermophilus","sthermophilus_2","streptomyces","suberis","szooepidemicus","taylorella","tenacibaculum","tvaginalis","vcholerae","vibrio","vparahaemolyticus","vtapetis","vvulnificus","wolbachia","xfastidiosa","yersinia","ypseudotuberculosis","yruckeri"], raw,3,2)
#app.addLabelOptionBox("Species", ["ecoli_2","ecoli","kpneumoniae","paeruginosa","pfluorescens","abaumannii","abaumannii_2","senterica","saureus","spneumoniae","spyogenes","efaecalis","efaecium","----------------------------","abaumannii","abaumannii_2","achromobacter","aeromonas","afumigatus","aphagocytophilum","arcobacter","bbacilliformis","bcc","bcereus","bhampsonii","bhenselae","bhyodysenteriae","bintermedia","blicheniformis","bordetella","borrelia","bpilosicoli","bpseudomallei","brachyspira","brucella","bsubtilis","calbicans","campylobacter","cbotulinum","cconcisus","cdifficile","cdiphtheriae","cfetus","cfreundii","cglabrata","chelveticus","chlamydiales","chyointestinalis","cinsulaenigrae","ckrusei","clanienae","clari","cmaltaromaticum","cronobacter","csepticum","csinensis","csputorum","ctropicalis","cupsaliensis","dnodosus","ecloacae","edwardsiella","fpsychrophilum","ganatis","hcinaedi","hinfluenzae","hparasuis","hpylori","hsuis","kaerogenes","kkingae","koxytoca","kseptempunctata","leptospira","leptospira_2","leptospira_3","liberibacter","llactis","lmonocytogenes","lsalivarius","mabscessus","magalactiae","mbovis","mcanis","mcaseolyticus","mcatarrhalis","mhaemolytica","mhyopneumoniae","mhyorhinis","miowae","mmassiliense","mplutonius","mpneumoniae","msynoviae","mycobacteria","neisseria","orhinotracheale","otsutsugamushi","pacnes","pdamselae","pgingivalis","plarvae","pmultocida_multihost","pmultocida_rirdc","ppentosaceus","pputida","psalmonis","ranatipestifer","rhodococcus","sagalactiae","sbsec","scanis","sdysgalactiae","sepidermidis","sgallolyticus","shaemolyticus","shominis","sinorhizobium","slugdunensis","smaltophilia","soralis","sparasitica","spseudintermedius","ssuis","sthermophilus","sthermophilus_2","streptomyces","suberis","szooepidemicus","taylorella","tenacibaculum","tpallidum","tvaginalis","ureaplasma","vcholerae","vcholerae2","vibrio","vparahaemolyticus","vtapetis","vvulnificus","wolbachia","xfastidiosa","yersinia","ypseudotuberculosis","yruckeri"], raw,3,2)
app.setOptionBox("Species",0, value=True)
app.setOptionBoxHeight("Species", 1)
app.setOptionBoxWidth("Species", 5)
raw=raw+1

app.addLabel("emm_", "Emm typing", raw, 1,1,0)
app.setLabelAlign("emm_", "e")
app.setLabelBg("emm_", "white")
app.addNamedButton("?", "em_H",emm_help, raw, 2)

app.addRadioButton("emm_","yes",raw,4)
app.addRadioButton("emm_","no",raw,3)
app.setRadioTick("emm_", tick=True)
raw=raw+1
app.addHorizontalSeparator(raw,0,6, colour="black")

raw=raw+1
app.addLabel("plasmids_finder_", "PlasmidFinder", raw, 1,1,0)
app.setLabelAlign("plasmids_finder_", "e")
app.setLabelBg("plasmids_finder_", "white")
app.addNamedButton("?", "PF_H",Plasmindfinder_help, raw, 2)

app.addRadioButton("plasmids_finder_","yes",raw,4)
app.addRadioButton("plasmids_finder_","no",raw,3)
app.setRadioTick("plasmids_finder_", tick=True)
raw=raw+1
#app.addLabelOptionBox("Database ", ["enterobacteriaceae","Inc18","NT_Rep","Rep1","Rep2","Rep3","RepA_N","RepL","Rep_trans"], raw,3,2)
#app.setOptionBox("Database ",0, value=True)
#app.addLabel("Database ", "plasmids_database:", 11, 0)              # Row 2,Column 0
raw=raw+1
app.addLabel("plasmids_finder_thre", "Cutoff:", raw, 3)              # Row 2,Column 0
app.addEntry("plasmids_finder_thre", raw, 4)                     # Row 2,Column 1
app.setEntry("plasmids_finder_thre", "95")
app.setLabelAlign("plasmids_finder_thre", "w")
raw=raw+1
app.addHorizontalSeparator(raw,0,6, colour="black")

raw=raw+1
app.addLabel("resfinder_", "ResFinder", raw, 1,1,0)
app.setLabelAlign("resfinder_", "e")
app.setLabelBg("resfinder_", "white")
app.addNamedButton("?", "RF_H",resfinder_help, raw, 2)

app.addRadioButton("resfinder_","yes",raw,4)
app.addRadioButton("resfinder_","no",raw,3)
app.setRadioTick("resfinder_", tick=True)
raw=raw+1
app.addLabel("resfinder_thre", "Cutoff:", raw, 3)              # Row 2,Column 0
app.addEntry("resfinder_thre", raw, 4)                     # Row 2,Column 1
app.setEntry("resfinder_thre","95")
app.setLabelAlign("resfinder_thre", "w")

#raw=raw+1
#app.addLabel("Antibiotics", "Resistance database:", raw, 5)
#app.setLabelAlign("Antibiotics", "w")
#app.addEntry("Antibiotics", 15, 1)                     # Row 2,Column 1


resistance_DB={"aminoglycoside":True, "beta-lactam":True, "colistin":True,
 "fosfomycin":True, "fusidicacid":True, "macrolide":True, "nitroimidazole":True, "oxazolidinone":True,
 "phenicol":True, "quinolone":True, "rifampicin":True, "sulphonamide":True, "tetracycline":True,
 "trimethoprim":True, "glycopeptide":True}


app.startToggleFrame("Antibiotics",raw,5)
app.addProperties("Antibiotics", resistance_DB)
##app.setPropertiesFunction("Antibiotics", changed)
app.stopToggleFrame()

'''
raw=raw+1
app.addLabel("resfinder_min_lngth", "ResFinder min length:", raw, 3)              # Row 2,Column 0
app.addEntry("resfinder_min_lngth", raw, 4)                     # Row 2,Column 1
app.setEntry("resfinder_min_lngth","0.60")
app.setLabelAlign("resfinder_min_lngth", "w")

'''
raw=raw+1


app.addLabel("card_", "CARD", raw, 1,1,0)
app.setLabelAlign("card_", "e")
app.setLabelBg("card_", "white")
app.addNamedButton("?", "CA_H",CARD_help, raw, 2)

app.addRadioButton("card_","yes",raw,4)
app.addRadioButton("card_","no",raw,3)
app.setRadioTick("card_", tick=True)
raw=raw+1
app.addLabel("resfams_", "ResFams", raw, 1)
app.setLabelAlign("resfams_", "e")
app.setLabelBg("resfams_", "white")
app.addNamedButton("?", "rF_H",resFams_help, raw, 2)

app.addRadioButton("resfams_","yes",raw,4)
app.addRadioButton("resfams_","no",raw,3)

raw=raw+1
app.addHorizontalSeparator(raw,0,6, colour="black")

raw=raw+1
app.addLabel("virulencefinder_", "VirulenceFinder", raw, 1,1,0)
app.setLabelAlign("virulencefinder_", "e")
app.setLabelBg("virulencefinder_", "white")
app.addNamedButton("?", "VF_H",virulenceFinder_help, raw, 2)

app.addRadioButton("virulencefinder_","yes",raw,4)
app.addRadioButton("virulencefinder_","no",raw,3)
app.setRadioTick("virulencefinder_", tick=True)
raw=raw+1
app.addLabel("virulencefinder_thre", "Cutoff:", raw, 3)              # Row 2,Column 0
app.addEntry("virulencefinder_thre", raw, 4)                     # Row 2,Column 1
app.setEntry("virulencefinder_thre","95")
app.setLabelAlign("virulencefinder_thre", "w")

#raw=raw+1
#app.addLabel("virulencefinder_db", "virulencefinder database:", raw, 3)              # Row 2,Column 0
#app.addLabelOptionBox("DB: ", ["virulence_ecoli","virulence_ent","listeria","s.aureus_exoenzyme","s.aureus_hostimm","s.aureus_toxin","stx"], raw,5,2)
app.addLabelOptionBox("DB: ", ["eaec","listeria","s.aureus","s.aureus_exoenzyme","s.aureus_hostimm","s.aureus_toxin","stx","virulence_ecoli","virulence_ent"], raw,5,2)
app.setOptionBox("DB: ",7, value=True)


raw=raw+1
app.addLabel("VirDB_", "VirDB", raw, 1,1,0)
app.setLabelAlign("VirDB_", "e")
app.setLabelBg("VirDB_", "white")
app.addNamedButton("?", "Vi_H",VirDB_help, raw, 2)
app.addRadioButton("VirDB_","yes",raw,4)
app.addRadioButton("VirDB_","no",raw,3)
app.setRadioTick("VirDB_", tick=True)
raw=raw+1
app.addHorizontalSeparator(raw,0,6, colour="black")

raw=raw+1

app.addLabel("parSNP_", "ParSNP", raw, 1,1,0)
app.setLabelAlign("parSNP_", "e")
app.setLabelBg("parSNP_", "white")
app.addNamedButton("?", "BS_H",BarSNP_help, raw, 2)
app.addRadioButton("parSNP_","yes",raw,4)
app.addRadioButton("parSNP_","no",raw,3)
app.setRadioTick("parSNP_", tick=True)
raw=raw+1
app.addLabel("parSNP_reference_", "Reference genome", raw, 3,0)              # Row 1,Column 0
app.addButtons(["Genbank file"], browse_parSNP, raw, 4, 1) # Row 3,Column 0,Span 2
app.setLabelAlign("parSNP_reference_", "w")
#raw=raw+1
#app.addLabel("parSNP_reference_fsa", "or", raw, 4,0)              # Row 1,Column 0
app.addButtons(["Fasta file"], browse_parSNP_fsa, raw, 5, 1) # Row 3,Column 0,Span 2
#app.setLabelAlign("parSNP_reference_fsa", "e")



raw=raw+1

app.addHorizontalSeparator(raw,0,6, colour="black")

raw=raw+1
app.addLabel("prokka_", "Prokka", raw, 1)
app.setLabelAlign("prokka_", "e")
app.setLabelBg("prokka_", "white")
app.addNamedButton("?", "Pr_H",Prokka_help, raw, 2)

app.addRadioButton("prokka_","yes",raw,4)
app.addRadioButton("prokka_","no",raw,3)


#app.addLabel("prokka_", "Output Path:", 20, 0)              # Row 2,Column 0
#app.addEntry("prokka_", 20, 1)                     # Row 2,Column 1
#raw=raw+1
#app.addLabel("prokka_path", "Prokka path:", raw, 3)              # Row 1,Column 0
app.addButtons(["Prokka path"], browse_prokka, raw, 5, 1) # Row 3,Column 0,Span 2
#app.setLabelAlign("prokka_path", "w")

raw=raw+1
app.addHorizontalSeparator(raw,0,6, colour="black")

#app.addLabel("prokka_path", "prokka_path:", 21, 0)              # Row 2,Column 0
#app.addEntry("prokka_path", 21, 1)                     # Row 2,Column 1

raw=raw+1

#app.addHorizontalSeparator(raw,0,5, colour="black")

#app.addLabel("resfams_", "Output Path:", 23, 0)              # Row 2,Column 0
#app.addEntry("resfams_", 23, 1)
'''
raw=raw+1
app.addLabel("output_", "Summary", raw, 1)
app.setLabelAlign("output_", "e")
app.setLabelBg("output_", "white")
app.addNamedButton("?", "Su_H",Summarize_help, raw, 2)

app.addRadioButton("output_","yes",raw,4)
app.addRadioButton("output_","no",raw,3)

#app.addLabel("output_", "Output Path:", 22, 0)              # Row 2,Column 0
#app.addEntry("output_", 22, 1)
raw=raw+1
#app.addHorizontalSeparator(raw,0,5, colour="black")
'''
raw=raw+1
app.addLabel("Processors", "Processors", raw, 1)
app.setLabelAlign("Processors", "e")
app.setLabelBg("Processors", "white")
app.addNamedButton("?", "Prc_H",Proc_help, raw, 2)

app.addEntry("Processors", raw, 3)                     # Row 2,Column 1
app.setEntry("Processors","4")

raw=raw+1
app.addHorizontalSeparator(raw,0,6, colour="black")

raw=raw+1
app.addButtons(["Start", "Cancel"], press, raw, 2, 2) # Row 3,Column 0,Span 2

raw=raw+1

#app.startLabelFrame("  ", 1, 0,1,5)
app.startLabelFrame("  ", 1, 0,1,50)

app.setLabelFrameRelief("  ","flat")#, "raised", "groove", "ridge", "flat"

#app.addImage("  ", "image_1.gif")
app.addImage("  ", "Presentation1.gif")
app.setLabelFrameRelief("  ","flat")#, "raised", "groove", "ridge", "flat"
app.stopLabelFrame()

'''
app.startLabelFrame("   ", 6, 0,1,6)
app.setLabelFrameRelief("   ","flat")#, "raised", "groove", "ridge", "flat"

app.addImage("   ", "logo.gif")
app.stopLabelFrame()

app.startLabelFrame("", 13, 0,1,30)
app.setLabelFrameRelief("","flat")#, "raised", "groove", "ridge", "flat"

app.addImage("logo1", "images_.gif")
app.stopLabelFrame()
'''


app.stopTab()
app.startTab("Progress")
#app.addLabel("l2", "Progress log will be shown after starting the analysis",0,0,0,0)

app.addScrolledTextArea("t1",2,0,20,20)
app.addMeter("progress_",1)


#app.addScrolledTextArea("t1_")


app.stopTab()

app.startTab("Results")
#app.addLabel("l3", "Tab 3 Label")
app.addLabel("l55", "Results")
app.addTextArea("tree")
app.setTextArea("tree", "")
	#app.startLabelFrame("     ", 13, 0,1,30)
##app.setLabelFrameRelief("     ","flat")
	#app.addImage("phy_tree_", "mytree.gif")#mytree.svg
	#app.shrinkImage("phy_tree_", 5)
##app.setImageSize("phy_tree_", 2, 5)
	#app.stopLabelFrame()

hh=[["Tool", "Results","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","",""]]
app.addGrid("Results grid",  hh,0,0,action=None, addRow=True)

#app.addButton("Resistance results", Resistance_detail,5,0) # Row 3,Column 0,Span 2


app.stopTab()
app.setTabbedFrameDisabledTab("TabbedFrame", "Results", disabled=True)

app.startTab("Help")
#app.addLabel("l4", "Tab 4 Label")
help_txt="""
 			 |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 			 ||                      Welcome To LMM WGS pipeline                ||
			 || Software for analysing whole genome sequencing data||
			 ||     for clinical diagnostics and outbreaks assessment     ||
			 ||                                 Copyright (C) 2017                            ||
			 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 
 This is free software, under GNU, also it includes other software also under GNU Copyright
 
 ####################################################################################################
 Settings and inputs:
 The pipeline can be run in whole or from intermediate/specific step(s), thus it is important to know the expected input for each step if you want to run specific step(s). It start with the raw reads (forward and reverse *.fastq.gz) and process them as following:
 a) Trim Galore: quality trimming
	INPUT: *_R1*.fastq.gz/*_R2*.fastq.gz
	OUTPUT: *_R1_001_val_1.fq.gz/*_R2_001_val_2.fq.gz
b) SPAdes: assembly and scaffolding
	INPUT: *_R1_001_val_1.fq.gz/*_R2_001_val_2.fq.gz
	OUTPUT: *.fasta
c) MLST typing:
	INPUT: *.fasta (DNA sequences)
d) emm typing (for Streptococcus):
	INPUT: *.fasta (DNA sequences)
e) Plasmids Finder:
	INPUT: *.fasta (DNA sequences)
f) ResFinder (Plasmid mediated Resistance only):
	INPUT: *.fasta (DNA sequences)
g) Virulence Finder:
	INPUT: *.fasta (DNA sequences)
h) Prokka (annotation):
	INPUT: *.fasta (DNA sequences)
i) CARD search (extensive resistance search):
	INPUT: *.faa (poteins sequences)
j) VirDB search:
	INPUT: *.faa (poteins sequences)
k) Outbreak assessment (parSNP)
	INPUT: *.fasta (DNA sequences, > 2 files)
l) ResFams search (extensive resistance search):
	INPUT: *.faa (poteins sequences)
m) Summarise output:
	OUTPUT: excel file for each sample (each selected tools shown in one sheet)
	
	The progress is shown in the Progress tab, where for each sample the percentage of completed steps are shown. When the run is finished the results tab will be accessible and the results will be shown.
	
####################################################################################################
Output:
There are three types of outputs, overall, summarised and details. All will be produced at the end of the analysis.
	A) Overall: A cross sample illustration of the various tools results. It is useful where more than one sample are analysed. It only covers MLST, resFinder, virulenceFinder and PlasmidFinder tools results. In cases where the identified gene is found in more than one location within a sample, their % identify to that gene are shown within the same cell (Comma separated). Additionally, a text format of sap-diversity tree will be illustrated (when selecting the outbreak option).
	
	B) Summarised: An excel excel file for each sample will be created. Within this sheet, the tools results will be shown (one sheet per tool). These excel files will be grouped in the â€œSummaryâ€ folder.
	
	C) Detailed:
	A folder will be created per sample within the output directory. Within this folder a sub-folders for each of the tools output will be created stating the detailed results for each tools.
	
	Additionally, a log-file will be created stating the step-wise update of the tools performed (pipeline.log) and a detailed log-file with their output (log.txt) for any possible errors.
 """
app.addTextArea("help_")
app.setTextArea("help_", help_txt)
app.stopTab()

app.startTab("Citation")
#app.addLabel("l5", "Citation")
cite_txt=str('VFDB: \n  Yang J, Chen LH, Sun LL, Yu J and Jin Q, 2008. VFDB 2008 release: an enhanced web-based resource for comparative pathogenomics. Nucleic Acids Res. \nEmm Typing: \nFacklam R, Beall B, Efstratiou A, Fischetti V, Johnson D, Kaplan E, et al. emm Typing and Validation of Provisional M Types for Group A Streptococci. Emerg Infect Dis. 1999;5(2):247-253. \nParSNP:\n Treangen TJ, Ondov BD, Koren S, Phillippy AM (2014) Rapid Core-Genome Alignment and Visualization for Thousands of Microbial Genomes.\nCARD:\n Jia et al. 2017. CARD 2017: expansion and model-centric curation of the Comprehensive Antibiotic Resistance Database. Nucleic Acids Research, 45, D566-573.\nProkka:\n Seemann T, Prokka: Rapid Prokaryotic Genome Annotation, Bioinformatics, 2014 Jul 15;30(14):2068-9.\nVirulenceFinder:\n Joensen KG, Scheutz F, Lund O, Hasman H, Kaas RS, Nielsen EM, Aarestrup FM, Real-time whole-genome sequencing for routine typing, surveillance, and outbreak detection of verotoxigenic Escherichia coli. J. Clin. Micobiol. 2014. 52(5): 1501-1510.\nResFinder:\n Zankari E, Hasman H, Cosentino S, Vestergaard M, Rasmussen S, Lund O, Aarestrup FM, Larsen MV, Identification of acquired antimicrobial resistance genes, J Antimicrob Chemother. 2012 Jul 10.\nPlasmidFinder:\n Carattoli A, Zankari E, Garcia-Fernandez A, Voldby Larsen M, Lund O, Villa L, Aarestrup FM, Hasman H.,PlasmidFinder and pMLST: in silico detection and typing of plasmids. Antimicrob. Agents Chemother. 2014. April 28th.\nMLST:\n Larsen MV, Cosentino S, Rasmussen S, Friis C, Hasman H, Marvig RL, Jelsbak L, Sicheritz-Pontén T, Ussery DW, Aarestrup FM and Lund O., Multilocus Sequence Typing of Total Genome Sequenced Bacteria. J. Clin. Micobiol. 2012. 50(4): 1355-1361.\nResfams:\n Gibson MK, Forsberg KJ, Dantas G. Improved annotation of antibiotic resistance functions reveals microbial resistomes cluster by ecology. The ISME Journal. 2014, doi:ISMEJ.2014.106')
app.addTextArea("cite")
app.setTextArea("cite", cite_txt)
app.stopTab()

app.stopTabbedFrame()





#show_results()



app.go()


###################
###################
###################
###################
###################
###################
###################
###################


