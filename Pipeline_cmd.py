# coding=utf-8
# coding=utf-8
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

global config_file

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


import argparse
import argparse as ap


import argparse

description="Welcome To BacPipe pipeline: \n a Software for analysing whole genome sequencing data  for clinical diagnostics and outbreaks assessment"
parser = argparse.ArgumentParser(prog='BacPipe',description=description)
parser.add_argument('--os',dest='os',choices=['unix','mac'],help='please select the operating system', required=True)
parser.add_argument('--config',dest='file',help='path to config file',nargs='*', required=True)
parser.add_argument('--processors',type=int,dest='processors',help='number of processors',default=4)
parser.add_argument('--version', action='version', version='%(prog)s 1.2.6')
args=parser.parse_args()
try:
	config_file=args.file[0]
	exists = os.path.isfile(config_file)
	print "using %s as config file" % (config_file)
		
except FileNotFoundError:
	print "%s file does not exist" % (config_file)

OSystem=args.os
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


global input_path
input_path=""
global output_path
output_path=""
global prokka_full_path
prokka_full_path=""
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




def pipeline(config_f, thread_):
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
			#print "python3.5 %s -i %s -o %s -p  %s -s %s -x" % (mlst_path,mlst_assembly,mlst_directory,mlst_DB_path,organism)
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
			#print "python3.5 %s -i  %s -o  %s -mp blastn -x -p  %s" % (virulencefinder_path,vir_assembly,vir_directory,virulencefinder_DB_path)
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
		logging.info("deactivating PlasmidsFinder  as specified")
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

def show_results(config_f):
	with open(config_f, 'r') as f:
		config = yaml.load(f)
	#deactivated tools
	if 'deactivate' in config['mlst_typing']:
		mlst_typing_flag= str(config['mlst_typing']['deactivate'])
	else:
		mlst_typing_flag= 'yes'
	if 'deactivate' in config['plasmids_finder']:
		plasmids_finder_flag= str(config['plasmids_finder']['deactivate'])
	else:
		plasmids_finder_flag= 'yes'
	if 'deactivate' in config['cardSearch']:
		cardSearch_flag= str(config['cardSearch']['deactivate'])
	else:
		cardSearch_flag= 'yes'
	if 'deactivate' in config['VirDBSearch']:
		VirDBSearch_flag= str(config['VirDBSearch']['deactivate'])
	else:
		VirDBSearch_flag= 'yes'
	if 'deactivate' in config['emmTyping']:
		emmTyping_flag= str(config['emmTyping']['deactivate'])
	else:
		emmTyping_flag= 'yes'
	if 'deactivate' in config['resfinder']:
		resfinder_flag= str(config['resfinder']['deactivate'])
	else:
		resfinder_flag= 'yes'
	if 'deactivate' in config['Resfams']:
		Resfams_flag= str(config['Resfams']['deactivate'])
	else:
		Resfams_flag= 'yes'
	if 'deactivate' in config['virulencefinder']:
		virulencefinder_flag= str(config['virulencefinder']['deactivate'])
	else:
		virulencefinder_flag= 'yes'


	output_path = config["directories"]["output"]
	
	if(mlst_typing_flag=="False"):
		try:
			global mlst_samples
			mlst_output_ = output_path+"/mlst_output.txt"
			try:
				os.remove(mlst_output_)
			except:
				print ""
			Output_organizer_path=os.path.join(pipeline_path+'/Output_organizer.pl')
			call(["perl",Output_organizer_path,mlst_samples,mlst_output_,"none"])
		except:
			print "Can't find MLST output"


	if(plasmids_finder_flag=="False"):

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
		except:
			print "Can't find Plasmidfinder output"

	if(virulencefinder_flag=="False"):
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
		except:
			print "Can't find virulencefinder output"

	if(resfinder_flag=="False"):
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
		except:
			print "Can't find resfinder output"


	if(emmTyping_flag=="False"):
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
		except:
			print "Can't find emmTyping output"

	if(Resfams_flag=="False"):
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
		except:
			print "Can't find resfams output"

	if(cardSearch_flag=="False"):
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
		except:
			print "Can't find cardSearch output"

	if(VirDBSearch_flag=="False"):
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
		except:
			print "Can't find VirDBSearch output"
#####################################
#####################################
#####################################

try:
	Processors_=str(args.processors)
	print "using %s thereads" % (Processors_)
		
except ValueError:
	print "using the default 4 threads"



	
pipeline(config_file,Processors_)
os.chdir(pipeline_path)
show_results(config_file)
