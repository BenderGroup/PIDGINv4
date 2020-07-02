#Authors : Maria-Anna Trapotsi, Layla Hosseini-Gerami and Lewis Mervin 
#Emails : mat64@cam.ac.uk and lh605@cam.ac.uk
#Supervisor : Dr. A. Bender
#All rights reserved 2020
#Protein Target Prediction Tool trained on SARs from PubChem (Mined 21/06/16) and ChEMBL21
#Molecular Descriptors : 2048bit Morgan Binary Fingerprints (Rdkit) - ECFP4
#Dependencies : rdkit, sklearn, numpy

#libraries
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import cPickle
import zipfile
import bz2
import glob
import os
import sys
import math
import numpy as np
from multiprocessing import Pool
import multiprocessing
multiprocessing.freeze_support()
from optparse import OptionParser
import json

def introMessage():
        print '=============================================================================================='
        print ' Authors: Maria-Anna Trapotsi, Layla Hosseini-Gerami and Lewis Mervin\n Email:  mat64@cam.ac.uk, lh605@cam.ac.uk\n Supervisor: Dr. A. Bender'
        print ' Address: Centre For Molecular Informatics, Dept. Chemistry, Lensfield Road, Cambridge CB2 1EW'
        print '==============================================================================================\n'
        return

#optionparser options
parser = OptionParser()
parser.add_option('-f', dest='inf', help='Input smiles or sdf file (required)', metavar='FILE')
parser.add_option('-b', '--bioactivity', default=None, type=str, dest='bioactivity', help='Bioactivity threshold (can use multiple split by ",". E.g. "100,10"')
parser.add_option('-n', '--ncores', default=1, type=int, dest='ncores', help='No. cores (default: 1)')
parser.add_option('--organism', dest='organism', default=None, type=str, help='Organism filter (multiple can be specified using commas ",")')
parser.add_option('--orthologues', action='store_true', default=False, dest='ortho', help='Set to use orthologue bioactivity data in model generation')

(options, args) = parser.parse_args()
#check smiles of sdf (if not quit)
def check_Input():
        global options
        extension = options.inf.split('.')[-1]
        if extension not in ['smi','smiles','sdf']:
                print ' Warning [-f]: File type not "smi", "smiles" or "sdf". Interpreting input as SMILES'
        return extension

#Read Dictionary to filter  any models with 0 actives
if options.ortho == False:
    with open('dict_uniprot_2_actives.json', 'r') as fp:
        dict_uniprot_2_actives = json.load(fp)
if options.ortho == True:
    with open('dict_uniprot_2_actives_ortho.json', 'r') as fp:
	dict_uniprot_2_actives = json.load(fp)

#calculate 2048bit morgan fingerprints, radius 2
def calcFingerprints(smiles):
        m1 = Chem.MolFromSmiles(smiles)
        fp = AllChem.GetMorganFingerprintAsBitVect(m1,2, nBits=2048)
        return fp

#calculate fingerprints for chunked array of smiles
def arrayFP(inp):
        outfp = []
        outsmi = []
        for i in inp:
                try:
                        outfp.append(calcFingerprints(i))
                        outsmi.append(i)
                except:
                        print 'SMILES Parse Error: ' + i
        return outfp,outsmi

#import user query
def importQuery(in_file):
        query = open(in_file).read().splitlines()
        #collect IDs, if present
        if len(query[0].split()) > 1:
                ids = [line.split()[1] for line in query]
                query = [line.split()[0] for line in query]
        else:
                ids = None
        smiles_per_core = int(math.ceil(len(query) / N_cores)+1)
        chunked_smiles = [query[x:x+smiles_per_core] for x in xrange(0, len(query), smiles_per_core)]
        pool = Pool(processes=N_cores)  # set up resources
        jobs = pool.imap(arrayFP, chunked_smiles)
        processed_fp = []
        processed_smi = []
        for i, result in enumerate(jobs):
                processed_fp += result[0]
                processed_smi += result[1]
        pool.close()
        pool.join()
        #if IDs aren't present, use SMILES as IDs
        if not ids:
                ids = processed_smi
        return processed_fp, processed_smi, ids

#get info for uniprots
#Change the function from v2 to v3 based on the location of information in columns of the file and also the 'sep' used in the file
def getUniprotInfo():
        if os.name == 'nt': sep = '\\'
        else: sep = '/'
        if options.ortho == True:
            model_info = [l.split('\t') for l in open(os.path.dirname(os.path.abspath(__file__)) + sep + 'classes_in_model_ortho.txt').read().splitlines()]
        if options.ortho == False:
            model_info = [l.split('\t') for l in open(os.path.dirname(os.path.abspath(__file__)) + sep + 'classes_in_model_no_ortho.txt').read().splitlines()]
        return_dict = {l[0] : l[0:7] for l in model_info}
        return return_dict

#sim worker
def doSimSearch(model_name):
        if os.name == 'nt': sep = '\\'
        else: sep = '/'
        mod = model_name.split(sep)[-1].split('.')[0]
        try:
            if options.ortho == False:
                with zipfile.ZipFile(os.path.dirname(os.path.abspath(__file__)) + sep + 'no_ortho/bioactivity_dataset' + sep + mod + '.smi.zip', 'r') as zfile:
                    #print(zfile)       
                    comps = [i.split('\t') for i in zfile.open(mod + '.smi', 'r').read().splitlines()]
            if options.ortho == True:
		with zipfile.ZipFile(os.path.dirname(os.path.abspath(__file__)) + sep + 'ortho/bioactivity_dataset' + sep + mod + '.smi.zip', 'r') as zfile:
		    comps = [i.split('\t') for i in zfile.open(mod + '.smi', 'r').read().splitlines()]
        except IOError: return
        comps2 = []
        afp = []
        #select comps[1:], because in PIDGINv3 .smi files the first row is the header (change comps to comps[1:])
        #comp[9], comp[10],comp[11] and comp[12] corresponds to the indication of activity in 100, 10, 1 and 0.1 bioactivity threshold respectively. 
        #comp[2] is the Compound ID and specify the comparison of query molecules only with actives from ChEMBL
        if options.bioactivity == '0.1':
            ind = 12
        if options.bioactivity == '1':
            ind = 11
        if options.bioactivity == '100':
            ind = 9
        if options.bioactivity == '10':
            ind = 10

        for comp in comps[1:]:
            #print(comps[0])
            if comp[3]!="Smiles" and comp[ind]=='1' and comp[2].startswith('CHEMBL'):
                try:
                        afp.append(calcFingerprints(comp[3]))
                        comps2.append(comp)
                except: pass
            

        ret = []
        for i,fp in enumerate(querymatrix):
                sims = DataStructs.BulkTanimotoSimilarity(fp,afp)
                if len(sims) > 0:
                    idx = sims.index(max(sims))
                    ret.append([sims[idx], mod] + comps2[idx] + [smiles[i]])
                else:
                    ret.append(["NaN","NaN"] + ["NaN"] + ["NaN"])
        return ret

#prediction runner
def performSimSearch(models):
        sims_results = []
        pool = Pool(processes=N_cores, initializer=initPool, initargs=(querymatrix,smiles))  # set up resources
        jobs = pool.imap_unordered(doSimSearch, models)
        out_file2.write('Uniprot\tPref_Name\tGene ID\tTarget_Class\tOrganism\tPDB_ID\tDisGeNET_Diseases_0.06\t' + '\t'.join(map(str,smiles)) + '\n')
        for i, result in enumerate(jobs):
                percent = (float(i)/float(len(models)))*100 + 1
                sys.stdout.write(' Calculating Sims for Query Molecules: %3d%%\r' % percent)
                sys.stdout.flush()
                if result is not None:
                        sims_results += result
                        try:
                            out_file2.write('\t'.join(map(str,model_info[result[0][1]])))
                            for sim in result:
                                    out_file2.write('\t' + str(round(sim[0],3)))
                        except:
                            continue
                        out_file2.write('\n')
        pool.close()
        pool.join()
        return sims_results

#initializer for the pool
def initPool(querymatrix_,smiles_):
        global querymatrix, smiles
        querymatrix = querymatrix_
        smiles = smiles_

#main
if __name__ == '__main__':
        if os.name == 'nt': sep = '\\'
        else: sep = '/'
        #input_name = sys.argv[1]
        input_name=options.inf
        #N_cores = int(sys.argv[2])
        N_cores = int(options.ncores)
        desired_organism = options.organism
        introMessage()
        print ' Calculating Near-Neighbors for ' + input_name
        print ' Using ' + str(N_cores) + ' Cores'
        if options.ortho == False:
	    models = [modelfile for modelfile in glob.glob(os.path.dirname(os.path.abspath(__file__)) + sep + 'no_ortho/bioactivity_dataset' + sep + '*.zip')]
        if options.ortho == True:
	    models = [modelfile for modelfile in glob.glob(os.path.dirname(os.path.abspath(__file__)) + sep + 'ortho/bioactivity_dataset' + sep + '*.zip')]
	if options.organism == None:
		models_filtered=[]
        	for model_id in models:
                	model_uniprot = model_id.split('/')[-1].split('.')[0]
			if options.bioactivity=='100':
				index_val=0
			if options.bioactivity=='10':
                        	index_val=1
			if options.bioactivity=='1':
				index_val=2
			if options.bioactivity=='0.1':
				index_val=3
                	n_actives=dict_uniprot_2_actives[model_uniprot][index_val]
                	if n_actives>=1:
                    	    models_filtered.append(model_id)
		models=models_filtered
		model_info = getUniprotInfo()
                output_name = input_name + '_out_similarity_details.txt'
                output_name2 = input_name + '_out_similarity_matrix.txt'

        if options.organism is not None:
		model_info = getUniprotInfo()
                models = [mod for mod in models if model_info[mod.split(sep)[-1].split('.')[0]][4] == options.organism]
                models_filtered=[]
                for model_id in models:
                    model_uniprot = model_id.split('/')[-1].split('.')[0]
                    if options.bioactivity=='100':
                        index_val=0
                    if options.bioactivity=='10':
                        index_val=1
                    if options.bioactivity=='1':
                        index_val=2
                    if options.bioactivity=='0.1':
                        index_val=3
                    n_actives=dict_uniprot_2_actives[model_uniprot][index_val]
                    if n_actives>=1:
			models_filtered.append(model_id)
		models=models_filtered
		print(n_actives)
                output_name = input_name + '_out_similarity_details' + '_'+options.organism.replace(' ', '_') + '.txt'
                output_name2 = input_name + '_out_similarity_matrix' + '_'+options.organism.replace(' ', '_') + '.txt'
                print ' Predicting for organism : ' + desired_organism
        print ' Total Number of Classes : ' + str(len(models_filtered))
        out_file = open(output_name, 'w')
        out_file2 = open(output_name2, 'w')
        querymatrix,smiles,ids = importQuery(input_name)
        print ' Total Number of Query Molecules : ' + str(len(querymatrix))
        sims_results = performSimSearch(models)
        out_file.write('Uniprot\tPref_Name\tGene ID\tTarget_Class\tOrganism\tNear_Neighbor_ChEMBLID\tNear_Neighbor_Smiles\tNear_Neighbor_Bioactive_organism\tNear_Neighbor_conf_score\tNN_activity\tNN_Units\tInput_Compound\tSimilarity\n')
        for row in sorted(sims_results,reverse=True):
            try:
                out_file.write('\t'.join(map(str,model_info[row[1]][0:5]))+ '\t' + '\t'.join(map(str,row[4:6]))+'\t'+str(row[3])+'\t'+str(row[7])+'\t'+str(row[6])+'\t'+str(row[8])+'\t'+str(row[15])+'\t'+str(row[0])+ '\n')
            except:
                continue
        print '\n Wrote Results to: ' + output_name
        print ' Wrote Results to: ' + output_name2
        out_file.close()


