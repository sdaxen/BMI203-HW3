from Bio import SeqIO
"""
implementation of homework 3 for the class BMI203 @UCSF
smith waterman algorithm implementation repurposed from: http://forrestbao.blogspot.in/2007/09/smith-waterman-algorithm-in-process.html
licence as stated:
	"This software is a free software. Thus, it is licensed under GNU General Public License version 3 or later.
	Forrest Bao, Sept. 26, 2007 <http://fsbao.net> <forrest.bao aT gmail.com>"

"""
import sys, string
from numpy import *
import numpy as np
import random
import pickle


#read in pairs from input files
posPairs = []
postivesF = open("Pospairs.txt", "r")
for line in postivesF.readlines():
	pair = line.split(" ")
	posPairs.append((pair[0], pair[1].replace("\n", "")))

negPairs = []
negativesF = open("Negpairs.txt", "r")
for line in negativesF.readlines():
	pair = line.split(" ")
	negPairs.append((pair[0], pair[1].replace("\n", "")))


def run(file1, file2, scoringMatrix, penalty, extension_penalty,npMatrix=-1):
	"""
	my method for running the smith waterman algorithm
	file 1 and file 2 are sequence files, scoring matrix is a string indicating the txt file scoring matrix,
	penalty is the new gap penalty, and npMatrix is an alternative to the scoringMatrix and is of type np array
	"""
	#get first sequence
	for seq_record in SeqIO.parse(file1, "fasta"):
		seq1 = seq_record.seq

	#get second sequence
	for seq_record in SeqIO.parse(file2, "fasta"):
		seq2 = seq_record.seq

	if type(npMatrix) is not np.ndarray: #in case an np array is NOT passed in as scoring matrix
		f3 = open(scoringMatrix,"r")
		scoreMat=[]
		for line in f3.readlines():
			if "#" in line or "A" in line or "R" in line:
				continue
			intLines = line.split(" ")
			cleaned = []
			for element in intLines:
				if element == "\n" or element == " " or element == "":
					continue
				else:
					cleaned.append(int(element))
			scoreMat.append(cleaned)

		scoreMat = np.asarray(scoreMat)
	else:
		scoreMat = npMatrix

	m,n =  len(seq1),len(seq2) #length of two sequences

	#generate DP table and traceback path pointer matrix
	score=zeros((m+1,n+1))   #the DP table
	pointer=zeros((m+1,n+1))  #to store the traceback path

	def match_score(alpha,beta,scoreMat):
		"""
		return the score between matching character alpha, with character beta, using similar matrix scoreMat
		"""
		alphabet={}
		alphabet["A"] = 0
		alphabet["R"] = 1
		alphabet["N"] = 2
		alphabet["D"] = 3
		alphabet["C"] = 4
		alphabet["Q"] = 5
		alphabet["E"] = 6
		alphabet["G"] = 7
		alphabet["H"] = 8
		alphabet["I"] = 9
		alphabet["L"] = 10
		alphabet["K"] = 11
		alphabet["M"] = 12
		alphabet["F"] = 13
		alphabet["P"] = 14
		alphabet["S"] = 15
		alphabet["T"] = 16
		alphabet["W"] = 17
		alphabet["Y"] = 18
		alphabet["V"] = 19
		alphabet["B"] = 20
		alphabet["Z"] = 21
		alphabet["X"] = 22
		alphabet["x"] = 22
		alphabet["-"] = 23
		lut_x=alphabet[alpha]
		lut_y=alphabet[beta]
		return scoreMat[lut_x][lut_y]

	
	max_score=0  #initial maximum score in DP table
	
	#calculate DP table and mark pointers
	for i in range(1,m+1): #go left to right, bottom to top in dp table
		for j in range(1,n+1):
			
			#if the coordinate pointed to is a gap, then we will add an extension penalty
			#pointer = 3 indicates a match, 2 indicates a gap by tracing left, and 1 indicates a gap by tracing up 
			if pointer[i-1][j] != 3:
				score_up=score[i-1][j]+extension_penalty
			else:
				score_up=score[i-1][j]+penalty
			if pointer[i][j-1] != 3:
				score_down=score[i][j-1]+extension_penalty
			else:
				score_down=score[i][j-1]+penalty
			#if diagonal, find the match score
			score_diagonal=score[i-1][j-1]+match_score(seq1[i-1],seq2[j-1],scoreMat)
			#tale the max score from the 4 options we have to set the pointer for coordinate i, j
			#0 introduced as an option to allow for local alignment instead of restricting to full global
			score[i][j]=max(0,score_up,score_down,score_diagonal)
			if score[i][j]==0:
			  pointer[i][j]=0 #0 means end of the path
			if score[i][j]==score_up:
			  pointer[i][j]=1 #1 means trace up
			if score[i][j]==score_down:
			  pointer[i][j]=2 #2 means trace left
			if score[i][j]==score_diagonal:
			  pointer[i][j]=3 #3 means trace diagonal
			if score[i][j]>=max_score:
			  max_i=i
			  max_j=j
			  max_score=score[i][j] #keep track of the max score in our DP table
	#END of DP table

	align1,align2='','' #initial sequences

	i,j=max_i,max_j #indices of path starting point

	#traceback, follow pointers
	while pointer[i][j]!=0:
		if pointer[i][j]==3: #go diagonal 
			align1=align1+seq1[i-1]
			align2=align2+seq2[j-1]
			i=i-1
			j=j-1
		elif pointer[i][j]==2: #go left
			align1=align1+'-'
			align2=align2+seq2[j-1]
			j=j-1
		elif pointer[i][j]==1:#go right 
			align1=align1+seq1[i-1]
			align2=align2+'-'
			i=i-1
	#END of traceback

	align1=align1[::-1] #reverse sequence 1
	align2=align2[::-1] #reverse sequence 2

	def calculateScore(align1,align2,scoreMat):
		"""
		calculate the final score of aligning align1 with align2 and using the scoring function defined by scoreMat
		"""
		score = 0
		for i in range(0,len(align1)):
			score += match_score(align1[i],align2[i],scoreMat)
		return score 
	return calculateScore(align1, align2, scoreMat) #/ float(min(len(seq1), len(seq2)))

def getAverageScores(matrix, gap, extension_penalty):
	"""
	average the scores across all positive pairs and compare this to the average score across all negative pairs
	"""
	avgPositives = 0
	for pair in posPairs:
		avgPositives += run(pair[0], pair[1], matrix, gap, extension_penalty)
	avgPositives = avgPositives / float(len(posPairs))
	avgNegatives = 0
	for pair in negPairs:
		avgNegatives += run(pair[0], pair[1], matrix,gap, extension_penalty)
	avgNegatives = avgNegatives / float(len(negPairs))
	print("+: " + str(avgPositives), "-:" + str(avgNegatives))
	return avgPositives, avgNegatives
def getCombinations():
	"""
	find all possible combinations of having an extension penalty [-1,-5] and a gap penalty of [-1.-20]
	"""
	combinations = []
	for i in range(-1, -20,-1):
		for j in range(-1,-5,-1):
			combinations.append((i,j))
	return combinations 

def findBestCombination(matrix):
	"""
	find the best combination of gap and extension penalty by trying all of them
	"""
	combinations = getCombinations()
	combinationResults = []
	for c in combinations:
		gap_penalty, extension_penalty = c[0], c[1]
		

		results = []
		allPairs = posPairs + negPairs
		for i in range(0, len(allPairs)):
			results.append((i, run(allPairs[i][0], allPairs[i][1], matrix, gap_penalty, extension_penalty)))
		results.sort(key=lambda x: x[1])
		TPcount = 0 #True positives, find the score at which 70% of positive pairs score higher than this score
		threshold = 0
		for j in range(len(results) - 1, 0, -1):
			if results[j][0] < 50:
				TPcount += 1
			if TPcount == 35: #70% of 50 pairs
				threshold = results[j][1]
		FPcount = 0 #false positives
		for j in range(len(results) - 1, 0, -1):
			if results[j][0] >= 50 and results[j][1] >= threshold: #if a negative pair and greater than threshold
				FPcount += 1


		true_positives = 0
		for i in range(50,100): #only look at top 50 highest scoring 
			if results[i][0] < 50: #if in positive pairs
				true_positives += 1
		combinationResults.append((true_positives, c, FPcount))
		print((true_positives, c, FPcount))
	combinationResults.sort(key=lambda x: x[0])
	print(combinationResults)

def generateROCs(penalty, extension_penalty, npMatrix=-1):
	"""
	function to generate the ROC curves given the penalty parameters and the scoring matrix npMatrix
	"""
	allPairs = posPairs + negPairs
	if type(npMatrix) is not np.ndarray:
		matrices = ["BLOSUM50", "BLOSUM62", "MATIO", "PAM100", "PAM250"]
	else: 
		matrices = ["optimized matrix"] #just to iterate once
	for mat in matrices:
		results = []
		for i in range(0, len(allPairs)):
			results.append((i, run(allPairs[i][0], allPairs[i][1], mat, penalty, extension_penalty, npMatrix=npMatrix)))
			results.sort(key=lambda x: x[1])
		print(results)
		predictions = []
		for j in range(0, len(results)):
			if results[j][0] < 50:
				predictions.append(1)
			else:
				predictions.append(-1)
		actual = [-1] * 50 + [1] * 50 #bottom scoring 50 should in theory be from neg pairs
		print(predictions)
		plotROC(actual, predictions, mat)

def plotROC(actual, predictions, mat):
	"""
	helper function to plot ROC curve
	"""

	from sklearn.metrics import roc_curve, auc
	import matplotlib.pyplot as plt
	import random
	fig, ax = plt.subplots()
	false_positive_rate, true_positive_rate, thresholds = roc_curve(actual, predictions)
	roc_auc = auc(false_positive_rate, true_positive_rate)
	plt.plot(false_positive_rate, true_positive_rate, 'b',
	label='AUC = %0.2f'% roc_auc)
	plt.legend(loc='lower right')
	plt.plot([0,1],[0,1],'r--')
	plt.xlim([0,1.0])
	plt.ylim([0,1.0])
	plt.ylabel('True Positive Rate')
	plt.xlabel('False Positive Rate')
	fig.savefig(mat + "_ROC_curve")
	plt.show()

def getFalsePosiveRate(truePositive, gap_penalty, extension_penalty):
	"""
	generates false positive rates for all of the given scoring matrices,
	truePositive is a percentage 0->1.0 to specify the true positive rate
	"""
	matrices = ["BLOSUM50", "BLOSUM62", "MATIO", "PAM100", "PAM250"]
	for matrix in matrices:
		results = []
		allPairs = posPairs + negPairs
		for i in range(0, len(allPairs)):
			results.append((i, run(allPairs[i][0], allPairs[i][1], matrix, gap_penalty, extension_penalty)))
		results.sort(key=lambda x: x[1])
		print(results)
		TPcount = 0 #True positives, find the score at which 70% of positive pairs score higher than this score
		threshold = 0
		for j in range(len(results) - 1, 0, -1):
			if results[j][0] < 50:
				TPcount += 1
			if TPcount == truePositive * 50: #70% of 50 pairs
				threshold = results[j][1]
		FPcount = 0 #false positives
		for j in range(len(results) - 1, 0, -1):
			if results[j][0] >= 50 and results[j][1] >= threshold: #if a negative pair and greater than threshold
				FPcount += 1
		print( matrix, FPcount / float(50))

def optimizationScore(matrix, gap_penalty, extension_penalty, np_matrix=-1):
	"""
	scoring function to optimize for part 2
	"""
	results = []
	totalScore = 0
	allPairs = posPairs + negPairs
	for i in range(0, len(allPairs)):
		results.append((i, run(allPairs[i][0], allPairs[i][1], matrix, gap_penalty, extension_penalty, npMatrix=np_matrix)))

	results.sort(key=lambda x: x[1])

	##specified in spec
	# falsePositives = [0, .1, .2, .3]
	# for falsePos in falsePositives:
	# 	posCount = 0 
	# 	threshold = -100
	# 	for j in range(0, len(results)):
	# 		if results[j][0] < 50: #if pos pair 
	# 			posCount += 1
	# 		if posCount == falsePos * 50: #70% of 50 pairs
	# 			threshold = results[j][1]
	# 			break
	# 	#get true positive
	# 	TPcount = 0 #false positives
	# 	for j in range(0, len(results)):
	# 		if results[j][0] < 50 and results[j][1] >= threshold: #if a positive pair and greater than threshold
	# 			TPcount += 1
	# 	TPrate = TPcount / float(50)
	# 	totalScore += TPrate
	# return totalScore

	##my optimization function to match part 1
	TPcount = 0 #True positives, find the score at which 70% of positive pairs score higher than this score
	threshold = 0
	for j in range(len(results) - 1, 0, -1):
		if results[j][0] < 50:
			TPcount += 1
		if TPcount == 35: #70% of 50 pairs
			threshold = results[j][1]
	FPcount = 0 #false positives
	for j in range(len(results) - 1, 0, -1):
		if results[j][0] >= 50 and results[j][1] >= threshold: #if a negative pair and greater than threshold
			FPcount += 1
	return 50 - FPcount


def optimize(scoringMatrix):
	"""
	runs a genetic optimization algorithm to construct an optimized scoring matrix, plots the ROC curve for this matrix
	"""
	def mutate(matrix):
		"""
		method to introduce mutations into the scoring matrix, randomly selects coordinate and adds random offset to it
		while prserving symmetry of matrix 
		"""
		for c in range(0, 2):
			i, j = random.randint(0,23), random.randint(0,23)
			newVal = matrix[i][j] + random.randint(1, 3) * random.choice([-1, 1])
			matrix[i][j] = newVal
			matrix[j][i] = newVal
		return matrix

	def replaceDiagonal(matrix):
		"""
		replaces matrix diagonal with that of BLOSUM50
		"""
		f3 = open("BLOSUM50","r")
		B50 =[]
		for line in f3.readlines():
			if "#" in line or "A" in line or "R" in line:
				continue
			intLines = line.split(" ")
			cleaned = []
			for element in intLines:
				if element == "\n" or element == " " or element == "":
					continue
				else:
					cleaned.append(int(element))
			B50.append(cleaned)
		for i in range(0, 24):
			matrix[i][i] = B50[i][i]
		return matrix

	f3 = open(scoringMatrix,"r")
	scoreMat=[]
	for line in f3.readlines():
		if "#" in line or "A" in line or "R" in line:
			continue
		intLines = line.split(" ")
		cleaned = []
		for element in intLines:
			if element == "\n" or element == " " or element == "":
				continue
			else:
				cleaned.append(int(element))
		scoreMat.append(cleaned)
	scoreMat = np.asarray(scoreMat)
	bestScore = optimizationScore(scoreMat, -7,-3, np_matrix=scoreMat)
	print("init best score", bestScore)
	
	for i in range(0, 35):
		if i == 0:
			newMat = replaceDiagonal(scoreMat)	
		else:
			newMat = mutate(scoreMat)
		newScore = optimizationScore(newMat, -7,-3, np_matrix=newMat)
		print("new score", newScore)
		if newScore > bestScore:
			scoreMat = newMat
			bestScore = newScore
		print("best score", bestScore)
	
	print("final matrix", scoreMat)
	generateROCs(-7,-3,npMatrix=scoreMat)
	np.savetxt('optimizedMatioMatrix.txt', scoreMat)
	pickle.dump( scoreMat, open( "optimizedMatioScoreMat.p", "wb" ))




##calls to run methods
# run("sequences/prot-0004.fa", "sequences/prot-0018.fa", "myBLOSOM62", -4, -2)
# getAverageScores("BLOSUM62", -7, -3)
# getAverageScores("BLOSUM50", -7, -3)
# getAverageScores("optimizedMatrix.txt", -7, -3)

getAverageScores("MATIO", -7, -3)
getAverageScores("optimizedMatio", -7, -3)

# scoreMat = pickle.load(open("optimizedMatioScoreMat.p", "rb"))
# f = open("blah", "w")
# for ar in scoreMat:
# 	for element in ar:
# 		f.write(str(element) + " ")
# 	f.write('\n')




# findBestCombination("BLOSUM50")
# generateROCs(-7, -3)
# getFalsePosiveRate(.7, -7, -3)

# optimizationScore("BLOSUM50", -7, -3)
# optimizationScore("BLOSUM62", -7, -3)
# optimizationScore("PAM100", -7, -3)
# optimizationScore("MATIO", -7, -3)

# optimize("BLOSUM50")
# optimize("MATIO")
# optimize("PAM100")



