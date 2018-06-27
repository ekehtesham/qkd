#!/usr/bin/python

"""

@title		: QuTiP simulation of KMB09 for paper 'A Survey of QKD Protocols Based upon Network Implementation'
@author		: ekhan (Ehtesham Khan)
@contact	: ehtesham_khan@live.com
@github		: @ekehtesham
@created on	: Sat 20 Mar 2018
@updated on	: Thu 14 Jun 2018

"""

from qutip import *
from collections import OrderedDict
import qutip as qt
import numpy as np
import statistics as st
import matplotlib.pyplot as plt
import pandas as pd
import multiprocessing
import pickle
import os
import datetime

"""
Standard Basis |0> 	   	 & |1>       for the rectilinear basis.
Hadamard Basis |0> + |1> & |0> - |1> for the diagonal basis.
0 0 -> |0>
0 1 -> |1>
1 0 -> |0> + |1>
1 1 -> |0> - |1>
"""


#This class is used for the qubit methods: initialization, measurements of both standard and hadamard basis & representation of qubit
#class start "Qubit"
class Qubit():
	
	#initialization
	def __init__(self,initial_state):
		if initial_state:
			self.__state = qt.basis(2,1)
		else:
			self.__state = qt.basis(2,0)
		self.__isMeasured = False
		self.__hadamard = snot()
		self.__standard = basis(2,0).dag()

	#standard measurement
	def StandardMeasurement(self):
		if self.__isMeasured:
			raise Exception("Qubit already measured!")
		M = 1000000
		m = np.random.randint(0,M)
		self.__isMeasured = True
		if m < round(pow(np.dot(self.__standard,self.__state).norm(),2),2)*M:
			return 0
		else:
			return 1

	#hadamard measurement
	def HadamardMeasurement(self):
		if self.__isMeasured:
			raise Exception("Qubit already measured!")
		self.__state = self.__hadamard*self.__state

	#representation of qubit
	def Show(self):
		var = ""
		
		if round(np.dot(basis(2,0).dag(),self.__state).norm(),2):
			var += "{0}|0>".format(str(round(np.dot(basis(2,0).dag(),self.__state).norm(),2)) if round(np.dot(basis(2,0).dag(),self.__state).norm(),2) != 1.0 else '')
		if round(np.dot(basis(2,1).dag(),self.__state).norm(),2):
			if var:
				var += " + "
			var += "{0}|1>".format(str(round(np.dot(basis(2,1).dag(),self.__state).norm(),2)) if round(np.dot(basis(2,1).dag(),self.__state).norm(),2) != 1.0 else '')

		return var.ljust(17)[:17]
#class end "Qubit"

#This class is used for the user communication methods: initialization, sending qubits & receiving qubits
#class start "User"
class User():

	#initialization
	def __init__(self,name):
		self.name = name

	#sending qubits
	def SendQubits(self,listOfBits,listOfBasis,bitsIndeces):
		
		assert len(listOfBits) == len(listOfBasis), "Basis and Bits must be the same length!"

		qubits = list()
		
		for i in range(len(listOfBits)):
			#for basis = 1, lets consider its standard basis
			if not listOfBasis[i]:
				#Standard Basis
				if not listOfBits[i]:
					qubits.append(Qubit(0))
					bitsIndeces.append('1')
				else:
					qubits.append(Qubit(0))
					bitsIndeces.append('2')

			#for basis = 0, lets consider its hadamard basis
			else:
				#Hadamard Basis
				if not listOfBits[i]:
					var = Qubit(1)
					bitsIndeces.append('1')
				else:
					var = Qubit(1)
					bitsIndeces.append('2')
					
				var.HadamardMeasurement()
				qubits.append(var)

		return qubits

	#receving qubits
	def ReceiveQubits(self,listOfQubits,listOfBasis,bitsIndeces):

		assert len(listOfQubits) == len(listOfBasis), "Basis and Bits must be the same length!"
		
		bits = list()
		
		for i in range(len(listOfQubits)):
			#for basis = 1, lets consider its standard basis
			if not listOfBasis[i]:
				#Standard Measurement
				bits.append(listOfQubits[i].StandardMeasurement())
				if not bits[i]:
					bitsIndeces.append('1')
				else:
					bitsIndeces.append('2')
			#for basis = 0, lets consider its hadamard basis
			else:
				listOfQubits[i].HadamardMeasurement()
				bits.append(listOfQubits[i].StandardMeasurement())
				#Hadamard Measurement
				if not bits[i]:
					bitsIndeces.append('1')
				else:
					bitsIndeces.append('2')

		return bits

	#normalize qubits
	def NormalizeQubits(self,listOfBits,listOfBasis):

		assert len(listOfBits) == len(listOfBasis), "Basis and Bits must be the same length!"
		
		qubits = list()
		
		for i in range(len(listOfBits)):
			#for basis = 1, lets consider its standard basis
			if not listOfBasis[i]:
				#Standard Basis
				if not listOfBits[i]:
					qubits.append(Qubit(0))
				else:
					qubits.append(Qubit(0))

			#for basis = 0, lets consider its hadamard basis
			else:
				#Hadamard Basis
				if not listOfBits[i]:
					var = Qubit(1)
				else:
					var = Qubit(1)
				var.HadamardMeasurement()
				qubits.append(var)

		return qubits
#class end "User"

#static method for the generation of array of random bits i.e. 0 & 1
def GenerateRandomBits(no_of_qubits):
	var = list()
	for i in range(no_of_qubits):
		rnd = np.random.randint(0,2)
		var.append(rnd)
	return var	

#static method for converting bits and basis into classical symbols
def ConvertToSymbols(listOfBits, listOfBasis, isBasis):
	var = list()
	if isBasis:
		for i in range(NO_OF_QUBITS):
			if not listOfBits[i]:
				var.append("+")
			else:
				var.append("X")
	else:
		for i in range(NO_OF_QUBITS):

			#for basis = 1, lets consider its standard basis
			if not listOfBasis[i]:
				#Standard Basis
				if not listOfBits[i]:
					var.append("|")
				else:
					var.append("━")

			#for basis = 0, lets consider its hadamard basis
			else:
				#Hadamard Basis
				if not listOfBits[i]:
					var.append("╱")
				else:
					var.append("╲")
	return var	

def ClearScreen():
	os.system('cls' if os.name == 'nt' else 'clear')

#===================================================================

#This class is used for the qkd protocol execution methods: initialization & execute for no. of qubits
#class start "QKDProtocol"
class QKDProtocol(multiprocessing.Process):

	def __init__(self):
		multiprocessing.Process.__init__(self)

	def exportDataToFile(self,fileName,data):
		with open(fileName+'.txt','ab') as textFile:
			np.savetxt(textFile, ["%s\n" % data], fmt='%s')

	def Execute(self):
		alice = User("Alice")
		aliceIndeces = list()
		aliceBasis = GenerateRandomBits(NO_OF_QUBITS)
		aliceBits = GenerateRandomBits(NO_OF_QUBITS)
		aliceSymBasis = ConvertToSymbols(aliceBasis, aliceBasis, True)
		aliceSymBits = ConvertToSymbols(aliceBits, aliceBasis, False)
		qubits = alice.SendQubits(aliceBits, aliceBasis, aliceIndeces)
		aliceQubits = qubits

		if EVE_EXIST:
			eve = User("Eve")
			eveIndeces = list()
			eveBasis = GenerateRandomBits(NO_OF_QUBITS)
			eveBits = eve.ReceiveQubits(qubits, eveBasis, eveIndeces)
			eveSymBasis = ConvertToSymbols(eveBasis, eveBasis, True)
			eveSymBits = ConvertToSymbols(eveBits, eveBasis, False)
			qubits = eve.SendQubits(eveBits, eveBasis, eveIndeces)
			eveQubits = qubits

		bob = User("Bob")
		bobIndeces = list()
		bobBasis = GenerateRandomBits(NO_OF_QUBITS)
		bobBits = bob.ReceiveQubits(qubits, bobBasis, bobIndeces)
		bobSymBasis = ConvertToSymbols(bobBasis, bobBasis, True)
		bobSymBits = ConvertToSymbols(bobBits, bobBasis, False)
		
		aliceKey = list()
		bobKey = list()
		mismatchIndeces = list()
		correctQubits = list()

		for i in range(NO_OF_QUBITS):
			if aliceIndeces[i] != bobIndeces[i]:
				mismatchIndeces.append("1")

				if aliceIndeces[i] == '1' or aliceIndeces[i] == '1':
					if bobSymBits[i] == "━":
						aliceKey.append(aliceBits[i])
						bobKey.append(1)
						correctQubits.append('1')
					elif bobSymBits[i] == "╲":
						aliceKey.append(aliceBits[i])
						bobKey.append(0)
						correctQubits.append('0')
					else:
						correctQubits.append('x')

				else:
					if bobSymBits[i] == "|":
						aliceKey.append(aliceBits[i])
						bobKey.append(1)
						correctQubits.append('1')
					elif bobSymBits[i] == "╱":
						aliceKey.append(aliceBits[i])
						bobKey.append(0)
						correctQubits.append('0')
					else:
						correctQubits.append('x')

			else:
				mismatchIndeces.append("0")
				correctQubits.append('x')


		if aliceKey != bobKey:
			key = False
			length = None

			if not SILENT:
				print("Basis was correct but key mismatch, eve is present.")
				print("Alice Key 	:", aliceKey)
				print("Bob Key  	:", bobKey)
			
			errorQubit = 0
			
			for i in range(len(aliceKey)):
				if aliceKey[i] != bobKey[i]:
					errorQubit += 1
			
			qber = np.round((errorQubit / len(aliceKey)),5)*100
			QBERs.append(qber)
	
			if not SILENT:
				print("QBER		:", qber, "%")

		else:
			if len(bobKey) == 0:
				key = False
				length = None

				if not SILENT:
					print("No qubit is successfully transffered.")
			else:
				key = True
				length = len(bobKey)
				QBERs.append(0)

				if not SILENT:
					print("QBER		: 0 %")
				
				if not SILENT:
					print("Successfully exchanged key!")
					print("Key Length 	: " + str(length))
					print("Key 		:", aliceKey)
			
		if LOG and not SILENT:
			print(" ")
			print(" ")
			print("--------BIT REPRESENTATION--------")
			print("Alice Basis        :", aliceBasis)
			print("Qubits |ei>        :", aliceBits)
			if EVE_EXIST:
				print("Eve Basis |gk>     :", eveBasis)
				print("Qubits |<gk|ei>|^2 :", eveBits)
				print("Bob Basis |ej>     :", bobBasis)
				print("Qubits |<ej|gk>|^2 :", bobBits)
			else:
				print("Bob Basis |ej>     :", bobBasis)
				print("Qubits |<ej|ei>|^2 :", bobBits)
			print(" ")
			print("Alice Indeces      :", aliceIndeces)
			print("Bob Indeces        :", bobIndeces)
			print(" ")
			print("Mismatch Index     :", mismatchIndeces)
			print(" ")
			print("Selected Bits      :", correctQubits)
			print(" ")
			print(" ")

			print("--------SYMBOL REPRESENTATION--------")
			print("Alice Basis        :", aliceSymBasis)
			print("Qubits |ei>        :", aliceSymBits)
			if EVE_EXIST:
				print("Eve Basis |gk>     :", eveSymBasis)
				print("Qubits |<gk|ei>|^2 :", eveSymBits)
				print("Bob Basis |ej>     :", bobSymBasis)
				print("Qubits |<ej|gk>|^2 :", bobSymBits)
			else:
				print("Bob Basis |ej>     :", bobSymBasis)
				print("Qubits |<ej|ei>|^2 :", bobSymBits)
			print(" ")
			print("Mismatch Index     :", mismatchIndeces)
			print(" ")
			print(" ")
			
			print("-------QUBIT REPRESENTATION-------")
			var = "Alice Qubits |ei>      : "
			for qubit in aliceQubits:
				var += qubit.Show() + "   "
			print(var)

			if EVE_EXIST:
				var = "Bob Qubits |<ej|ei>|^2 : "
				for qubit in eveQubits:
					var += qubit.Show() + "   "
				print(var)

			print(" ")
			print(" ")
			
			print("Alice generates", format(str(NO_OF_QUBITS)), "random Basis")
			print("Alice sends to Bob", format(str(NO_OF_QUBITS)), "encoded Qubits")
			
			if EVE_EXIST:
				print("Eve generates", format(str(NO_OF_QUBITS)), "random Basis")
				print("Eve intercepts and decode Alice's", format(str(NO_OF_QUBITS)), "encoded Qubits")
				print("Eve sends to Bob as Alice's", format(str(NO_OF_QUBITS)), "encoded Qubits")
			
			print("Bob generates", format(str(NO_OF_QUBITS)), "random Basis")
			print("Bob receives and decode Alice's", format(str(NO_OF_QUBITS)), "encoded Qubits")

			print(" ")
			print(" ")


#===================================================================

QBERs = list()

#number of cycles in a QKD simulation
NO_OF_CYCLES = 1000

#number of qubits in a QKD simulation
NO_OF_QUBITS = 10

#scenario of eve
EVE_EXIST = True

#logging
LOG = False

#graph only
SILENT = True

ClearScreen()

#execution of protocol
qkd = QKDProtocol()

if SILENT:
	print("Executing", NO_OF_CYCLES, "Cycle(s) of", NO_OF_QUBITS, "Qubit(s)... Please wait!", datetime.datetime.now())

for i in range(NO_OF_CYCLES):
	if not SILENT:
		print("==================================")
		print("            CYCLE(S)", i + 1, "           ")
		print("==================================")
	#else:
		#ClearScreen()
	#	print("Executing", NO_OF_CYCLES, "Cycle(s) of", NO_OF_QUBITS, "Qubit(s)... Please wait!", np.round((i+1)/NO_OF_CYCLES*100,1), "% Completed")
	qkd.Execute()

if SILENT:
	print(len(QBERs), "Cycle(s) successfully executed! Generating Plot...", datetime.datetime.now())

avg = np.round(st.mean(QBERs),2)
print("Avg. QBER =", avg, "≈", int(np.round(avg,0)), datetime.datetime.now())

#exporting results into text file
qkd.exportDataToFile('kmb09_qkd_results', avg)

#plotting results
df = pd.DataFrame({'cycle': range(len(QBERs)), 'qber': QBERs})

for x in range(len(QBERs)):
	scatterPlot = plt.scatter(x, QBERs[x], c=(0,0,1), s=7, label='QBER')

if NO_OF_CYCLES > 1:
	x1 = np.array(df['cycle'])
	y1 = np.array(df['qber'])
	z = np.polyfit(x1, y1, 1)
	p = np.poly1d(z)
	linePlot = plt.plot(x1, p(x1), "r--", label='Avg. QBER')

plt.xlabel("CYCLE(S)")
plt.ylabel("QBER (%)")
plt.title("QBER of KMB09 for " + str(len(QBERs)) + " Cycle(s) of " + str(NO_OF_QUBITS) + " Qubit(s) Per Cycle")

handles, labels = plt.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))

plt.legend(by_label.values(), by_label.keys())
plt.yticks(np.arange(0, 100 + 1, 5))
plt.grid(axis='y', linestyle='-')
plt.show()

#===================================================================
