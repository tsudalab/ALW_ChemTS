#!/sw/bin/python3.4
import sys, math

def AtomicWeight(Element):

	AtomicWeight = {'H': 1.008,'He': 4.003, \
	'Li': 6.938, 'Be': 9.012, 'B': 10.81, 'C': 12.01, 'N': 14.01, 'O': 16.00, 'F': 19.00, 'Ne': 20.18, \
	'Na': 22.99, 'Mg': 24.80, 'Al':26.98, 'Si':28.08, 'P': 30.97, 'S': 32.06, 'Cl': 35.45, 'Ar': 39.95, \
	'K': 39.10, 'Ca': 40.08, 'Sc': 44.96, 'Ti': 47.87, 'V': 50.94, 'Cr': 52.00, 'Mn': 54.94,'Fe': 55.85, 'Co': 58.93,\
	'Ni': 58.69, 'Cu': 63.55, 'Zn': 65.38,'Ga': 69.72, 'Ge': 72.63, 'As': 74.92,'Se': 78.97, 'Br': 79.90, 'Kr': 83.80,\
	'Rb': 85.47, 'Sr': 87.62, 'Y': 88.91, 'Zr': 91.22, 'Nb': 92.21, 'Mo': 95.95, 'Tc': 99, 'Ru': 101.07, 'Rh': 102.91,\
	'Pd': 106.42, 'Ag': 107.87, 'Cd': 112.41, 'In': 114.82, 'Sn': 118.71, 'Sb': 121.76, 'Te': 127.60, 'I': 126.9, 'Xe': 131.29}

	if(Element in AtomicWeight):
		return AtomicWeight[Element]
	else:
		#print ("We don't have the information about %-s!" % (Element)) 
		exit()


def AtomicNumElec(Element):

	AtomicNumElec = {'H': 1,'He': 2, \
	'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10, \
	'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, \
	'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25,'Fe': 26, 'Co': 27,\
	'Ni': 28, 'Cu': 29, 'Zn': 30,'Ga': 31, 'Ge': 32, 'As': 33,'Se': 34, 'Br': 35, 'Kr': 36,\
	'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45,\
	'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54}

	if(Element in AtomicNumElec):
		return AtomicNumElec[Element]
	else:
		#print ("We don't have the information about %-s!" % (Element)) 
		exit()
	
def One_Atom_Energy(Element, Functional, Basis):

	#print (Functional)
	#print (Basis)

	if Functional == 'b3lyp':
#		#print (Functional)
		if Basis == '3-21g*':
#			#print (Basis)

			Atom_Energies = {'H': -0.497311436764,'He': -2.88600130979, \
			'Li': -7.43894355058, 'Be': -14.5838411624, 'B': -24.5186288965, \
			'C': -37.6426964829, 'N': -54.2954627446, 'O': -74.6602935932, \
			'F': -99.1821673337, 'Ne': -128.203678601}
	
		elif Basis == '6-31g*':
#			#print (Basis)

			Atom_Energies = {'H': -0.500272784186,'He': -2.90704897442, \
			'Li': -7.49098472494, 'Be': -14.6684425431, 'B': -24.6543539563, \
			'C': -37.8462799747, 'N': -54.5844893898, 'O': -75.0606214291, \
			'F': -99.7155354580, 'Ne': -128.894359950 }

	#print (Atom_Energies)
	
	if(Element in Atom_Energies):
		return Atom_Energies[Element]
	else:
		#print ("We don't have the information about %-s!" % (Element)) 
		exit()
 
 

 
#var = input()
#print (AtomicNumElec(var))
#print (One_Atom_Energy(var,'B3LYP', '6-31G*'))

