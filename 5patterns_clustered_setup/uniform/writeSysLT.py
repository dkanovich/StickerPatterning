# arguments - n (segment count), NA (count of A-polymer), NB, Box Length (L), fileName
import sys 

_, n, NA, NB, L, fName = sys.argv 
n, NA, NB, L = int(n), int(NA), int(NB), float(L)

L_bound = L + 20

with open(fName, 'w') as tf:
	tf.write(f'import "polyA_n{n}.lt" \n')
	tf.write(f'import "polyB_n{n}.lt" \n\n')

	tf.write(f'p1 = new polyA_n{n}[{NA}]\n')
	tf.write(f'p2 = new polyB_n{n}[{NB}]\n\n')

	tf.write('write_once("Data Boundary")\n')
	tf.write('{\n')
	tf.write(f'{-L_bound} {L_bound} xlo xhi\n')
	tf.write(f'{-L_bound} {L_bound} ylo yhi\n')
	tf.write(f'{-L_bound} {L_bound} zlo zhi\n')
	tf.write('}')
	



