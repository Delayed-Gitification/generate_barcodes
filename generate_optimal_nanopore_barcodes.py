import random
import argparse
from rapidfuzz import fuzz

nts = ["A", "T", "C", "G"]

def rev_c(seq):
	"""
	simple function that reverse complements a given sequence
	"""
	tab = str.maketrans("ACTGN", "TGACN")
	# first reverse the sequence
	seq = seq[::-1]
	# and then complement
	seq = seq.translate(tab)
	return seq


def fast_fuzz(s1, s2):
	if (s1 in s2) or (s2 in s1):
		return 100
	else:
		return fuzz.partial_ratio(s1, s2)


def find_barcode(seq, barcodes, min_score, max_ambiguity, stored_results):
	scores = []
	names = []

	# # print(seq)
	# # print(barcodes)

	if seq in stored_results.keys():
		return stored_results[seq], stored_results

	for name, bc in barcodes.items():
		score = fast_fuzz(seq, bc)
		scores.append(score)
		names.append(name)

		# print(' '.join([name, bc, str(score)]))

	# # print(scores)

	# Reject if similar to >1 barcode (ambiguous)
	# # print(max_ambiguity)
	if len([a for a in scores if a > max_ambiguity]) > 1:
		# # print([a for a in scores if a > max_ambiguity])
		stored_results[seq] = -1
	# # print("ambiguous")

	elif max(scores) >= min_score:
		stored_results[seq] = names[scores.index(max(scores))]
	# # print("yeh!")

	else:
		stored_results[seq] = -1
	# # print("nah")

	### print(stored_results[seq])

	return stored_results[seq], stored_results


def mutate_nt(nt, rate):

	if random.random() <= rate:
		return random.choice(nts)
	else:
		return nt


def mutate_seq(seq, rate):
	return [mutate_nt(i, rate) for i in seq]


def insertion(seq, rates):
	rand = random.random()
	cumulative_rate = 0
	for l, r in enumerate(rates):
		cumulative_rate += r
		if rand <= cumulative_rate:
			# chose a random insertion
			ins = random.choices(nts, k=l+1)

			# chose a random position
			pos = random.choice(list(range(len(seq))))

			new_seq = seq[0:pos] + ins + seq[pos:-1]
			return new_seq
	return seq  # if didn't return new seq


def deletion(seq, rates):
	rand = random.random()
	cumulative_rate = 0
	for l, r in enumerate(rates):
		cumulative_rate += r
		if rand <= cumulative_rate:

			# chose a random position
			pos = random.choice(list(range(l+2, len(seq)-1)))

			new_seq = seq[0:pos-l-1] + seq[pos:-1]
			return new_seq
	return seq  # if didn't return new seq


def mutate_and_indel(seq, times, error_rate, indel_rate):
	seq = list(seq)
	for _ in range(times):
		seq = mutate_seq(seq, error_rate)
		seq = insertion(seq, [0.25*indel_rate, 0.125*indel_rate, 0.06125*indel_rate])
		seq = deletion(seq, [0.25 * indel_rate, 0.125 * indel_rate, 0.06125 * indel_rate])

	return ''.join(seq)




def are_these_bcs_good(forward_bcs, reverse_bcs, padding, times, error_rate, indel_rate,
                       length, min_score, max_ambiguity):
	# Step 1 - chose forward and reverse barcode
	forward_bc = random.choice(list(forward_bcs.values()))
	reverse_bc = random.choice(list(reverse_bcs.values()))

	# Step 2 - add random sequence to all of these
	seq_f = ''.join(random.choices(nts, k=padding)) + forward_bc + ''.join(random.choices(nts, k=padding))
	seq_r = ''.join(random.choices(nts, k=padding)) + reverse_bc + ''.join(random.choices(nts, k=padding))

	# Step 3 - add mutations and indels
	seq_f = mutate_and_indel(seq_f, times, error_rate, indel_rate)
	seq_r = mutate_and_indel(seq_r, times, error_rate, indel_rate)

	full_seq = seq_f + ''.join(random.choices(nts, k=50)) + rev_c(seq_r)

	# print(full_seq)

	# Step 5 - determine whether these sequences can be demultiplexed
	can_demulti = analyse_demulti(full_seq, length, forward_bcs, reverse_bcs,
	                              min_score, max_ambiguity)


	return can_demulti




def analyse_demulti(full_seq, length, forward_bcs, reverse_bcs, min_score, max_ambiguity):
	can_demulti = False  # assume false

	reverse_bcs = {key: rev_c(value) for key, value in reverse_bcs.items()}

	forward_results = []
	reverse_results = []
	stored_results1 = {}
	stored_results2 = {}

	for rc in [True, False]:
		if rc:
			seq = rev_c(str(full_seq))
		else:
			seq = str(full_seq)

		s1 = seq[0:length]
		s2 = seq[-length:]

		# print("")
		# print(rc)
		# print("barcode1")
		bc1, stored_results1 = find_barcode(s1, forward_bcs, min_score, max_ambiguity, stored_results1)
		# print("barcode2")
		bc2, stored_results2 = find_barcode(s2, reverse_bcs, min_score, max_ambiguity, stored_results2)
		### print(bc2)
		# bc2 = -1

		if rc:
			reverse_results = [bc1, bc2]
		else:
			forward_results = [bc1, bc2]

	if forward_results == [-1, -1] and reverse_results != [-1, -1]:
		# reverse results are good
		bc1 = reverse_results[0]
		bc2 = reverse_results[1]
	elif forward_results != [-1, -1] and reverse_results == [-1, -1]:
		bc1 = forward_results[0]
		bc2 = forward_results[1]

	if bc1 != -1 and bc2 != -1:
		can_demulti = True

	return can_demulti

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("--initial_barcodes")
	parser.add_argument("-e", "--error_rate", type=float, default=0.05)
	parser.add_argument("-d", "--indel_rate", type=float, default=0.03)
	parser.add_argument("-p", "--padding", type=int, default=15)
	parser.add_argument("-a", "--attempts", type=int, default=100)
	parser.add_argument("--barcode_length", type=int, default=40)
	parser.add_argument("-n", "--barcode_n", type=int, default=6)
	parser.add_argument("--length", type=int, default=80, help="how long to consider")
	parser.add_argument("--times", type=int, help="How many attempts to do mutations and indels", default=2)
	parser.add_argument("-s", "--min_score", type=float, default=85, help="When matching barcodes")
	parser.add_argument("--max_ambiguity", type=float, default=75, help="If another barcode has this score or higher, then ignore this read because it's ambiguous")
	parser.add_argument("--iterations", type=int, default=1000, help="Number of times to check each "
	                                                                 "barcode set")
	args = parser.parse_args()

	# print(args)

	forward_bcs = {"F1": "CTCTAATCGAGAATTCCGAT",
"F2": "TAAGCATTGAACATGTGTCC",
"F3": "GGTGTGTTGTGTGTTGTATT",
"F4": "AACCTCGTATTATAACACGA",
"F5": "TAGAACGCATAGCCATAATA",
"F6": "ACTTATATAGTGGCGTAACG"}


	reverse_bcs = {"R1": "AACAGTTAATGCGTTCATCA",
"R2": "TATCCTAGACTGCTAGTTCT",
"R3": "TCTTCTTATTAGTAGCCAAC",
"R4": "TCCATTCATAGAATCGATGG",
"R5": "ATAGTATTGGCCTCAGAAGT",
"R6": "GATATATTGAGAAGGCTTGG"}

	forward_bcs = {'1': 'GAGCGTCTGGAATTGGAGGTAAGATCGAGTCTCTAATTCC', 2: 'ATTATATATCAACATTGGCTCCGGCGACCTGTATTCACAA',
	 3: 'CGGTTAACTGAGAGAAGACCATGATGAATACAGATAAGCT', 4: 'CTGTACGTTCTTGTGGCATCAGTTGGTCCTTCCAACGTGA',
	 5: 'ATGCCGGAATGAACCTAATCTACTTATGCCGCTGTGGTGT', 6: 'CGGATAACCTTCGTTCCAACTTCTTCTGAACACATATCCT'}
	reverse_bcs = {'1': 'ATATGTCGCGTGCAGAGGTATGTGGACGTGATAATCAAGG', 2: 'GATAACTGATCCTTCCACTACAGAACTCCTTCGCTAATTG',
	 3: 'AGCTGCTTAGGTCACCGATTCGAGCGAGCGTCTCAATTAA', 4: 'TCTCTGCAACCTACTGGTATACATATTATATACTTCAACC',
	 5: 'CGGTGACCGTTGTTATATTGTTATCCAGGTAAGGATCTCC', 6: 'CAATCCACTTGACATGGTTGGTGCTCACGAACGTGAGGAT'}

	successful = 0
	for _ in range(args.iterations):
		if are_these_bcs_good(forward_bcs, reverse_bcs, args.padding, args.times,
		                      args.error_rate, args.indel_rate, args.length,
		                      args.min_score, args.max_ambiguity):
			successful += 1
	success_rate = successful / args.iterations

	print(success_rate)

	# Iteratively add barcodes and choose best ones
	for _ in range(5):
		forward_bcs = {str(1): make_barcode(args.barcode_length)}
		reverse_bcs = {str(1): make_barcode(args.barcode_length)}
		for n in range(2, args.barcode_n+1):
			results = {}
			for _ in range(args.attempts):
				# Try a barcode
				forward_bcs[n] = make_barcode(args.barcode_length)
				reverse_bcs[n] = make_barcode(args.barcode_length)
				successful = 0
				for _ in range(args.iterations):
					if are_these_bcs_good(forward_bcs, reverse_bcs, args.padding, args.times,
					                      args.error_rate, args.indel_rate, args.length,
					                      args.min_score, args.max_ambiguity):
						successful += 1
				success_rate = successful / args.iterations

				results[forward_bcs[n] + "_" + reverse_bcs[n]] = success_rate

				#print("Success rate: " + str(success_rate*100) + "%")
			best_key = max(results, key=results.get)
			print(results[best_key])
			forward_bcs[n] = best_key.split("_")[0]
			reverse_bcs[n] = best_key.split("_")[1]

		print(forward_bcs)
		print(reverse_bcs)



def make_barcode(length, GC_min=0.3, GC_max=0.5, max_run=2):
	while True:
		bad = False  # assume good
		possible_bc = ''.join(random.choices(nts, k=length))

		GC = len([i for i in possible_bc if i =="C" or i =="G"])/length

		if GC < GC_min:
			bad = True
		if GC > GC_max:
			bad = True
		if "A"*(max_run+1) in possible_bc:
			bad = True
		if "C"*(max_run+1) in possible_bc:
			bad = True
		if "T"*(max_run+1) in possible_bc:
			bad = True
		if "G"*(max_run+1) in possible_bc:
			bad = True

		if not bad:
			break
	return possible_bc


if __name__ == "__main__":
	main()







