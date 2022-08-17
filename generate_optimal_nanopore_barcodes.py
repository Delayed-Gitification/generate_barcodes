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
	parser.add_argument("-e", "--error_rate", type=float, default=0.01)
	parser.add_argument("-d", "--indel_rate", type=float, default=0.01)
	parser.add_argument("-p", "--padding", type=int, default=15)
	parser.add_argument("--length", type=int, default=50, help="how long to consider")
	parser.add_argument("--times", type=int, help="How many attempts to do mutations and indels", default=2)
	parser.add_argument("-s", "--min_score", type=float, default=92, help="When matching barcodes")
	parser.add_argument("--max_ambiguity", type=float, default=82, help="If another barcode has this score or higher, then ignore this read because it's ambiguous")
	parser.add_argument("--iterations", type=int, default=100, help="Number of times to check each "
	                                                                 "barcode set")
	args = parser.parse_args()

	# print(args)

	forward_bcs = {'f1': "TGCAGCACGACATCAGCACT", 'f2':"AGCCGATCGATCGATCAAGC"}
	reverse_bcs = {'r1':"AATTCGATCGGGACTAACAC", 'r2':"TCGACTACGCACTAGCATTA"}

	successful = 0
	for _ in range(args.iterations):
		if are_these_bcs_good(forward_bcs, reverse_bcs, args.padding, args.times,
		                      args.error_rate, args.indel_rate, args.length,
		                      args.min_score, args.max_ambiguity):
			successful += 1
	success_rate = successful / args.iterations

	print("Success rate: " + str(success_rate*100) + "%")


if __name__ == "__main__":
	main()






