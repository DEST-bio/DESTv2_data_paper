import sys

def main():
	file = "moments_output_Transatlantic.tsv"
	num_cols = int(sys.argv[1]) #32 + 2#- 3 * 4 # Desired number of columns per line
	
	with open(file, 'r') as f:
		lines = f.readlines()
	
	with open(file, 'w') as f:
		for i, line in enumerate(lines):
			fields = line.split()

			# Check that the intended number of columns exists
			if len(fields) != num_cols:
				continue

			# Check output errors by confirming that the first column is a valid
			# model name
			if fields[0] not in ["admixture", "split", "two_epoch", "twosplits"]:
				continue
			
			f.write(line)

if __name__ == "__main__":
	main()
