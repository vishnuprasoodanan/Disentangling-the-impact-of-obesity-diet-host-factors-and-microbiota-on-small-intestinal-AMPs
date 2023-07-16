#This python Script is written by Vishnu Prasoodanan P K
#It takes the input file with inter-sample distances (with 3 columns). First column is "Sample1" second column is "Sample2" and third column is the inter sample distance (between Sample1 and Sample2)
#Inter-sample distance can be Weighted or unweighted unifrac distances or Bray-Curtis distance

import sys

# Read the file path from command line argument
file_path = sys.argv[1]
file2_path = sys.argv[2]
output_file = sys.argv[3]
output_file2 = sys.argv[4]

# Initialize dictionaries
dict1 = {}
dict2 = {}

# Read the file line by line and split the content based on tab delimiter
with open(file_path, 'r') as file:
	for line in file:
		elements = line.strip().split('\t')
		if len(elements) >= 8:
			key = elements[0]
			value1 = elements[6]
			value2 = elements[7]
			dict1[key] = value1
			dict2[key] = value2

# Print the dictionaries as standard output
print("Dictionary 1:")
for key, value in dict1.items():
	print(f"{key}: {value}")

print("\nDictionary 2:")
for key, value in dict2.items():
	print(f"{key}: {value}")

file.close()

with open(file2_path, 'r') as file2:
	lines2 = [line.strip().split('\t') for line in file2]

with open(output_file, 'w') as outfile:
	values = ["Sample1", "Sample2", "Value", "Gender1", "Gender2"]
	# Join the values with tabs
	tab_delimited_string = '\t'.join(values)
	outfile.write(f"{tab_delimited_string}\n")
	for line2 in lines2:
		print(line2[0])
		search_key = line2[0]
		if search_key in dict1:
			print("yes")
			sam = '\t'.join(line2)
			outfile.write(f"{sam}\t{dict1[line2[0]]}\t{dict1[line2[1]]}\n")
		else:
			print("No match found.")
		

outfile.close()
with open(output_file2, 'w') as outfile2:
	values2 = ["Sample1", "Sample2", "Value", "Diet1", "Diet2"]
	# Join the values with tabs
	tab_delimited_string2 = '\t'.join(values2)
	outfile2.write(f"{tab_delimited_string2}\n")
	for line3 in lines2:
		print(line3[1])
		search_key = line3[1]
		if search_key in dict2:
			print("yes")
			sam2 = '\t'.join(line3)
			outfile2.write(f"{sam2}\t{dict2[line3[0]]}\t{dict2[line3[1]]}\n")
		else:
			print("No match found.")
