import sys


if __name__ == "__main__":
    input_file = open(sys.argv[1], "r")
    output_file = open("out.txt", "w")
    lines = input_file.readlines()
    start_left = 0
    start_right = 0
    for line in lines:
        line = line.split(" ")
        if start_left == 0 or start_right == 0:
            start_left = int(line[2])
            start_right = int(line[6])
            continue
        out_line = line[1] + "\t" + str(int(line[2])-start_left) + "\t" + str(int(line[6]) - start_right)
        output_file.write(out_line + "\n")
        start_left = int(line[2])
        start_right = int(line[6])
    input_file.close()
    output_file.close()
