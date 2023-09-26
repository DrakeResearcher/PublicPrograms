import sys

def main(args):
    if len(args) < 3:
        print("INCORRECT USE")
        print("  To use, pass two arguments to this program.")
        print("  args: two files two compare.")
        exit()

    compare(args[1], args[2], args[3:])

def compare(file1, file2, ignore):
    print("ERIC'S SIMPLE FILE COMPARATOR")

    diffs = []

    fd1 = open(file1, "r")
    fd2 = open(file2, "r")

    lines1 = fd1.readlines()
    lines2 = fd2.readlines()

    for i in range(min(len(lines1), len(lines2))):
        if lines1[i] != lines2[i] and str(i) not in ignore: # disregard line(s) because that line contains the current date.
            diffs.append(i + 1)

    if len(diffs) != 0:
        print("  >Test failed. diffs on lines: ", diffs)
    else:
        print("  >Test passed!")

if __name__ == "__main__":
    main(sys.argv)