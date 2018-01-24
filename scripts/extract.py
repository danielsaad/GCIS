import sys

def extract_from_file(f,l,r):
    f.seek(l,0)
    chunk = f.read(r-l+1)
    print(l,r,'=',chunk)


def extract(file_path, extract_file_path):
    with open(file_path) as f, open(extract_file_path) as f_ext:
        for line in f_ext:
            interval = [int(x) for x in line.split()] # read the numbers
            extract_from_file(f,interval[0],interval[1])

# argv[1] = file path
# argv[2] = extract input
if __name__ == "__main__":
    extract(sys.argv[1],sys.argv[2])
