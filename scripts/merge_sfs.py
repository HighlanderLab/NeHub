import numpy as np
import glob
import sys


def read_fsc_sfs(filepath):
    with open(filepath) as f:
        lines = f.readlines()
    header_line1 = lines[0].strip()      # "1 observations"
    header_line2 = lines[1].strip()      # column labels
    data_line = lines[2].strip().split()
    data = np.array([float(x) for x in data_line])
    return header_line1, header_line2, data


def merge_fsc_sfs(sfs_files, output_path):
    total = None
    header1, header2 = None, None

    for f in sfs_files:
        h1, h2, data = read_fsc_sfs(f)
        if total is None:
            header1, header2 = h1, h2
            total = data
        else:
            if len(data) != len(total):
                raise ValueError(
                    f"Dimension mismatch: {f} has {len(data)} bins, "
                    f"but previous files have {len(total)} bins"
                )
            total += data

    with open(output_path, "w") as out:
        out.write(header1 + "\n")
        out.write(header2 + "\n")
        out.write("\t".join(map(str, total.astype(int))) + "\n")

    print(f"Merged {len(sfs_files)} files")
    print(f"Total sites after merging: {total.sum():.0f}")


if __name__ == "__main__":
    # Usage: python merge_fsc_sfs.py file1.obs file2.obs ... -o merged.obs
    args = sys.argv[1:]
    if "-o" not in args:
        print("Usage: python merge_fsc_sfs.py file1.obs file2.obs ... -o output.obs")
        sys.exit(1)
    out_idx = args.index("-o")
    output_path = args[out_idx + 1]
    sfs_files = args[:out_idx]
    merge_fsc_sfs(sfs_files, output_path)