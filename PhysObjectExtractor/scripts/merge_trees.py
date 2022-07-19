import pathlib
import sys

import uproot


def merge_trees(f_in_path: pathlib.Path, f_out_path: pathlib.Path) -> None:
    all_branches = {}

    with uproot.open(f_in_path) as f_in:
        # filter out cycle number, skip top-level directories
        trees = sorted(set([k.split(";")[0] for k in f_in.keys() if "/" in k]))

        for tree in trees:
            tree_content = f_in[tree].arrays(library="ak") #, entry_stop=1000)
            tree_content = dict([(f, tree_content[f]) for f in tree_content.fields])
            all_branches.update(tree_content)

    with uproot.recreate(f_out_path) as f_out:
        f_out["events"] = all_branches


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("provide path to file as argument")
        raise SystemExit

    f_in_path = pathlib.Path(sys.argv[-1])
    # add _merged to name for output
    f_out_path = f_in_path.with_name(f_in_path.stem + "_merged" + f_in_path.suffix)

    merge_trees(f_in_path, f_out_path)
