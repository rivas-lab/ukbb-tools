import glob
import os

def do_symlinks(src_dir):
    dst_dir = "/oak/stanford/groups/mrivas/private_data/ukbb/24983/phenotypedata"
    phes = glob.glob(src_dir + "phe/*phe")
    phe_names = [os.path.splitext(os.path.basename(phe))[0] for phe in phes]
    print("Symlinking phe/info/log files in source directory " + src_dir + "...")
    for phe_name, phe in zip(phe_names, phes):
        log = os.path.join(
            os.path.dirname(os.path.dirname(phe)), "logs/{0}.log".format(phe_name)
        )
        info_dir = os.path.join(os.path.dirname(os.path.dirname(phe)), "info")
        info_file = os.path.join(info_dir, "{}.info".format(phe_name))
        # Remove existing symlinks if they exist
        if os.path.exists(
            os.path.join(dst_dir, "current", "phe/{}.phe".format(phe_name))
        ):
            print("Removing:")
            print(os.path.join(dst_dir, "current", "phe/{}.phe".format(phe_name)))
            for folder, filetype in zip(
                ["phe", "info", "logs"], ["phe", "info", "log"]
            ):
                path = os.path.join(
                        dst_dir, "current", folder, phe_name + "." + filetype
                    )
                if os.path.exists(path):
                    os.remove(path)
        # Add symlinks if the paths exist
        for path, folder, filetype in zip(
            [phe, info_file, log], ["phe", "info", "logs"], ["phe", "info", "log"]
        ):
            if os.path.exists(path):
                print(path, os.path.join(
                    dst_dir, "current", folder, phe_name + "." + filetype
                    )
                )
                os.symlink(
                    path,
                    os.path.join(
                    dst_dir, "current", folder, phe_name + "." + filetype
                    ),
                )
            else:
                print("FUCK")

# Clinical breathing indices, High-confidence QC, Cancer, Family History
cbi_dir = "/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/clinical_breathing_indices/"
hc_dir = "/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/highconfidenceqc/current/"
cancer_dir = "/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/cancer/current/"
fh_dir = "/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/family_history/"
bio_dir = "/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/adjusted_biomarkers/"
iop_dir = "/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/iop/"

for dirname in cbi_dir, hc_dir, cancer_dir, fh_dir, bio_dir, iop_dir:
    do_symlinks(dirname)
