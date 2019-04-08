import os
from annotate_phe import make_phe_info
import glob

if __name__ == "__main__":
    with open('../../../../wiki/ukbb/icdinfo/icdinfo.txt', 'r') as icd:
        icdinfo = {line.split()[0]:line.split()[2] for line in icd}
    with open('activity_filename_to_gbe_id.tsv', 'r') as new_ids:
        ac_gbe_map = {line.split()[1]:line.split()[0] for line in new_ids}
    activity_phenos = glob.glob('/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/physical_activity/phe/*.phe')
    
    for pheno in activity_phenos:
        ac = os.path.splitext(os.path.basename(pheno))[0] 
        if ac in ac_gbe_map:
            gbe_id = ac_gbe_map[ac]
            os.rename(pheno, '/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/physical_activity/phe/' + gbe_id + '.phe')
    with open('activity_filename_to_gbe_id.tsv', 'r') as new_ids:
        gbe_ac_map = {line.split()[0]:line.split()[1] for line in new_ids}
    
    activity_phenos = glob.glob('/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/physical_activity/phe/*.phe')
    activity_phenos = [x for x in activity_phenos if (os.path.splitext(os.path.basename(x))[0].startswith('HC') or os.path.splitext(os.path.basename(x))[0].startswith('INI'))]
    activity_names  = []

    for pheno in activity_phenos:
        gbe_id = os.path.splitext(os.path.basename(pheno))[0]
        if gbe_id in icdinfo:
            activity_names.append(icdinfo[gbe_id])
        elif gbe_id in gbe_ac_map:
            activity_names.append(gbe_ac_map[gbe_id])
    
    make_phe_info(in_phe   = activity_phenos,
                  out_path = '/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/physical_activity/info/',
                  name     = activity_names,
                  field    = 'NA',
                  table    = 'NA',
                  basket   = 'NA',
                  app_id   = '24983',
                  source   = 'ashley_lab',
                  one_file = False)
