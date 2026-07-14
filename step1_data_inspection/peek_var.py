
import h5py
base="/storage/home/hcoda1/9/ykim3030/scratch/hackathon/0-Data"
for tag,p in [("DE", base+"/perturbseq_marson2025/GWCD4i.DE_stats.h5ad"),
              ("COSMX", base+"/ST_CosMx_NSCLC/cosmx_human_nsclc_clustered.h5ad")]:
    print("====",tag,p)
    with h5py.File(p,"r") as f:
        print("root keys:", list(f.keys()))
        for grp in ["var","obs"]:
            if grp in f:
                v=f[grp]
                print(grp,"keys:", list(v.keys()) if hasattr(v,"keys") else str(v))
                print(grp,"attrs:", dict(v.attrs))
