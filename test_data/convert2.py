# transfer file to vcf format
import tskit
import gzip
for j in range(10):
    ts=tskit.load("./DEC_cor/dec{}".format(j))
    n_dip_indv = int(ts.num_samples / 2)
    indv_names = [f"tsk_{i}indv" for i in range(n_dip_indv)]


    with gzip.open("./DEC_cor/cdec{}.vcf.gz".format(j), "wt") as vcf_file:
        ts.write_vcf(vcf_file, individual_names=indv_names,contig_id=j)

with open("./DEC_cor/cdec.popLabel","w") as f:
    f.write("sample_id\tpopulation\n")
    for i in ts.individuals():
        f.write("{}\t{}\n".format(f"tsk_{i.id}indv",i.population))