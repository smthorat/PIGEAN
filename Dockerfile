FROM ubuntu:plucky
RUN apt -y update && apt -y upgrade && apt -y install python3 python3-scipy
WORKDIR /app
# Reference data
ADD data/gene_set_list_mouse_2024.txt data/
ADD data/gene_set_list_msigdb_nohp.txt data/
ADD data/portal_gencode.gene.map  data/
ADD data/NCBI37.3.plink.gene.loc  data/
ADD data/refGene_hg19_TSS.subset.loc  data/
ADD data/NCBI37.3.plink.gene.exons.loc data/
# Engine (unchanged)
ADD engine/priors.py engine/priors.py
# Wrapper (Phase 1-6; ADD recursively includes pigean/modes/)
ADD run_pigean.py run_pigean.py
ADD pigean/ pigean/
ADD docs/ docs/

ENTRYPOINT ["python3", "-u", "/app/run_pigean.py"]
