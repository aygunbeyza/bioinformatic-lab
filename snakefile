rule all:
    input:"Seurat.rds"




rule multi2:
    input:
        h5= "GSM3489182_Donor_01_filtered_gene_bc_matrices_h5.h5"
    output:
        out= "Seurat.rds"
    shell:
        """
        Rscript deneme.R\
        {input.h5} \
        {output.out}
        """

#"Sample.Process_code.rds"
