import pandas as pd

prefix = config["prefix"]

rule all:
    input:
        "{prefix}/gwas_lifted.tsv".format(prefix=prefix)
        # "{prefix}/gwas_catalogue.tsv.gz".format(prefix=config["prefix"]),
        # "{prefix}/chain/38_to_19.chain.gz".format(prefix=config["prefix"])


rule download_gwas_catalogue:
    output:
        "{prefix}/gwas_catalogue.tsv"
    shell:
        """curl https://www.ebi.ac.uk/gwas/api/search/downloads/full > {output[0]}"""
        # | tail -n +2 | cut -f 12,13,22,28 | py -x 'print(x) if float(x.split("\\t")[3]) <= 5e-8 else ""' | grep -v '^$' | awk '{{print "chr"$1, $2, $2+1, $3, $4, "."}}' | tr " " "\t" | gzip -9 > {output[0]}"""


rule parse_gwas_catalogue:
    input:
        "{prefix}/gwas_catalogue.tsv"
    output:
        "{prefix}/gwas_catalogue_parsed.tsv"
    run:
        df = pd.read_table("data/gwas_catalogue.tsv",
                           sep="\t",
                           usecols=[11,12,21,27],
                           dtype=object,
                           header=None, skiprows=1, names="Chromosome Start ID PValue".split())
        print(df.shape)
        print(df.head())
        df = df.dropna()
        print(df.shape)

        df = df[~df.Chromosome.str.contains("x|;").fillna(False)]
        df.loc[:, "ID"] = df.ID.str.replace(", ", ",")
        print(df.shape)

        df.loc[:, "PValue"] = df.PValue.astype(float)
        df.loc[:, "Start"] = df.Start.astype(int)
        df.insert(2, "End", df.Start + 1)
        df.loc[:, "Chromosome"] = "chr" + df.Chromosome
        # df = df[df.PValue <= 5e-8]

        df.to_csv(output[0], sep="\t", index=False, header=False)


rule download_chain:
    output:
        "{prefix}/chain/38_to_19.chain.gz"
    shell:
        "wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz -O {output[0]}"


        #    liftOver oldFile map.chain newFile unMapped
rule liftover:
    input:
        chain = "{prefix}/chain/38_to_19.chain.gz",
        bed = "{prefix}/gwas_catalogue_parsed.tsv"
    output:
        lifted = "{prefix}/gwas_lifted.tsv",
        unmapped = "{prefix}/gwas_unmapped.tsv"
    shell:
        "liftOver {input.bed} {input.chain} {output.lifted} {output.unmapped}"
