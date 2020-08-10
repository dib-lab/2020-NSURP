## NSURP IBD Analysis ##
import pandas

SAMPLES = ["SRR1211157", "SRR1211428", "SRR1211440", "SRR1211568", "SRR1757110", "SRR1765025"]

#samples_csv= "ibd_samples.csv"
# to do: read the samples from csv instead

out_dir = "outputs"
logs_dir = os.path.join(out_dir, "logs")

rule all:
    input:
        expand(os.path.join(out_dir, "gather", "{alphabet}", "{sample}_x_genbank-k{ksize}.csv"), sample=SAMPLES, ksize=[21,31,51], alphabet="dna"), #[21,31,51])
        expand(os.path.join(out_dir, "gather", "{alphabet}", "{sample}_x_gtdb_pep.rep_genus-k{ksize}.csv"), sample=SAMPLES, alphabet="protein", ksize=33),
        expand(os.path.join(out_dir, "gather", "{alphabet}", "{sample}_x_gtdb_pep-k{ksize}.csv"), sample=SAMPLES, alphabet="protein", ksize=33),
        #expand(os.path.join(out_dir, "gather_to_tax", "{alphabet}", "{sample}_x_gtdb_pep.rep_genus-k{ksize}.gather_summary.csv"), sample=SAMPLES, alphabet="protein", ksize=33),
        expand(os.path.join(out_dir, "gather", "{alphabet}", "{sample}_x_gtdb_pep.rep_genus-k{ksize}.csv"), sample=SAMPLES, alphabet="dayhoff", ksize=57),
        expand(os.path.join(out_dir, "gather", "{alphabet}", "{sample}_x_gtdb_pep-k{ksize}.csv"), sample=SAMPLES, alphabet="dayhoff", ksize=57),
        #expand(os.path.join(out_dir, "gather_to_tax", "{alphabet}", "{sample}_x_gtdb_pep.rep_genus-k{ksize}.gather_summary.csv"), sample=SAMPLES, alphabet="dayhoff", ksize=57),
        #expand(os.path.join(out_dir, "compare","{sample}.{alphabet}-k{ksize}.compare.np.matrix.pdf"), sample=SAMPLES, alphabet="dna", ksize=[21,31,51]),
        #expand(os.path.join(out_dir, "compare","{sample}.{alphabet}-k{ksize}.compare.np.matrix.pdf"), sample=SAMPLES, alphabet="protein", ksize=[33]),
        #expand(os.path.join(out_dir, "compare","{sample}.{alphabet}-k{ksize}.compare.np.matrix.pdf"), sample=SAMPLES, alphabet="dayhoff", ksize=[57]),
        

#rule download_hostseqs:
## http://seqanswers.com/forums/archive/index.php/t-42552.html
# https://drive.google.com/file/d/0B3llHR93L14wd0pSSnFULUlhcUk/edit?usp=sharing
#    output: "databases/hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz"
#    params:
#        download_link="https://osf.io/p6dg4/download"
#    log: os.path.join(logs_dir, "download_hostseqs.log")
#    shell:
#        """
#        wget {params.download_link} > {output} 2> {log}
#        """

#rule download_reads:
#    output:
#        r1=os.path.join(out_dir, "raw_data", "{sample}_1.fastq.gz"),
#        r2=os.path.join(out_dir, "raw_data", "{sample}_2.fastq.gz"),
#    params:
#        download_link = lambda wildcards: samples_csv[wildcards.sample][""]
#    shell:
#        """
#        wget {params.download_link} > {output}
#        """

rule trim_reads:
    input:
        in1=os.path.join(out_dir, "raw_data", "{sample}_1.fastq.gz"),
        in2=os.path.join(out_dir, "raw_data", "{sample}_2.fastq.gz"),
    output:
        out1=os.path.join(out_dir, "trim", "{sample}_1.trim.fastq.gz"),
        out2=os.path.join(out_dir, "trim", "{sample}_2.trim.fastq.gz"),
        json=os.path.join(out_dir, "trim", "{sample}.fastp.json"),
        html=os.path.join(out_dir, "trim", "{sample}.fastp.html")
    resources:
        mem_mb=lambda wildcards, attempt: attempt *10000,
        runtime=600,
    shell:
        """
        fastp --in1 {input.in1}  --in2 {input.in2}  \
        --out1 {output.out1} --out2 {output.out2}  \
        --detect_adapter_for_pe  --qualified_quality_phred 4 \
        --length_required 31 --correction \
        --json {output.json} --html {output.html}
        """

rule remove_host:
    input:
        r1=os.path.join(out_dir, "trim", "{sample}_1.trim.fastq.gz"),
        r2=os.path.join(out_dir, "trim", "{sample}_2.trim.fastq.gz"),
        human='databases/hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz'
    output:
        r1 = os.path.join(out_dir, "bbduk", "{sample}_1.nohost.fq.gz"),
        r2 = os.path.join(out_dir, "bbduk", "{sample}_2.nohost.fq.gz"),
        human_r1=os.path.join(out_dir, "bbduk", "{sample}_1.human.fq.gz"),
        human_r2=os.path.join(out_dir, "bbduk", "{sample}_2.human.fq.gz")
    conda: 'envs/bbduk-env.yml'
    resources:
        mem_mb=lambda wildcards, attempt: attempt *64000,
        runtime=600,
    shell:
        """
        bbduk.sh -Xmx64g t=3 in={input.r1} in2={input.r2} out={output.r1} out2={output.r2} outm={output.human_r1} outm2={output.human_r2} k=31 ref={input.human}
        """

rule kmer_trim:
    input:
        r1 = os.path.join(out_dir, "bbduk", "{sample}_1.nohost.fq.gz"),
        r2 = os.path.join(out_dir, "bbduk", "{sample}_2.nohost.fq.gz"),
    output:
        os.path.join(out_dir, "kmer-trim", "{sample}.nohost.kmer-trim.pe.fastq.gz"),
    conda: 'envs/bbduk-env.yml'
    resources:
        mem_mb=lambda wildcards, attempt: attempt *5000,
        runtime=600,
    shell:
        """
        interleave-reads.py {input.r1} {input.r2} |
        trim-low-abund.py --gzip -C 3 -Z 18 -M 20e9 -V - -o {output}
        """

ksizesD={"dna": "21,31,51", "protein": "33", "dayhoff": "57"}
scaledD={"dna": "2000", "protein": "100", "dayhoff": "100"}

rule sourmash_compute:
    input: rules.kmer_trim.output
    output: os.path.join(out_dir, "sourmash_signatures", "{alphabet}", "{sample}.nohost.kmer-trim.pe.sig")
    params:
        scaled= lambda w: scaledD[w.alphabet], #2000,
        k=lambda w: ksizesD[w.alphabet], #"21,31,51,33,57",
        alpha_cmd= lambda w: "--" + w.alphabet, # --dna --protein --dayhoff",
        abund_cmd= "--track-abundance"
    log: os.path.join(logs_dir, "sourmash", "{alphabet}", "{sample}.nohost.kmer-trim.pe.compute.log")
    conda: "envs/sourmash-env.yml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt *5000,
        runtime=600,
    shell:
        """
        sourmash compute -k {params.k} --scaled={params.scaled}  \
        {input} -o {output} {params.alpha_cmd} {params.abund_cmd} --merge={wildcards.sample:q} 2> {log}
        """

rule sourmash_compare:
    input: expand(os.path.join(out_dir, "sourmash_signatures", "{{alphabet}}", "{sample}.nohost.kmer-trim.pe.sig"), sample=SAMPLES)
    output:
        csv = os.path.join(out_dir, "compare", "{sample}.{alphabet}-k{ksize}.compare.csv"),
        np = os.path.join(out_dir, "compare", "{sample}.{alphabet}-k{ksize}.compare.np")
    params:
        alpha_cmd = lambda wildcards: "--" + wildcards.alphabet 
    conda: "envs/sourmash-env.yml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt *5000,
        runtime=60,
    shell:
        """
        sourmash compare -k {wildcards.ksize} -o {output.np} --csv {output.csv} --ignore-abundance {params.alpha_cmd}
        """

localrules: sourmash_plot

rule sourmash_plot:
    input: rules.sourmash_compare.output.np
    output:
        matrix = os.path.join(out_dir, "compare","{sample}.{alphabet}-k{ksize}.compare.np.matrix.pdf")
    conda: "envs/sourmash-env.yml"
    #resources:
    #    mem_mb=lambda wildcards, attempt: attempt *1000,
    #    runtime=60,
    shell:
        """
        sourmash plot --labels {input}
        """

rule sourmash_gather_genbank:
    input: 
        sig=rules.sourmash_compute.output,
        ref="databases/genbank-k{ksize}.sbt.zip"
    output:
        csv=os.path.join(out_dir, "gather", "{alphabet}", "{sample}_x_genbank-k{ksize}.csv")
    conda: "envs/sourmash-env.yml"
    params:
        #threshold="10000000",
        alpha_cmd=lambda wildcards: "--" + wildcards.alphabet
    resources:
        mem_mb=lambda wildcards, attempt: attempt *10000,
        runtime=600,
    shell:
        """
        sourmash gather {input.sig} {input.ref} -o {output.csv} {params.alpha_cmd} -k {wildcards.ksize}
        """
        #  --threshold-bp {params.threshold}


pepD_full={"protein": "databases/gtdb_pep.protein_scaled100_k11.index.sbt.zip", "dayhoff":"databases/gtdb_pep.dayhoff_scaled100_k19.index.sbt.zip"}
pepD={"protein": "databases/gtdb_pep.rep_genus.protein_scaled100_k11.sbt.zip", "dayhoff":"databases/gtdb_pep.rep_genus.dayhoff_scaled100_k19.sbt.zip"}

rule gather_gtdb_pep_rep:
    input:
        query=os.path.join(out_dir, "sourmash_signatures", "{alphabet}", "{sample}.nohost.kmer-trim.pe.sig"),
        db=lambda w: pepD[w.alphabet]
    output:
        csv = os.path.join(out_dir, "gather", "{alphabet}", "{sample}_x_gtdb_pep.rep_genus-k{ksize}.csv"),
        matches = os.path.join(out_dir, "gather", "{alphabet}", "{sample}_x_gtdb_pep.rep_genus-k{ksize}.matches"),
        unassigned = os.path.join(out_dir, "gather", "{alphabet}", "{sample}_x_gtdb_pep.rep_genus-k{ksize}.unassigned")
    params:
        alpha_cmd = lambda w: "--" + w.alphabet,
        scaled = 100,
    resources:
        mem_mb=lambda wildcards, attempt: attempt *10000,
        runtime=600,
    log: os.path.join(logs_dir, "gather", "{alphabet}-k{ksize}",  "{sample}_x_gtdb_pep.rep_genus.{alphabet}-k{ksize}.gather.log")
    benchmark: os.path.join(logs_dir, "gather","{alphabet}-k{ksize}", "{sample}_x_gtdb_pep.rep_genus.{alphabet}-k{ksize}.gather.benchmark")
    conda: "envs/sourmash-env.yml"
    shell:
        # --ignore-abundance to turn abund off
        """
        sourmash gather {input.query} {input.db} -o {output.csv} {params.alpha_cmd} \
        --save-matches {output.matches} \
        --output-unassigned {output.unassigned} \
        --scaled {params.scaled} \
        -k {wildcards.ksize} 2> {log}
        touch {output}
        """

rule gather_gtdb_pep_full:
    input:
        query=os.path.join(out_dir, "sourmash_signatures", "{alphabet}", "{sample}.nohost.kmer-trim.pe.sig"),
        db=lambda w: pepD_full[w.alphabet]
    output:
        csv = os.path.join(out_dir, "gather", "{alphabet}", "{sample}_x_gtdb_pep-k{ksize}.csv"),
        matches = os.path.join(out_dir, "gather", "{alphabet}", "{sample}_x_gtdb_pep-k{ksize}.matches"),
        unassigned = os.path.join(out_dir, "gather", "{alphabet}", "{sample}_x_gtdb_pep-k{ksize}.unassigned")
    params:
        alpha_cmd = lambda w: "--" + w.alphabet,
        scaled = 100,
    resources:
        mem_mb=lambda wildcards, attempt: attempt *10000,
        runtime=600,
    log: os.path.join(logs_dir, "gather", "{alphabet}-k{ksize}",  "{sample}_x_gtdb_pep.{alphabet}-k{ksize}.gather.log")
    benchmark: os.path.join(logs_dir, "gather","{alphabet}-k{ksize}", "{sample}_x_gtdb_pep.{alphabet}-k{ksize}.gather.benchmark")
    conda: "envs/sourmash-env.yml"
    shell:
        # --ignore-abundance to turn abund off
        """
        sourmash gather {input.query} {input.db} -o {output.csv} {params.alpha_cmd} \
        --save-matches {output.matches} \
        --output-unassigned {output.unassigned} \
        --scaled {params.scaled} \
        -k {wildcards.ksize} 2> {log}
        touch {output}
        """

rule gather_to_tax_pep_rep:
    input:
        gather_csv = os.path.join(out_dir, "gather", "{alphabet}", "{sample}_x_gtdb_pep.rep_genus-k{ksize}.csv"),
        lineages_csv = "databases/gtdb-lineages.protein-filenames.representative-at-genus.csv"
    output:
        gather_tax = os.path.join(out_dir, "gather_to_tax", "{alphabet}", "{sample}_x_gtdb_pep.rep_genus-k{ksize}.gather_summary.csv"),
        top_matches = os.path.join(out_dir, "gather_to_tax", "{alphabet}", "{sample}_x_gtdb_pep.rep_genus-k{ksize}.gather_tophits.csv"),
    log: os.path.join(logs_dir, "gather_to_tax", "{alphabet}-k{ksize}", "{sample}_x_gtdb_pep.rep_genus.{alphabet}-k{ksize}.gather-to-tax.log")
    benchmark: os.path.join(logs_dir, "gather_to_tax", "{alphabet}-k{ksize}","{sample}_x_gtdb_pep.rep_genus.{alphabet}-k{ksize}.gather-to-tax.benchmark")
    group: "gather"
    resources:
        mem_mb=lambda wildcards, attempt: attempt *10000,
        runtime=200,
    conda: "envs/sourmash-env.yml"
    shell:
        """
        python scripts/gather-to-tax.py {input.gather_csv} {input.lineages_csv} --tophits-csv {output.top_matches} > {output.gather_tax} 2> {log}
        """


#rule gather_to_tax_dna:
#    input:
#        gather_csv = os.path.join(out_dir, "gather", "{sample}_x_genbank-k{ksize}.csv"
#        #lineages_csv = lambda w: refInfo[moltype_map[w.alphabet]][w.db_name]["lineages_csv"]
#    output:
#        gather_tax = os.path.join(gather_dir, "{db_name}.{alphabet}-k{ksize}", "{genome}_x_{db_name}.{alphabet}-k{ksize}.gather_summary.csv"),
#        top_matches = os.path.join(gather_dir, "{db_name}.{alphabet}-k{ksize}", "{genome}_x_{db_name}.{alphabet}-k{ksize}.gather_tophits.csv")
#    log: os.path.join(logs_dir, "gather_to_tax", "{db_name}.{alphabet}-k{ksize}", input_type, "{genome}_x_{db_name}.{alphabet}-k{ksize}.gather-to-tax.log")
#    benchmark: os.path.join(logs_dir, "gather_to_tax", "{db_name}.{alphabet}-k{ksize}","{genome}_x_{db_name}.{alphabet}-k{ksize}.gather-to-tax.benchmark")
#    resources:
#        mem_mb=lambda wildcards, attempt: attempt *10000,
#        runtime=200,
#    conda: "envs/sourmash-env.yml"
#    shell:
#        """
#        python scripts/gather-to-tax.py {input.gather_csv} {input.lineages_csv} --tophits-csv {output.top_matches} > {output.gather_tax} 2> {log}
#        """

#rule megahit_assemble:
#    input:
#        os.path.join(out_dir, "kmer-trim", "{sample}.nohost.kmer-trim.pe.fastq.gz") 
#    output:
#        os.path.join(out_dir, "megahit", "{sample}_contigs.fa" 
#    message:
#        """### Assembling read data with MEGAHIT ### """
#    params:
#        megahit_dir=os.path.join(out_dir, "megahit")
#    threads: 10
#    resources:
#        mem_mb=20000,
#        runtime=600
#    log: os.path.join(logs_dir, "megahit", "{sample}_megahit.log")
#    benchmark: os.path.join(logs_dir, "megahit", "{sample}_megahit.benchmark")
#    conda: "envs/megahit-env.yaml"
#    shell:
#        """
#        megahit --12 {input} -o {params.megahit_dir} --out_prefix {wildcards.sample} -t {threads}
#        """
