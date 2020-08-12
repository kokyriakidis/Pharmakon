from sys import platform as _platform
genome_build = config["ref"]["build"]

if config["ref"]["build"] == "hg38":

    rule get_genome:
        output:
            fasta = "{output_dir}/{genome_build}/seq/{genome_build}.fa".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
            fasta_fai = "{output_dir}/{genome_build}/seq/{genome_build}.fa.fai".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
            fasta_gz = "{output_dir}/{genome_build}/seq/{genome_build}.fa.gz".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
            fasta_gz_fai = "{output_dir}/{genome_build}/seq/{genome_build}.fa.gz.fai".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
            fasta_gz_gzi = "{output_dir}/{genome_build}/seq/{genome_build}.fa.gz.gzi".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
            fasta_dict = "{output_dir}/{genome_build}/seq/{genome_build}.dict".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build)
        log:
            "{output_dir}/logs/reference/{genome_build}/get_{genome_build}_genome.log".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build)
        params:
            output_dir=config["ref"]["output_dir"],
            genome_build=config["ref"]["build"]
        cache: True
        shell:
            """
            url=http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome
            base=GRCh38_full_analysis_set_plus_decoy_hla
            new=hg38
            mkdir -p {params.output_dir}/{params.genome_build}/seq
            for suffix in .dict .fa.fai
            do
                [[ -f {params.output_dir}/{params.genome_build}/seq/$new$suffix ]] || wget -c -O {params.output_dir}/{params.genome_build}/seq/$new$suffix $url/$base$suffix > {log} 2>&1
            done
            # gzipped references for use with ensembl-vep HGVS
            urlgz=https://s3.amazonaws.com/biodata/hg38_bundle
            for suffix in .fa.gz .fa.gz.fai .fa.gz.gzi
            do
                [[ -f {params.output_dir}/{params.genome_build}/seq/$new$suffix ]] || wget --no-check-certificate -c -O {params.output_dir}/{params.genome_build}/seq/$new$suffix $urlgz/$new$suffix > {log} 2>&1
            done
            gunzip -c {params.output_dir}/{params.genome_build}/seq/hg38.fa.gz > {params.output_dir}/{params.genome_build}/seq/hg38.fa
            touch {params.output_dir}/{params.genome_build}/seq/hg38.fa.fai
            touch {params.output_dir}/{params.genome_build}/seq/hg38.dict
            """

    rule get_dbsnp:
        input:
            fasta = "{output_dir}/{genome_build}/seq/{genome_build}.fa".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build)
        output:
            dbsnp = "{output_dir}/{genome_build}/variation/dbsnp-153.vcf.gz".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
            dbsnp_csi = "{output_dir}/{genome_build}/variation/dbsnp-153.vcf.gz.csi".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
            dbsnp_tbi = "{output_dir}/{genome_build}/variation/dbsnp-153.vcf.gz.tbi".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build)
        log:
            "{output_dir}/logs/reference/{genome_build}/get_{genome_build}_dbsnp-153.log".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build)
        params:
            output_dir=config["ref"]["output_dir"],
            genome_build=config["ref"]["build"]
        cache: True
        shell:
            """
            build=153
            version=GCF_000001405.38
            url=http://ftp.ncbi.nih.gov/snp/archive/b$build/VCF/$version.gz
            remap_url=https://gist.githubusercontent.com/matthdsm/f833aedd2d67e28013ff1d171c70f4ee/raw/442a45ed3ddc6e85c66c5e58e0fa78e16a0821c8/refseq2ucsc.tsv
            ref={params.output_dir}/{params.genome_build}/seq/hg38.fa
            mkdir -p {params.output_dir}/{params.genome_build}/variation
            wget -c -O {params.output_dir}/{params.genome_build}/variation/dbsnp-$build-orig.vcf.gz $url
            wget -c -O {params.output_dir}/{params.genome_build}/variation/dbsnp-$build-orig.vcf.gz.tbi $url.tbi
            [[ -f {params.output_dir}/{params.genome_build}/variation/dbsnp-$build.vcf.gz ]] || bcftools annotate -Ou --rename-chrs $remap_url {params.output_dir}/{params.genome_build}/variation/dbsnp-$build-orig.vcf.gz |\
            bcftools sort -m 1G -Oz -T . -o {params.output_dir}/{params.genome_build}/variation/dbsnp-$build.vcf.gz && \
            tabix -f -p vcf -C {params.output_dir}/{params.genome_build}/variation/dbsnp-$build.vcf.gz
            tabix -f -p vcf {params.output_dir}/{params.genome_build}/variation/dbsnp-$build.vcf.gz
            """

elif config["ref"]["build"] == "hg19":

    rule get_genome:
        output:
            fasta = "{output_dir}/{genome_build}/seq/{genome_build}.fa".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
            fasta_fai = "{output_dir}/{genome_build}/seq/{genome_build}.fa.fai".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
            fasta_gz = "{output_dir}/{genome_build}/seq/{genome_build}.fa.gz".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
            fasta_gz_fai = "{output_dir}/{genome_build}/seq/{genome_build}.fa.gz.fai".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
            fasta_gz_gzi = "{output_dir}/{genome_build}/seq/{genome_build}.fa.gz.gzi".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
            fasta_dict = "{output_dir}/{genome_build}/seq/{genome_build}.dict".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build)
        log:
            "{output_dir}/logs/reference/{genome_build}/get_{genome_build}_genome.log".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build)
        params:
            output_dir=config["ref"]["output_dir"],
            genome_build=config["ref"]["build"]
        cache: True
        shell:
            """
            mkdir -p {params.output_dir}/{params.genome_build}
            wget --no-check-certificate -c -O {params.output_dir}/{params.genome_build}/hg19-seq.tar.gz https://s3.amazonaws.com/biodata/genomes/hg19-seq.tar.gz > {log} 2>&1
            tar -xzvpf {params.output_dir}/{params.genome_build}/hg19-seq.tar.gz -C {params.output_dir}/{params.genome_build}
            gunzip -c {params.output_dir}/{params.genome_build}/seq/hg19.fa.gz > {params.output_dir}/{params.genome_build}/seq/hg19.fa
            touch {params.output_dir}/{params.genome_build}/seq/hg19.fa.fai
            touch {params.output_dir}/{params.genome_build}/seq/hg19.dict
            #mv -v {params.output_dir}/{params.genome_build}/seq/* {params.output_dir}/{params.genome_build}
            #rm -rf {params.output_dir}/{params.genome_build}/seq
            rm {params.output_dir}/{params.genome_build}/hg19-seq.tar.gz
            """

    rule get_dbsnp:
        input:
            fasta = "{output_dir}/{genome_build}/seq/{genome_build}.fa".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build)
        output:
            dbsnp = "{output_dir}/{genome_build}/variation/dbsnp-151.vcf.gz".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
            dbsnp_tbi = "{output_dir}/{genome_build}/variation/dbsnp-151.vcf.gz.tbi".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build)
        log:
            "{output_dir}/logs/reference/{genome_build}/get_{genome_build}_dbsnp-153.log".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build)
        params:
            output_dir=config["ref"]["output_dir"],
            genome_build=config["ref"]["build"]
        cache: True
        shell:
            """
            build=151
            org=human_9606_b${version}_GRCh37p13            
            release=20180423
            url=http://ftp.ncbi.nih.gov/snp/pre_build152/organisms/$org/VCF/GATK/All_${release}.vcf.gz
            mkdir -p {params.output_dir}/{params.genome_build}/variation
            wget -c -O {params.output_dir}/{params.genome_build}/variation/dbsnp-$version-orig.vcf.gz $url
            [[ -f {params.output_dir}/{params.genome_build}/variation/dbsnp-$version.vcf.gz ]] || zcat {params.output_dir}/{params.genome_build}/variation/dbsnp-$version-orig.vcf.gz | bgzip -c > {params.output_dir}/{params.genome_build}/variation/dbsnp-$version.vcf.gz
            tabix -f -p vcf {params.output_dir}/{params.genome_build}/variation/dbsnp-$version.vcf.gz
            """


if _platform == "darwin" or config["variant_tool"] == "gatk":

    if config["ref"]["build"] == "hg38":

        rule get_hapmap:
            output:
                hapmap = "{output_dir}/{genome_build}/variation/hapmap_3.3.vcf.gz".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
                hapmap_gz = "{output_dir}/{genome_build}/variation/hapmap_3.3.vcf.gz.tbi".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build)
            log:
                "{output_dir}/logs/reference/{genome_build}/get_{genome_build}_hapmap.log".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build)
            params:
                output_dir=config["ref"]["output_dir"],
                genome_build=config["ref"]["build"]
            cache: True
            shell:
                """
                url=ftp://gsapubftp-anonymous:none@ftp.broadinstitute.org/bundle/hg38/
                new=hapmap_3.3
                base=$new.hg38
                mkdir -p {params.output_dir}/{params.genome_build}/variation
                for suffix in .vcf.gz .vcf.gz.tbi
                do
                    [[ -f {params.output_dir}/{params.genome_build}/variation/$new$suffix ]] || wget --no-check-certificate -c -O {params.output_dir}/{params.genome_build}/variation/$new$suffix $url/$base$suffix
                done
                """

        rule get_1000g_omni_snps:
            output:
                omni_snps = "{output_dir}/{genome_build}/variation/1000G_omni2.5.vcf.gz".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
                omni_snps_gz = "{output_dir}/{genome_build}/variation/1000G_omni2.5.vcf.gz.tbi".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build)
            log:
                "{output_dir}/logs/reference/{genome_build}/get_{genome_build}_1000G_omni_snps.log".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build)
            params:
                output_dir=config["ref"]["output_dir"],
                genome_build=config["ref"]["build"]
            cache: True
            shell:
                """
                url=ftp://gsapubftp-anonymous:none@ftp.broadinstitute.org/bundle/hg38/
                new=1000G_omni2.5
                base=$new.hg38
                mkdir -p {params.output_dir}/{params.genome_build}/variation
                for suffix in .vcf.gz .vcf.gz.tbi
                do
                    [[ -f {params.output_dir}/{params.genome_build}/variation/$new$suffix ]] || wget --no-check-certificate -c -O {params.output_dir}/{params.genome_build}/variation/$new$suffix $url/$base$suffix
                done
                """

        rule get_1000g_snps:
            output:
                onek_snps = "{output_dir}/{genome_build}/variation/1000G_phase1.snps.high_confidence.vcf.gz".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
                onek_snps_gz = "{output_dir}/{genome_build}/variation/1000G_phase1.snps.high_confidence.vcf.gz.tbi".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build)
            log:
                "{output_dir}/logs/reference/{genome_build}/get_{genome_build}_1000G_phase1.snps.high_confidence.log".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build)
            params:
                output_dir=config["ref"]["output_dir"],
                genome_build=config["ref"]["build"]
            cache: True
            shell:
                """
                url=ftp://gsapubftp-anonymous:none@ftp.broadinstitute.org/bundle/hg38/
                new=1000G_phase1.snps.high_confidence
                base=$new.hg38
                mkdir -p {params.output_dir}/{params.genome_build}/variation
                for suffix in .vcf.gz .vcf.gz.tbi
                do
                    [[ -f {params.output_dir}/{params.genome_build}/variation/$new$suffix ]] || wget --no-check-certificate -c -O {params.output_dir}/{params.genome_build}/variation/$new$suffix $url/$base$suffix
                done
                """

        rule get_mills_indels:
            output:
                mills_indels = "{output_dir}/{genome_build}/variation/Mills_and_1000G_gold_standard.indels.vcf.gz".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
                mills_indels_gz = "{output_dir}/{genome_build}/variation/Mills_and_1000G_gold_standard.indels.vcf.gz.tbi".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build)
            log:
                "{output_dir}/logs/reference/{genome_build}/get_{genome_build}_Mills_and_1000G_gold_standard.indels.log".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build)
            params:
                output_dir=config["ref"]["output_dir"],
                genome_build=config["ref"]["build"]
            cache: True
            shell:
                """
                url=ftp://gsapubftp-anonymous:none@ftp.broadinstitute.org/bundle/hg38/
                new=Mills_and_1000G_gold_standard.indels
                base=$new.hg38
                mkdir -p {params.output_dir}/{params.genome_build}/variation
                for suffix in .vcf.gz .vcf.gz.tbi
                do
                    [[ -f {params.output_dir}/{params.genome_build}/variation/$new$suffix ]] || wget --no-check-certificate -c -O {params.output_dir}/{params.genome_build}/variation/$new$suffix $url/$base$suffix
                done
                """


    elif config["ref"]["build"] == "hg19":

        rule get_hapmap:
            output:
                hapmap = "{output_dir}/{genome_build}/variation/hapmap_3.3.vcf.gz".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
                hapmap_gz = "{output_dir}/{genome_build}/variation/hapmap_3.3.vcf.gz.tbi".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build)
            log:
                "{output_dir}/logs/reference/{genome_build}/get_{genome_build}_hapmap.log".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build)
            params:
                output_dir=config["ref"]["output_dir"],
                genome_build=config["ref"]["build"]
            cache: True
            shell:
                """
                baseurl=ftp://gsapubftp-anonymous:none@ftp.broadinstitute.org/bundle/hg19/hapmap_3.3.hg19.sites.vcf.gz
                mkdir -p {params.output_dir}/{params.genome_build}/variation
                cd {params.output_dir}/{params.genome_build}/variation
                wget -O - $baseurl | gunzip -c | bgzip -c > hapmap_3.3.vcf.gz
                tabix -f -p vcf hapmap_3.3.vcf.gz
                """

        rule get_1000g_omni_snps:
            output:
                omni_snps = "{output_dir}/{genome_build}/variation/1000G_omni2.5.vcf.gz".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
                omni_snps_gz = "{output_dir}/{genome_build}/variation/1000G_omni2.5.vcf.gz.tbi".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build)
            log:
                "{output_dir}/logs/reference/{genome_build}/get_{genome_build}_1000G_omni_snps.log".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build)
            params:
                output_dir=config["ref"]["output_dir"],
                genome_build=config["ref"]["build"]
            cache: True
            shell:
                """
                baseurl=ftp://gsapubftp-anonymous:none@ftp.broadinstitute.org/bundle/hg19/1000G_omni2.5.hg19.sites.vcf.gz
                mkdir -p {params.output_dir}/{params.genome_build}/variation
                cd {params.output_dir}/{params.genome_build}/variation
                wget -O - $baseurl | gunzip -c | bgzip -c > 1000G_omni2.5.vcf.gz
                tabix -f -p vcf 1000G_omni2.5.vcf.gz
                """

        rule get_1000g_snps:
            output:
                onek_snps = "{output_dir}/{genome_build}/variation/1000G_phase1.snps.high_confidence.vcf.gz".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
                onek_snps_gz = "{output_dir}/{genome_build}/variation/1000G_phase1.snps.high_confidence.vcf.gz.tbi".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build)
            log:
                "{output_dir}/logs/reference/{genome_build}/get_{genome_build}_1000G_phase1.snps.high_confidence.log".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build)
            params:
                output_dir=config["ref"]["output_dir"],
                genome_build=config["ref"]["build"]
            cache: True
            shell:
                """
                baseurl=ftp://gsapubftp-anonymous:none@ftp.broadinstitute.org/bundle/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz
                mkdir -p {params.output_dir}/{params.genome_build}/variation
                cd {params.output_dir}/{params.genome_build}/variation
                wget -O - $baseurl | gunzip -c | bgzip -c > 1000G_phase1.snps.high_confidence.vcf.gz
                tabix -f -p vcf 1000G_phase1.snps.high_confidence.vcf.gz
                """

        rule get_mills_indels:
            output:
                mills_indels = "{output_dir}/{genome_build}/variation/Mills_and_1000G_gold_standard.indels.vcf.gz".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
                mills_indels_gz = "{output_dir}/{genome_build}/variation/Mills_and_1000G_gold_standard.indels.vcf.gz.tbi".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build)
            log:
                "{output_dir}/logs/reference/{genome_build}/get_{genome_build}_Mills_and_1000G_gold_standard.indels.log".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build)
            params:
                output_dir=config["ref"]["output_dir"],
                genome_build=config["ref"]["build"]
            cache: True
            shell:
                """
                baseurl=ftp://gsapubftp-anonymous:none@ftp.broadinstitute.org/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
                mkdir -p {params.output_dir}/{params.genome_build}/variation
                cd {params.output_dir}/{params.genome_build}/variation
                wget -O - $baseurl | gunzip -c | bgzip -c > Mills_and_1000G_gold_standard.indels.vcf.gz
                tabix -f -p vcf Mills_and_1000G_gold_standard.indels.vcf.gz
                """

    
if _platform == "darwin":    
    
    rule bwa_index:
        input:
            "{output_dir}/{genome_build}/seq/{genome_build}.fa".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build)
        output:
            "{output_dir}/{genome_build}/bwa/{genome_build}.fa.amb".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
            "{output_dir}/{genome_build}/bwa/{genome_build}.fa.ann".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
            "{output_dir}/{genome_build}/bwa/{genome_build}.fa.bwt".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
            "{output_dir}/{genome_build}/bwa/{genome_build}.fa.pac".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
            "{output_dir}/{genome_build}/bwa/{genome_build}.fa.sa".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build)
        log:
            "{output_dir}/logs/bwa_index/{genome_build}/bwa_index.log".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build)
        params:
            prefix="{output_dir}/{genome_build}/bwa/{genome_build}.fa".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
            algorithm="bwtsw"
        resources:
            mem_mb=369000
        cache: True
        wrapper:
            "0.64.0/bio/bwa/index"

    rule relocate_bwa_index_files:
        input:
            "{output_dir}/{genome_build}/seq/{genome_build}.fa.amb".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
            "{output_dir}/{genome_build}/seq/{genome_build}.fa.ann".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
            "{output_dir}/{genome_build}/seq/{genome_build}.fa.bwt".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
            "{output_dir}/{genome_build}/seq/{genome_build}.fa.pac".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
            "{output_dir}/{genome_build}/seq/{genome_build}.fa.sa".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build)        
        output:
            "{output_dir}/{genome_build}/bwa/{genome_build}.fa.amb".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
            "{output_dir}/{genome_build}/bwa/{genome_build}.fa.ann".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
            "{output_dir}/{genome_build}/bwa/{genome_build}.fa.bwt".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
            "{output_dir}/{genome_build}/bwa/{genome_build}.fa.pac".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
            "{output_dir}/{genome_build}/bwa/{genome_build}.fa.sa".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build)
        params:
            output_dir=config["ref"]["output_dir"],
            genome_build=config["ref"]["build"]
        shell:
            """
            mkdir -p {params.output_dir}/{params.genome_build}/bwa
            mv {params.output_dir}/{params.genome_build}/seq/{genome_build}.fa.amb {params.output_dir}/{params.genome_build}/bwa
            mv {params.output_dir}/{params.genome_build}/seq/{genome_build}.fa.ann {params.output_dir}/{params.genome_build}/bwa
            mv {params.output_dir}/{params.genome_build}/seq/{genome_build}.fa.bwt {params.output_dir}/{params.genome_build}/bwa
            mv {params.output_dir}/{params.genome_build}/seq/{genome_build}.fa.pac {params.output_dir}/{params.genome_build}/bwa
            mv {params.output_dir}/{params.genome_build}/seq/{genome_build}.fa.sa {params.output_dir}/{params.genome_build}/bwa
            """


elif _platform == "linux" or _platform == "linux2":

    rule bwa_mem2_index:
        input:
            "{output_dir}/{genome_build}/seq/{genome_build}.fa".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build)
        output:
            "{output_dir}/{genome_build}/bwa/{genome_build}.fa.0123".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
            "{output_dir}/{genome_build}/bwa/{genome_build}.fa.amb".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
            "{output_dir}/{genome_build}/bwa/{genome_build}.fa.ann".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
            "{output_dir}/{genome_build}/bwa/{genome_build}.fa.bwt.2bit.64".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
            "{output_dir}/{genome_build}/bwa/{genome_build}.fa.bwt.8bit.32".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
            "{output_dir}/{genome_build}/bwa/{genome_build}.fa.pac".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build)
        log:
            "{output_dir}/logs/bwa_index/{genome_build}/bwa_mem2_index.log".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build)
        params:
            prefix="{output_dir}/{genome_build}/bwa/{genome_build}.fa".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build)
        resources:
            mem_mb=369000
        cache: True
        wrapper:
            "0.64.0/bio/bwa-mem2/index"

    rule relocate_bwa_index_files:
        input:
            "{output_dir}/{genome_build}/seq/{genome_build}.fa.0123".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
            "{output_dir}/{genome_build}/seq/{genome_build}.fa.amb".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
            "{output_dir}/{genome_build}/seq/{genome_build}.fa.ann".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
            "{output_dir}/{genome_build}/seq/{genome_build}.fa.bwt.2bit.64".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
            "{output_dir}/{genome_build}/seq/{genome_build}.fa.bwt.8bit.32".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
            "{output_dir}/{genome_build}/seq/{genome_build}.fa.pac".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build)
        output:
            "{output_dir}/{genome_build}/bwa/{genome_build}.fa.0123".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
            "{output_dir}/{genome_build}/bwa/{genome_build}.fa.amb".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
            "{output_dir}/{genome_build}/bwa/{genome_build}.fa.ann".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
            "{output_dir}/{genome_build}/bwa/{genome_build}.fa.bwt.2bit.64".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
            "{output_dir}/{genome_build}/bwa/{genome_build}.fa.bwt.8bit.32".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build),
            "{output_dir}/{genome_build}/bwa/{genome_build}.fa.pac".format(output_dir=config["ref"]["output_dir"], genome_build=genome_build)
        params:
            output_dir=config["ref"]["output_dir"],
            genome_build=config["ref"]["build"]
        shell:
            """
            mkdir -p {params.output_dir}/{params.genome_build}/bwa
            mv {params.output_dir}/{params.genome_build}/seq/{genome_build}.0123 {params.output_dir}/{params.genome_build}/bwa
            mv {params.output_dir}/{params.genome_build}/seq/{genome_build}.amb {params.output_dir}/{params.genome_build}/bwa
            mv {params.output_dir}/{params.genome_build}/seq/{genome_build}.ann {params.output_dir}/{params.genome_build}/bwa
            mv {params.output_dir}/{params.genome_build}/seq/{genome_build}.bwt.2bit.64 {params.output_dir}/{params.genome_build}/bwa
            mv {params.output_dir}/{params.genome_build}/seq/{genome_build}.bwt.8bit.32 {params.output_dir}/{params.genome_build}/bwa
            mv {params.output_dir}/{params.genome_build}/seq/{genome_build}.pac {params.output_dir}/{params.genome_build}/bwa
            """