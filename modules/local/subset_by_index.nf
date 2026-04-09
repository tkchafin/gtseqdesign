process SUBSET_BY_INDEX {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5acacb55c52bec97c61fd34ffa8721fce82ce823005793592e2a80bf71632cd0/data' :
        'community.wave.seqera.io/library/bcftools:1.21--4335bec1d7b44d11' }"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(tbi)
    tuple val(meta3), path(top_loci)

    output:
    tuple val(meta), path("${meta.id}.selected.vcf.gz"), emit: vcf
    tuple val(meta), path("${meta.id}.selected.vcf.gz.tbi"), emit: tbi, optional: true
    tuple val(meta), path("targets.txt"), emit: targets
    path "versions.yml", emit: versions

    script:
    """
    echo "📋 Mapping VCF records to CHROM:POS..."
    bcftools query -f '%CHROM\t%POS\\n' ${vcf} > all_sites.txt

    echo "🧮 Selecting top loci from top_loci.txt..."
    awk 'NR==FNR{keep[\$1]; next} (FNR) in keep' ${top_loci} all_sites.txt > targets.txt

    echo "🔍 Filtering VCF using bcftools view with targets..."
    bcftools view \\
        --targets-file targets.txt \\
        --output-type z \\
        --output ${meta.id}.selected.vcf.gz \\
        --threads ${task.cpus} \\
        ${vcf}

    echo "📦 Indexing selected VCF..."
    bcftools index ${meta.id}.selected.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
