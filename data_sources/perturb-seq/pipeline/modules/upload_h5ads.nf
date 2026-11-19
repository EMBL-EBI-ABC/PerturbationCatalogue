process UPLOAD_H5ADS {
    tag "${dataset_id}"
    container null
    cpus 2
    memory '4 GB'
    time '12h'
    errorStrategy 'terminate'

    input:
    tuple val(dataset_id), path(raw_h5ad), path(filtered_h5ad), path(summary)

    output:
    val dataset_id

    script:
    def destination = params.h5ad_release_uri.replaceAll('/+$', '')
    """
    gcloud storage cp '${raw_h5ad}' '${destination}/${dataset_id}.h5ad'
    gcloud storage cp '${filtered_h5ad}' '${destination}/${dataset_id}.filtered.h5ad'
    """
}
