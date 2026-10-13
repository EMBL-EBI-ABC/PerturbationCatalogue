-- Read-only evidence queries for the 2026-10-13 handoff deck.
-- Run against project prj-ext-dev-pertcat-437314, location europe-west2.

WITH
  wanted AS (
    SELECT dataset_id
    FROM UNNEST([
      'adamson_2016_pilot',
      'adamson_2016_upr_perturb_seq',
      'norman_2019_raw',
      'nadig_2025_hepg2',
      'nadig_2025_jurkat',
      'replogle_2022_k562_essential_normalized',
      'replogle_2022_rpe1_essential_normalized',
      'replogle_2022_k562_gw_normalized',
      'zhu_2025_D1_rest_cl',
      'zhu_2025_D1_stim48hr_cl',
      'zhu_2025_D1_stim8hr_cl',
      'zhu_2025_D2_rest_cl',
      'zhu_2025_D2_stim48hr_cl',
      'zhu_2025_D2_stim8hr_cl',
      'zhu_2025_D3_rest_cl',
      'zhu_2025_D3_stim48hr_cl',
      'zhu_2025_D3_stim8hr_cl',
      'zhu_2025_D4_rest_cl',
      'zhu_2025_D4_stim48hr_cl',
      'zhu_2025_D4_stim8hr_cl'
    ]) AS dataset_id
  ),
  dea AS (
    SELECT dataset_id, COUNT(*) AS dea_rows
    FROM `prj-ext-dev-pertcat-437314.perturb_seq.pertpy_dea`
    WHERE dataset_id IN (SELECT dataset_id FROM wanted)
    GROUP BY dataset_id
  ),
  gsea AS (
    SELECT dataset_id, COUNT(*) AS gsea_rows
    FROM `prj-ext-dev-pertcat-437314.perturb_seq.pertpy_gsea`
    WHERE dataset_id IN (SELECT dataset_id FROM wanted)
    GROUP BY dataset_id
  )
SELECT
  w.dataset_id,
  m.study_title,
  m.study_year,
  m.first_author,
  m.data_modalities,
  m.perturbation_type_labels,
  m.tissue_labels,
  m.cell_type_labels,
  m.cell_line_labels,
  m.number_of_perturbed_targets,
  m.number_of_perturbed_samples,
  IFNULL(dea.dea_rows, 0) AS dea_rows,
  IFNULL(gsea.gsea_rows, 0) AS gsea_rows
FROM wanted AS w
LEFT JOIN `prj-ext-dev-pertcat-437314.unified_data_ensg.dataset_summary` AS m USING (dataset_id)
LEFT JOIN dea USING (dataset_id)
LEFT JOIN gsea USING (dataset_id)
ORDER BY w.dataset_id;

-- Representative metadata JSON fields for adamson_2016_pilot.
SELECT
  dataset_id,
  study_title,
  study_year,
  first_author,
  data_modalities,
  perturbation_type_labels,
  tissue_labels,
  cell_type_labels,
  cell_line_labels,
  number_of_perturbed_targets,
  number_of_perturbed_samples
FROM `prj-ext-dev-pertcat-437314.unified_data_ensg.dataset_summary`
WHERE dataset_id = 'adamson_2016_pilot';

-- Representative DEA export row (dataset_id is used as the file key).
SELECT
  perturbed_target_ensg,
  perturbed_target_symbol AS perturbed_target_name,
  effect_gene_ensg,
  effect_gene_symbol AS effect_gene_name,
  log2foldchange,
  padj,
  score_name,
  score_value,
  cell_type
FROM `prj-ext-dev-pertcat-437314.perturb_seq.pertpy_dea`
WHERE dataset_id = 'adamson_2016_pilot'
  AND perturbed_target_ensg = 'ENSG00000134107'
  AND effect_gene_ensg = 'ENSG00000094914'
  AND score_name = 'Wilcoxon Score'
LIMIT 1;

-- Representative GSEA export row; the array is shortened only for display.
SELECT
  perturbed_target_ensg,
  perturbed_target_symbol,
  term,
  es,
  nes,
  pval,
  sidak,
  fdr,
  geneset_size,
  ARRAY_SLICE(leading_edge, 0, 4) AS leading_edge_example,
  cell_type
FROM `prj-ext-dev-pertcat-437314.perturb_seq.pertpy_gsea`
WHERE dataset_id = 'adamson_2016_pilot'
  AND perturbed_target_ensg = 'ENSG00000066336'
  AND term = 'HALLMARK_ALLOGRAFT_REJECTION'
LIMIT 1;

-- Source types for the data fields selected by the release export.
SELECT table_name, column_name, data_type, is_nullable
FROM `prj-ext-dev-pertcat-437314.perturb_seq.INFORMATION_SCHEMA.COLUMNS`
WHERE table_name IN ('pertpy_dea', 'pertpy_gsea')
ORDER BY table_name, ordinal_position;

-- Representative CRISPR release row (score_value is exported as STRING).
SELECT
  dataset_id,
  perturbed_target_ensg,
  perturbed_target_symbol AS perturbed_target_name,
  score_name,
  score_value,
  significant,
  significance_criteria
FROM `prj-ext-dev-pertcat-437314.crispr.data`
WHERE dataset_id = 'biogrid_1161'
LIMIT 1;

-- Representative MAVE release row with the same derived fields as the export.
SELECT
  dataset_id,
  perturbed_target_ensg,
  perturbed_target_symbol AS perturbed_target_name,
  perturbation_name,
  SAFE_CAST(REGEXP_EXTRACT(perturbation_name, r'p\.[A-Za-z]+(\d+)') AS INT64)
    AS perturbation_position,
  REGEXP_EXTRACT(perturbation_name, r'p\.([A-Za-z]+)\d+')
    AS perturbation_aa_wt,
  REGEXP_EXTRACT(perturbation_name, r'p\.[A-Za-z]+\d+([A-Za-z=]+)')
    AS perturbation_aa_change,
  score_name,
  score_value
FROM `prj-ext-dev-pertcat-437314.mavedb.data`
WHERE dataset_id = 'urn:mavedb:00000001-a-1'
LIMIT 1;
