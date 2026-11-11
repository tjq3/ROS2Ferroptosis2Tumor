# External directional consistency

This analysis compares the directional CTMC prediction with independent FSP1 lipidomics and orthotopic tumour growth data. It uses the commitment and clearance scan output with the original parameters and thresholds.

1. Generate `Macro_Q1_TargetedRefinement_FINAL.mat` by running the [CTMC commitment and clearance scan](../CTMC/). The script looks for the MAT file in the CTMC output directory.
2. Download the following workbooks from Wu et al., *Targeting FSP1 triggers ferroptosis in lung cancer*, *Nature* 649, 487–495 (2026), and place them in this directory with their original filenames. [Article and data access](https://doi.org/10.1038/s41586-025-09710-8)
   - `41586_2025_9710_MOESM3_ESM.xlsx`. Use sheet `Tumor Lipidomics Raw`, cells `B4:L27`. [Supplementary Data 3](https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-025-09710-8/MediaObjects/41586_2025_9710_MOESM3_ESM.xlsx)
   - `41586_2025_9710_MOESM14_ESM.xlsx`. Use sheet `4f`, cells `B4:O7`. [Source Data Extended Data Fig. 4](https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-025-09710-8/MediaObjects/41586_2025_9710_MOESM14_ESM.xlsx)
3. Run `ExternalValidation_FSP1Consistency`. It writes the figure, the endpoint summary, the lipid-species table and the mouse-level area-under-the-tumour-curve table under `Fig6_ExternalConsistency_FSP1_LIGHT_FINAL_TOPJOURNAL_OUT/`.

The study copies `LPO.xlsx` and `tumor.xlsx` are also accepted, with the same cell ranges on their first sheets.
