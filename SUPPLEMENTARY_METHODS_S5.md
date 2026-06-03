# Supplementary Methods S5 — Custom Pfam HMM panel (12 profiles)

A custom HMM panel comprising 12 Pfam profiles relevant to mobile-element machinery and cell-wall hydrolysis was scanned against the Prokka-predicted proteome of PL25 contig_11 with HMMER v3.4 (Eddy, 2011) using `hmmscan` at an E-value threshold of ≤ 10⁻⁵. Pfam profiles were obtained from the Pfam database release 37.0 (Mistry et al., 2021), downloaded 2026-05-15.

The 12 profiles, grouped by functional category, are:

| Category | Pfam ID | Profile name | Function |
|---|---|---|---|
| **Cell-wall hydrolase — M23** | PF01551 | Peptidase_M23 | Zn²⁺-dependent endopeptidase (cleaves D-Ala–Gly cross-bridges) |
| **Cell-wall hydrolase — SLT** | PF01464 | SLT | Soluble lytic transglycosylase (Slt35/Slt70-like glycan cleavage) |
| **Cell-wall hydrolase — GH25** | PF01183 | Glyco_hydro_25 | GH25 muramidase fold (CwlT-type) |
| **Cell-wall hydrolase — NlpC/P60** | PF00877 | NLPC_P60 | NlpC/P60 cysteine-peptidase fold |
| **Relaxase / replication-initiator** | PF03389 | MobA_MobL | Plasmid replication-relaxation protein family |
| **T4SS coupling protein** | PF02534 | T4SS_TraD-like | DNA-transfer coupling protein |
| **VirB4-like ATPase** | PF03135 | TrwB_AAD_bind | VirB4-family AAA-ATPase associated with T4SS energetics |
| **AAA-domain ATPases (generic)** | PF13245 | AAA_19 | RecD2-like helicase (matches PL25_00094) |
|  | PF13191 | AAA_16 | Generic AAA ATPase |
|  | PF12846 | AAA_10 | AAA-domain ATPase, VirB4-associated |
| **DNA-repair / branch-resolution** | PF05866 | RusA | Holliday-junction resolvase (mobile-element marker) |
| **Cell-wall-binding** | PF08239 | SH3_3 | SH3b cell-wall-targeting domain (M23-associated) |

Hits with E ≤ 10⁻⁵ were retained as supporting evidence for the functional annotation of PL25_00076 (replication-relaxation), PL25_00082 (VirB4-like ATPase), PL25_00094 (RecD2-like helicase), and PL25_00148 (T4SS-DNA transfer coupling protein), as reported in the main text. The full panel HMM file is provided at `panels/ICE_M23_panel.hmm` in the project Zenodo archive.

## References (S5)

- Eddy SR. Accelerated profile HMM searches. *PLoS Comput Biol*. 2011;7:e1002195.
- Mistry J, Chuguransky S, Williams L, et al. Pfam: The protein families database in 2021. *Nucleic Acids Res*. 2021;49:D412–D419.
