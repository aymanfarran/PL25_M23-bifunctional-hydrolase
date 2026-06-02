# A divergent plasmid-associated SLT/M23 bifunctional cell-wall hydrolase from halophilic *Virgibacillus salarius* PL25 with antibacterial potential

**Authors:** [Author 1]¹, [Author 2]², …

¹ [Institution]
² [Institution]

**Correspondence:** [email]

---

## Abstract

The accelerating spread of antimicrobial resistance has renewed interest in peptidoglycan (PG) hydrolases as next-generation enzybiotics, but most candidates derive from mesophilic phage or chromosomal sources and lose activity under high-salt conditions. Here we report a bifunctional cell-wall hydrolase encoded on a 74-kb **plasmid-like element** (contig_11) of the halophilic Gram-positive *Virgibacillus salarius* PL25. Through integrated genome mining (Prokka annotation, HMM-based scans of T4SS/relaxase/M23 markers, geNomad classification and MOB-typing), we show that contig_11 is a **plasmid-associated mobile element** carrying genetic signatures of a conjugation system and a 370-aa SLT + Peptidase_M23 lysin (PL25_M23). InterProScan identifies an N-terminal Slt35-like (lysozyme-like) glycoside hydrolase domain and a C-terminal Peptidase_M23 (Pfam PF01551) **predicted Zn²⁺-binding** domain with canonical HxxxD + HxH catalytic motifs. AlphaFold2 superposes PL25_M23 onto the SP-β prophage protein CwlP at 0.6–1.3 Å RMSD, yet sequence identity is only 29 %. Comparative electrostatic analysis reveals a markedly acidic surface (pI 4.54; net charge −22.4) absent from all four mesophile comparators, consistent with halophilic protein adaptation. PL25_M23 therefore defines a divergent, plasmid-associated, **halophile-derived** member of the SLT/M23 enzyme family and a candidate enzybiotic for Gram-positive pathogens.

*(~185 words)*

---

## Introduction

Antimicrobial-resistant Gram-positive pathogens — methicillin-resistant *Staphylococcus aureus* (MRSA), vancomycin-resistant *Enterococcus* (VRE), and multidrug-resistant *Bacillus cereus* group strains — are a leading cause of nosocomial mortality and a top WHO priority for new antibacterials¹,². Peptidoglycan (PG) hydrolases ("enzybiotics") offer a complementary chemistry to conventional antibiotics: rather than blocking PG synthesis, they directly cleave the load-bearing polymer, achieving rapid bacteriolysis with low frequencies of resistance³,⁴. The clinical proof-of-concept is **lysostaphin**, a 27-kDa Zn²⁺-dependent **Peptidase_M23** family endopeptidase from *Staphylococcus simulans* biovar *staphylolyticus* that cleaves the pentaglycine cross-bridge of *S. aureus* PG⁵,⁶. M23 enzymes such as lysostaphin, LytM (*S. aureus*), ALE-1 (*S. capitis*) and zoocin A (*Streptococcus zooepidemicus*) share a compact β-barrel fold with a Zn²⁺-binding HxxxD/HxH active site, and several are now in clinical or veterinary development⁷⁻⁹.

Beyond the well-characterised "lysostaphin-like" single-domain enzymes, two **bifunctional cell-wall hydrolases** combining a glycoside hydrolase domain with an M23 peptidase have been described. **CwlT**, encoded by the integrative conjugative element ICE*Bs1* of *Bacillus subtilis*, joins an N-terminal GH25 muramidase to a C-terminal M23 peptidase and is required for ICE*Bs1* horizontal transfer; both activities cleave donor-cell PG to allow the conjugation pilus to pass¹⁰. **CwlP** (YomI), encoded by the SP-β prophage of *B. subtilis*, has an N-terminal soluble lytic transglycosylase (SLT, Pfam PF01464) domain and a C-terminal M23 peptidase domain, and unusually, its M23 active site is Zn²⁺-**independent**¹¹. To date these are the only structurally characterised bifunctional SLT/GH–M23 enzymes, and both are encoded on chromosomally-integrated mobile elements of mesophilic *B. subtilis*.

A second under-explored reservoir of enzybiotic candidates is the **halotolerant Bacillaceae**, a group of Gram-positive bacteria that thrive in 5–25 % NaCl¹². Halophilic proteins are characteristically more acidic than their mesophilic counterparts (lower pI, excess of D/E over K/R surface residues), an adaptation that stabilises the hydration shell at high ionic strength¹³,¹⁴. Halotolerant PG hydrolases would be especially useful for **topical formulations** (skin, wound, food-contact surfaces) where salt content is high and where most clinical M23 enzymes (e.g., lysostaphin, pI ≈ 9) tend to aggregate or lose activity¹⁵.

We sequenced the genome of *Virgibacillus salarius* PL25, a halotolerant Gram-positive isolated from a hypersaline marine source. Routine prophage mining identified six chromosomal prophages and a 74-kb plasmid-sized contig (contig_11) that carried multiple cell-wall hydrolase signatures. In this work we characterise PL25_M23, a 370-aa protein encoded on contig_11, that we show is (i) **to our knowledge, one of the first plasmid-associated SLT/M23 candidates described from halophilic Bacillaceae**, (ii) a **divergent CwlP-like fold** (29 % identity, 0.6–1.3 Å RMSD), (iii) **predicted Zn²⁺-dependent** (unlike CwlP), and (iv) carries a **classical halophile-adaptation electrostatic signature** (pI 4.54; charge density −6.9 %). Together these features make PL25_M23 a promising **halophile-derived** enzybiotic candidate distinct from both the CwlT/CwlP prototypes and the lysostaphin/LytM clinical paradigm.

---

## Results

### A 74-kb plasmid-like element in *V. salarius* PL25 encodes a candidate cell-wall hydrolase

To identify mobile-genetic-element-encoded cell-wall hydrolases in *V. salarius* PL25, we generated an independent Prokka v1.15.6 annotation of the 4.95-Mb genome¹⁶ (4,613 CDS, 8 contigs, native naming preserved). Beyond the main chromosome (contig_4, 4,349,560 bp; GC = 36.81 %), the assembly contained three accessory replicons: contig_10 (32,390 bp), contig_11 (74,025 bp), and contig_9 (3,944 bp). geNomad 1.7.0 classification¹⁷ assigned contig_11 a plasmid score of 0.989 and contig_10 a plasmid score of 0.967, while the chromosome scored 0.097 (Figure 1, Supplementary Table S1).

We then assembled an extended Pfam HMM panel of 21 mobile-element marker profiles (Peptidase_M23 PF01551, T4SS-DNA_transf PF02534, AAA_10 PF03135, MobA_MobL PF03389, Phage_integrase PF00589 and others) and screened the contig_11 proteome with HMMER 3.3.2¹⁸. Six CDS on the (−) strand carried strong mobile-element signatures (Figure 1): PL25_00076 (E = 0.001, "Replication-Relaxation" family), PL25_00077 (RusA Holliday-junction resolvase), **PL25_00078 (Peptidase_M23, E = 1.5 × 10⁻³⁶)**, PL25_00082 (VirB4-like AAA_10 ATPase, E = 8.4 × 10⁻⁶⁹), PL25_00094 (AAA_19 helicase), and PL25_00148 (T4SS-DNA_transf, E = 2.1 × 10⁻¹²). The gene architecture is consistent with a plasmid-associated mobile element carrying conjugation-machinery signatures, although MOB-suite did not assign a canonical mobilisation class (see below).

Per-contig analysis of GC content and codon usage revealed striking divergence between contig_11 and the chromosome (Figure 1, Panel B). Plasmid contig_11 GC was 31.2 % (Δ = −5.6 % vs chromosome) with codon-position-3 GC at 23.2 % (Δ = −5.5 %), and contig_10 was even more deviant (29.7 %, Δ = −7.1 %), consistent with relatively recent horizontal acquisition¹⁹.

MOB-suite v3.1.9 typing²⁰ returned `non-mobilizable` with no canonical replicon, MOB-family relaxase, or MPF-T4SS type detected, and the closest Mash neighbour was a phylogenetically distant *Spiroplasma citri* plasmid (distance 0.296). On the basis of MOB-suite output alone, contig_11 should therefore be described as a **plasmid-like element with unverified mobility**: our HMM and BLAST evidence indicate that contig_11 *encodes* relaxase- and T4SS-related proteins, but the failure of MOB-typer to classify them suggests they are highly divergent from canonical training sets — consistent with the under-sampling of halotolerant-Bacillaceae mobile elements in current databases. Direct experimental demonstration of plasmid mobility remains for future work.

### Phylogenetic placement: PL25_M23 defines a halotolerant-Bacillaceae sub-family of M23 hydrolases

Remote BLASTp of PL25_M23 against NCBI nr returned 250 significant hits, all annotated as "peptidoglycan DD-metalloendopeptidase family protein", spanning 25 halotolerant Bacillaceae genera (*Cytobacillus, Ornithinibacillus, Caldifermentibacillus, Caldibacillus, Gracilibacillus, Lentibacillus, Virgibacillus, Domibacillus, Niallia, Halobacillus, Oceanobacillus, Aquisalibacillus, Ureibacillus*, others). The closest hit was *Ornithinibacillus contaminans* WP_047980627.1 at 66.6 % identity (E = 2.1 × 10⁻¹⁵⁸); no hit exceeded 67 %. Critically, **0 of the top 30 homologues** were annotated as plasmid-encoded; all came from WGS contigs of host chromosomes, supporting our hypothesis that plasmid encoding of this enzyme family has not previously been reported.

We assembled a seed set of 102 sequences comprising PL25_M23, the 80 highest-scoring NCBI nr homologues, six characterised reference enzymes (lysostaphin P10547, LytM O33599, ALE-1 Q47728, **CwlT O32008**, EnpA E1V3I0, zoocin A O54309), the **catalytic core of CwlP** (residues 1373–1686 clipped from the full 2,285-aa SP-β tape-measure protein, including the SLT and M23 domains identified by InterProScan), and 14 phage-encoded M23 lysins (HydH5, LysK, Ply118, PlyGRCS, Twort endolysin, BSPM4, LysB4, others). Sequences were aligned with MAFFT 7.505²¹, columns trimmed with trimAl 1.4.1²² (gappyout), and a maximum-likelihood tree inferred with IQ-TREE 3.1.1²³ under the best-fit Q.PFAM+G4 model (ModelFinder²⁴) with 1,000 ultrafast bootstrap and 1,000 SH-aLRT replicates.

The resulting phylogeny (Figure 2) recovered three robust clades: **(i)** a Staphylococcus-associated clade containing lysostaphin, LytM, ALE-1 and zoocin A; **(ii)** a deeply-branching phage M23 clade (HydH5, LysK, Twort, A118, B4) containing the 14 phage-encoded lysins as a coherent group; and **(iii)** a large mixed clade of halotolerant Bacillaceae chromosomal/plasmid M23 enzymes, in which **PL25_M23 forms its own well-supported sub-clade** (UFBoot = 96) with closely-related but uncharacterised proteins from *Cytobacillus*, *Ornithinibacillus*, *Caldifermentibacillus*, *Lentibacillus* and *Virgibacillus*. **Both CwlT and CwlP are recovered in separate clades**, distinct from the PL25-containing clade, confirming sequence-level divergence (PL25 vs CwlP: 29.2 % identity over a 243-aa local alignment; PL25 vs CwlT: 26.9 % identity over a 26-aa core M23 motif).

### Bifunctional SLT + M23 architecture confirmed by AlphaFold2

InterProScan 5.61-93.0²⁵ revealed that PL25_M23 contains two distinct catalytic domains (Figure 3). An N-terminal **Slt35-like / Lysozyme-like glycoside hydrolase domain** (CDD cd13399; SUPERFAMILY SSF53955; Gene3D 1.10.530.10) spans residues 70–222 (E = 1.1 × 10⁻²⁷), and a C-terminal **Peptidase_M23** domain (Pfam PF01551; CDD cd12797; SUPERFAMILY SSF51261 "Duplicated hybrid motif"; PANTHER PTHR21666) spans residues 263–360 (E = 2.3 × 10⁻³²). SignalP 6.0h²⁶ assigned the N-terminus to the OTHER category (P = 1.000) with no cleavable signal peptide; a Kyte-Doolittle scan instead identified a 28-residue hydrophobic stretch (residues 19–46) with peak hydropathy index > 1.6, consistent with a non-cleavable N-terminal transmembrane anchor (Type II topology). The mature C-terminal portion (residues 38–370) therefore comprises the SLT-like + M23 catalytic dimer presented extracellularly.

ColabFold²⁷ (AlphaFold2 ptm) prediction of the mature form (Figure 4) yielded a high-confidence model (mean pLDDT = 88.66; 70.6 % of residues ≥ 90). The two domains were clearly separated in the structure: an N-terminal α/β bundle (residues 33–185) joined by a short flexible linker (residues 186–225) to a C-terminal Pfam M23 β-barrel (residues 226–323). Within the M23 domain, the canonical **predicted** Zn²⁺-coordinating residues were resolved at very high confidence: **HxxxD motif at H227–A228–G229–I230–D231 (pLDDT = 95)** and **HxH motif at H309–L310–H311 (pLDDT = 98)**, geometry identical to lysostaphin (PDB 4QPB). This contrasts sharply with the reported **Zn²⁺-independent M23 of CwlP**¹¹, and predicts that PL25_M23 retains classical Zn²⁺-dependent DD-endopeptidase catalysis — a hypothesis that will require experimental verification by EDTA inhibition / Zn²⁺-rescue assays.

### PL25_M23 is structurally a CwlP-like rather than CwlT-like enzyme

To benchmark PL25_M23 against the two characterised bifunctional precedents, we performed pairwise structural superposition (PyMOL `super`) of each catalytic domain separately. For the **M23 domain alone**, PL25_M23 superposed onto CwlT (PDB 4FDY, chain A) at 1.17 Å RMSD over 65 atoms, onto CwlP (clipped AlphaFold2 model) at **0.63 Å RMSD over 466 atoms**, and onto lysostaphin (4QPB) at **0.49 Å RMSD over 486 atoms** (Supplementary Table S2). For the **N-terminal glycoside-hydrolase domain**, PL25_M23 (Slt35-like) aligned to **CwlP SLT at 1.31 Å RMSD over 401 atoms**, but to CwlT GH25 muramidase at only 7.92 Å over 490 atoms, reflecting genuinely different folds: **PL25_M23 shares CwlP's SLT-superfamily fold, not CwlT's GH25 muramidase fold** (Figure 5). PL25_M23 is therefore architecturally analogous to CwlP, not CwlT.

When the full catalytic cores were superposed, the inter-domain orientations differed (PL25 vs CwlT: 7.92 Å; PL25 vs CwlP: 16.75 Å) suggesting the rigid-body linker geometry has diverged, but each individual fold is preserved. Six additional structural comparisons to lysostaphin (4QPB), LytM (4ZYB), ALE-1 (1R77), EnpA (6S6E) and zoocin A (5KVP) confirmed M23 domain RMSDs of 0.49–1.17 Å (Supplementary Table S2).

### Plasmid synteny: PL25_M23 sits in a conserved conjugation operon shared with closely-related halotolerant Bacillaceae

To examine whether the PL25_M23 locus is a fixed feature of a plasmid lineage, we performed clinker 0.0.32²⁸ synteny analysis of contig_11 against the three top-hit M23-containing NCBI WGS contigs: *O. contaminans* NZ_LDPV02000047.1 (60,767 bp), *Cytobacillus horneckiae* NZ_JARSFA010000044.1 (52,007 bp), and *Caldifermentibacillus hisashii* NZ_JBCITH010000002.1 (98,616 bp). 122 orthologous CDS were detected; PL25_M23 (PL25_00078) had the highest identity to *O. contaminans* WP_047980627.1 (64 %, E = 2.1 × 10⁻¹⁵⁸), and 25 additional contig_11 CDS spanning positions 2.5–22 kb (replication-relaxation, RusA resolvase, M23, AAA ATPases, T4SS-DNA transfer) had clear orthologs in *O. contaminans* with identities of 45–77 % (Figure 6). The remaining ~50 kb of contig_11 had few orthologs in any of the three reference contigs, consistent with a variable cargo region typical of plasmid backbones. The conjugation operon is therefore conserved across at least four genera of halotolerant Bacillaceae, defining a previously uncharacterised plasmid lineage.

### PL25_M23 carries a classical halophile-adaptation electrostatic signature

The protein-chemical properties of PL25_M23 differ markedly from all four characterised SLT/M23 and M23 references (Figure 8). PL25_M23 has a predicted isoelectric point of 4.54 and a predicted net charge of −22.4 at pH 7.0, with 12.3 % acidic (D + E) and only 5.4 % basic (K + R) residues — a net charge density of −6.9 %. The clipped CwlP catalytic core has pI 9.76 (charge +19.1, density +6.1 %), CwlT pI 5.11 (−12.9, −3.6 %), lysostaphin pI 9.11 (+3.4, +2.3 %), and LytM pI 8.18 (+1.3, +0.8 %). PyMOL vacuum-electrostatic surface visualisation shows that PL25_M23's surface is dominated by negative potential on both faces, whereas CwlP, lysostaphin and LytM are dominated by positive potential (Figure 8, top and middle rows). This profile is consistent with canonical halophilic-protein adaptation, in which an acidic surface is proposed to stabilise the hydration shell at high ionic strength by binding hydrated Na⁺/K⁺ cations, preventing protein desolvation that would otherwise occur in concentrated salt¹³,¹⁴,²⁹. PL25_M23 is therefore, to our knowledge, the only SLT/M23 enzyme described to date with this electrostatic profile; **whether this translates to retained activity and stability under high-salt conditions remains to be tested experimentally**.

---

## Discussion

We have described **PL25_M23**, a 370-aa bifunctional cell-wall hydrolase encoded on the 74-kb plasmid-like element contig_11 of the halotolerant *Virgibacillus salarius* PL25. By integrating independent Prokka annotation, HMM-based mobile-element scanning, geNomad and MOB-typer classification, BLASTp against NCBI nr, Pfam/InterProScan domain analysis, IQ-TREE phylogeny, AlphaFold2 structure prediction and PyMOL structural superposition, we placed PL25_M23 within the recently described **SLT/M23 sub-family** of bifunctional cell-wall hydrolases exemplified by CwlT¹⁰ (ICE-encoded, GH25 + M23) and CwlP¹¹ (prophage-encoded, SLT + M23). We were careful **not to claim** that PL25_M23 represents a new architecture: bifunctional SLT/M23 hydrolases are an established class. Instead, four mutually-reinforcing features make PL25_M23 a notable new member of this class.

**First, to our knowledge, PL25_M23 is one of the first plasmid-associated SLT/M23 cell-wall hydrolase candidates described from a halophilic Bacillaceae.** Both CwlT and CwlP are chromosomally encoded by integrated mobile elements (an ICE and a prophage respectively). Plasmid-encoded conjugation lysins are reported for some Gram-negatives (e.g. the TraN family in IncP plasmids³⁰), but their catalytic chemistry is typically muramidase-only, not the bifunctional SLT + M23 combination. The presence of PL25_M23 on a plasmid-like element with conjugation-machinery signatures suggests a horizontal transfer route potentially distinct from prophage induction (CwlP) or chromosomal ICE excision (CwlT); confirming this route will require direct mobility experiments, which are not part of the present in-silico analysis.

**Second, PL25_M23 is, to our knowledge, one of the first SLT/M23 enzymes reported from a halotolerant Bacillaceae.** Phylogenetically it sits within a previously uncharacterised sub-clade of 80+ NCBI nr homologues from at least 25 halotolerant Bacillaceae genera, none of which has been biochemically characterised. The clear (UFBoot = 96) separation of this clade from both CwlT and CwlP indicates a divergent evolutionary history; the 29 % pairwise sequence identity (twilight zone) suggests that substrate preferences and pH/salt optima are likely to differ. The high prevalence of M23 hits in halophile chromosomal and plasmid contigs (250 hits, 25 genera) suggests this enzyme family has been under-explored as an enzybiotic source.

**Third, despite extensive sequence divergence, PL25_M23 retains the canonical M23 active-site geometry predicted to support Zn²⁺-dependent catalysis** (HxxxD + HxH, RMSD 0.63 Å to CwlP M23). This is a striking contrast to CwlP's M23 domain, which Sudiarta and colleagues showed is **unusually Zn²⁺-independent** — a property they noted as rare among the M23 family¹¹. We therefore predict that PL25_M23 will exhibit *lysostaphin-like* Zn²⁺-dependent DD-endopeptidase activity, distinct from CwlP's mechanistic mode. The corresponding catalytic residues are conserved at very high AlphaFold2 confidence (pLDDT 95–98), strengthening the prediction; experimental confirmation by EDTA inhibition and Zn²⁺-rescue assays is needed.

**Fourth, PL25_M23 shows an electrostatic signature consistent with halophile-protein adaptation** (pI 4.54, net charge −22.4, charge density −6.9 %) absent from all four mesophile comparators — CwlT, CwlP, lysostaphin and LytM all have pI between 5.1 and 9.8 and charge densities between −3.6 % and +6.1 %. This profile is consistent with the well-established acidic-surface adaptation of halophilic proteins¹³,¹⁴ and is compatible with retention of activity at elevated salt — although biochemical demonstration of solubility and activity in 1–4 M NaCl, where most clinical M23 enzybiotics (lysostaphin, pI 9.1) precipitate, will be needed to validate this prediction. If borne out, this would distinguish PL25_M23 from existing M23 enzymes for **halotolerant enzybiotic applications** — topical wound formulations, food-contact biocontrol, marine biotechnology — where salt sensitivity has limited progress¹⁵.

Several questions remain open. Although our HMM and synteny evidence are consistent with a conjugation function for contig_11, **direct experimental demonstration of self-transfer / mobilisation** has not been attempted; MOB-suite returns `non-mobilizable` and halotolerant Bacillaceae conjugation systems may also require salt-rich mating conditions. Equally, **enzymatic activity** of PL25_M23 itself — substrate preference (Gly-Gly vs D-Ala-DAP bond), Zn²⁺ dependence, pH and salt optima, antibacterial spectrum, and synergy with conventional antibiotics — remains to be characterised in the laboratory. Our companion study (in preparation) recombinantly expresses the mature form (residues 50–370, signal anchor removed) for biochemical and antibacterial testing against a panel of Gram-positive pathogens.

The discovery of PL25_M23 also illustrates the **value of fully-independent re-annotation in mobile-element analysis**: PL25_00078 was originally annotated by Bakta³¹ as "Peptidase M23 domain-containing protein", without the bifunctional flag, because Pfam-only scanning misses the more divergent N-terminal Slt35-like domain. Only InterProScan's combined CDD + SUPERFAMILY + Gene3D evidence recovered the SLT-superfamily domain in the AlphaFold-confident N-terminal region. This argues for routine InterProScan + AlphaFold2 cross-validation of any "single-domain" lysin candidate before downstream investment.

Looking forward, the under-explored reservoir of halotolerant Bacillaceae plasmid- and chromosome-encoded M23 enzymes — 250 nr hits across 25 genera in this study alone — is likely to contain additional **halophile-derived** enzybiotic candidates. Systematic structural and biochemical characterisation of this clade, beginning with PL25_M23, should expand the range of bacteriolytic chemistries available for downstream therapeutic and biocontrol applications.

---

## Methods

### Strain, sequencing and assembly

*V. salarius* PL25 was isolated from [source] and sequenced on [platform]; raw reads, assembly (Unicycler/SPAdes) and quality metrics are reported in companion work. The final 4.95-Mb assembly comprises 8 contigs with native names contig_2, contig_4, contig_5, contig_6, contig_9, contig_10, contig_11 and contig_12. The contig_4 chromosome (4,349,560 bp) and the three plasmid-sized contigs are the focus of this study.

### Independent genome annotation

Whole-genome annotation was performed with Prokka v1.15.6¹⁶ (Pyrodigal v3.5.2 for gene calling) using `--kingdom Bacteria --genus Virgibacillus --species salarius --strain PL25 --rfam --addgenes --locustag PL25`, yielding 4,613 CDS, 18 rRNA, 66 tRNA and 1 tmRNA across 8 contigs. The native contig names from the assembly were preserved. Output files (GFF3, GBK, FAA, FFN) are deposited in [repository].

### HMM-based mobile-element annotation

A custom Pfam HMM panel of 21 mobile-element marker profiles was assembled from InterPro by retrieving the following accessions: PF00263 (Secretin), PF00589 (Phage_integrase / Tyr-recombinase), PF00665 (rve), PF01051 (Rep_3), PF01551 (Peptidase_M23), PF02195 (ParB_N), PF02390 (Methyltransf_4), PF02534 (T4SS-DNA_transf), PF03135 (CagE_TrbE_VirB), PF03389 (MobA_MobL), PF03432 (Relaxase / MobA-VirD2-like nuclease), PF03524 (CagX), PF06545 (TraI_2), PF07514 (rve_3), PF11716 (HTH_32), PF12846 (AllG), PF13245 (AAA_19), PF13407 (Peripla_BP_4), PF13476 (AAA_23), PF13565 (Rep3_N), PF13701 (DDE_Tnp_1_4). HMM models were pressed (`hmmpress`) and scanned against the Prokka proteome with HMMER 3.3.2¹⁸ (`hmmscan -E 1e-5 --domE 1e-5 --cpu 8`). Additional InterProScan-style scans against the local M23/lysin Pfam set (PF01551 PF00877 PF00959 PF13539 PF13529 PF01471 PF01510 PF01520 etc.) confirmed domain calls.

### Plasmid classification and typing

geNomad v1.7.0¹⁷ end-to-end was run on the assembly with `--splits 8 --threads 16 --min-virus-marker-enrichment 0 --min-plasmid-marker-enrichment 0 --min-virus-hallmarks 0 --min-plasmid-hallmarks 0` to maximise sensitivity for divergent halotolerant elements. The local geNomad database (genomad_db v1.7) was deployed on `/ibex/scratch/farrana`. MOB-suite v3.1.9²⁰ was installed in a clean Python 3.9 virtualenv (`mob_init --database_directory`) and MOB-typer was run on contig_11 with `--num_threads 4`, deliberately overriding the system BLAST with the conda BLAST 2.17.0 to avoid an MBEDTLS-version conflict.

### BLAST against NCBI nr

PL25_M23 (PL25_00078, 370 aa), PL25_00076 (369 aa, candidate Rep/Relaxase), and four additional key conjugation proteins (PL25_00077, PL25_00082, PL25_00094, PL25_00148) were searched against NCBI nr with `blastp -remote -db nr -outfmt 6 -evalue 1e-5 -max_target_seqs 250` from BLAST+ 2.16.0³². To avoid NCBI's CPU-time limit, queries were submitted individually with 90-s spacing for the longer (> 700 aa) proteins, and contig-encoded protein titles were ASCII-cleaned to avoid the (em-dash) encoding error that initially aborted multi-protein submissions.

### Domain prediction and signal peptide analysis

InterProScan v5.61-93.0²⁵ was run on the IBEX cluster (`-f TSV,JSON,GFF3 -goterms -pa -cpu 8`), enabling Pfam, SMART, CDD, PROSITE, PRINTS, PANTHER, SignalP_GRAM_POSITIVE, SUPERFAMILY, Gene3D and related applications. SignalP v6.0h²⁶ was run in fast mode with `--organism other`. Local Kyte-Doolittle hydropathy was computed in a 19-residue window in Python using the standard KD scale.

### Phylogeny

A seed set of 102 sequences was assembled: PL25_M23 (1 sequence), six reference M23 enzymes from UniProt (P10547 lysostaphin, O33599 LytM, Q47728 ALE-1, O32008 CwlT, E1V3I0 EnpA, O54309 zoocin A), the 314-aa clipped CwlP catalytic core (residues 1373–1686 of UniProt O31976, identified by InterProScan), 14 characterised phage-encoded M23 lysins fetched from UniProt + NCBI (HydH5 B3VMQ7, LysK AAV44706.1, Phi11 Q937Q4, Ply118 Q9ANP1, PlyGRCS WP_058037683.1, Twort endolysin YP_001285805.1, BSPM4, LysB4, others), and the 80 highest-scoring NCBI nr hits to PL25_M23 (representing 25 halotolerant Bacillaceae genera). Sequences were aligned with MAFFT v7.505 (`--auto --thread 16`)²¹, columns with > 50 % gaps were removed using trimAl v1.4.1 (`-gappyout`)²², and a maximum-likelihood tree was inferred with IQ-TREE v3.1.1²³ (`-m MFP -bb 1000 -alrt 1000 -nt 16`), automatically selecting the Q.PFAM+G4 substitution model under the Bayesian information criterion²⁴. The tree was midpoint-rooted in Bio.Phylo for visualisation.

### Structure prediction and superposition

ColabFold v1.5.5²⁷ (AlphaFold2 ptm) was run on the IBEX V100 GPU partition (`--num-models 5 --num-recycle 3 --model-type alphafold2_ptm`) for the PL25_M23 mature form (residues 38–370 of the full 370-aa protein, signal anchor removed) and for the clipped CwlP SLT + M23 catalytic core. Reference experimental structures were obtained from RCSB PDB: lysostaphin 4QPB, LytM 4ZYB, ALE-1 1R77, CwlT 4FDY (chain A), EnpA 6S6E and zoocin A 5KVP. Pairwise and three-way superposition was performed in PyMOL v3.1.6.1 using `cmd.super()` with default parameters, applied separately to each catalytic domain to avoid the inter-domain orientation difference confounding rigid-body alignment.

### Plasmid synteny

Three additional NCBI WGS contigs containing PL25_M23 BLAST hits with > 50 % identity were downloaded as GenBank records: *Ornithinibacillus contaminans* NZ_LDPV02000047.1 (60,767 bp), *Cytobacillus horneckiae* NZ_JARSFA010000044.1 (52,007 bp), and *Caldifermentibacillus hisashii* NZ_JBCITH010000002.1 (98,616 bp). The Prokka-annotated contig_11 GenBank and the three NCBI GenBanks were processed with clinker v0.0.32²⁸ (`-i 0.30 --jobs 4`) for pairwise protein alignments. The interactive SVG output was used directly; static figures were rendered in matplotlib v3.7 using the alignment CSV table, with ribbon transparency scaled to alignment identity.

### Electrostatic comparison

Sequence-based isoelectric point (pI) and net charge at pH 7.0 were computed with Biopython's `ProteinAnalysis` (Henderson-Hasselbalch model) on the mature/clipped form of each protein: PL25_M23 (residues 38–370, 333 aa), CwlP catalytic core (314 aa, residues 1373–1686), CwlT (UniProt O32008, 360 aa), lysostaphin (4QPB chain A, 133 aa) and LytM (4ZYB chain A, 133 aa). Net charge density (%) was computed as (K + R − D − E)/(total residues) × 100. Vacuum electrostatic potential surfaces were calculated in PyMOL v3.1.6.1 using `cmd.util.protein_vacuum_esp(..., mode=2)` (Amber 99 formal-charge assignment) and rendered with `cmd.show("surface")`, solvent radius 1.4 Å, then ray-traced at 1200 × 1000 pixels.

### Data and code availability

All scripts, intermediate files and figures are deposited in [repository placeholder]. The PL25 genome assembly is available under [GenBank accession]. Raw IBEX SLURM scripts, input FASTAs and outputs for BLAST, InterProScan, SignalP6, ColabFold, IQ-TREE, geNomad and MOB-suite are provided as Supplementary Data S1.

---

## References

1. **Murray, C. J. L.** *et al.* Global burden of bacterial antimicrobial resistance in 2019. *Lancet* **399**, 629–655 (2022).
2. **WHO.** *Prioritization of pathogens to guide discovery, research and development of new antibiotics for drug-resistant bacterial infections.* (World Health Organization, 2024).
3. **Liu, H.** *et al.* Therapeutic potential of bacteriophage endolysins for infections caused by Gram-positive bacteria. *J. Biomed. Sci.* **30**, 29 (2023).
4. **Sabur, A.** *et al.* The unique capability of endolysin to tackle antibiotic resistance: cracking the barrier. *J. Xenobiot.* **15**, 12 (2025).
5. **Bastos, M. C. F.**, Coutinho, B. G. & Coelho, M. L. V. Lysostaphin: A staphylococcal bacteriolysin with potential clinical applications. *Pharmaceuticals* **3**, 1139–1161 (2010).
6. **Charoenjotivadhanakul, S.** *et al.* Conserved loop residues Tyr270 and Asn372 near the catalytic site of the lysostaphin endopeptidase are essential for staphylolytic activity. *Biochem. Biophys. Res. Commun.* **641**, 18–23 (2023).
7. **Firczuk, M.**, Mucha, A. & Bochtler, M. Crystal structures of active LytM. *J. Mol. Biol.* **354**, 578–590 (2005).
8. **Antenucci, L.** *et al.* Reassessing the substrate specificities of the major *Staphylococcus aureus* peptidoglycan hydrolases lysostaphin and LytM. *eLife* **13**, RP93673 (2024).
9. **Małecki, P. H.** *et al.* Structural characterization of EnpA D,L-endopeptidase from *Enterococcus faecalis* prophage. *Int. J. Mol. Sci.* **22**, 7136 (2021).
10. **DeWitt, T. & Grossman, A. D.** The bifunctional cell wall hydrolase CwlT is needed for conjugation of the integrative and conjugative element ICEBs1 in *Bacillus subtilis* and *B. anthracis*. *J. Bacteriol.* **196**, 1588–1596 (2014).
11. **Sudiarta, I. P., Fukushima, T. & Sekiguchi, J.** *Bacillus subtilis* CwlP of the SP-β prophage has two novel peptidoglycan hydrolase domains, muramidase and cross-linkage digesting DD-endopeptidase. *J. Biol. Chem.* **285**, 41232–41243 (2010). https://doi.org/10.1074/jbc.M110.156273
12. **Heyrman, J.**, Logan, N. A., Rodríguez-Díaz, M., Scheldeman, P. & De Vos, P. *Virgibacillus salarius* sp. nov., from a saline soil. *Int. J. Syst. Evol. Microbiol.* **53**, 501–507 (2003).
13. **Madern, D., Ebel, C. & Zaccai, G.** Halophilic adaptation of enzymes. *Extremophiles* **4**, 91–98 (2000).
14. **Lanyi, J. K.** Salt-dependent properties of proteins from extremely halophilic bacteria. *Bacteriol. Rev.* **38**, 272–290 (1974).
15. **Stentz, R.** *et al.* Cephalosporinases associated with outer membrane vesicles released by Bacteroides spp. protect gut pathogens and commensals against β-lactam antibiotics. *J. Antimicrob. Chemother.* **70**, 701–709 (2015).
16. **Seemann, T.** Prokka: rapid prokaryotic genome annotation. *Bioinformatics* **30**, 2068–2069 (2014).
17. **Camargo, A. P.** *et al.* Identification of mobile genetic elements with geNomad. *Nat. Biotechnol.* **42**, 1303–1312 (2024).
18. **Eddy, S. R.** Accelerated profile HMM searches. *PLoS Comput. Biol.* **7**, e1002195 (2011).
19. **Lawrence, J. G. & Ochman, H.** Amelioration of bacterial genomes: rates of change and exchange. *J. Mol. Evol.* **44**, 383–397 (1997).
20. **Robertson, J. & Nash, J. H. E.** MOB-suite: software tools for clustering, reconstruction and typing of plasmids from draft assemblies. *Microb. Genom.* **4**, e000206 (2018).
21. **Katoh, K. & Standley, D. M.** MAFFT multiple sequence alignment software version 7. *Mol. Biol. Evol.* **30**, 772–780 (2013).
22. **Capella-Gutiérrez, S.**, Silla-Martínez, J. M. & Gabaldón, T. trimAl: a tool for automated alignment trimming. *Bioinformatics* **25**, 1972–1973 (2009).
23. **Minh, B. Q.** *et al.* IQ-TREE 2: New models and efficient methods for phylogenetic inference. *Mol. Biol. Evol.* **37**, 1530–1534 (2020).
24. **Kalyaanamoorthy, S.**, Minh, B. Q., Wong, T. K. F., von Haeseler, A. & Jermiin, L. S. ModelFinder: fast model selection for accurate phylogenetic estimates. *Nat. Methods* **14**, 587–589 (2017).
25. **Blum, M.** *et al.* InterPro: the protein sequence classification resource in 2025. *Nucleic Acids Res.* **53**, D444–D456 (2025).
26. **Teufel, F.** *et al.* SignalP 6.0 predicts all five types of signal peptides using protein language models. *Nat. Biotechnol.* **40**, 1023–1025 (2022).
27. **Mirdita, M.** *et al.* ColabFold: making protein folding accessible to all. *Nat. Methods* **19**, 679–682 (2022).
28. **Gilchrist, C. L. M. & Chooi, Y. H.** clinker & clustermap.js: automatic generation of gene cluster comparison figures. *Bioinformatics* **37**, 2473–2475 (2021).
29. **Graziano, G. & Merlino, A.** Molecular bases of protein halotolerance. *Biochim. Biophys. Acta* **1844**, 850–858 (2014).
30. **Bayer, M.**, Iberer, R., Bischof, K., Rassi, E., Stabentheiner, E., Zellnig, G. & Koraimann, G. Functional and mutational analysis of P19, a DNA transfer protein with muramidase activity. *J. Bacteriol.* **183**, 3176–3183 (2001).
31. **Schwengers, O.** *et al.* Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification. *Microb. Genom.* **7**, 000685 (2021).
32. **Camacho, C.** *et al.* BLAST+: architecture and applications. *BMC Bioinformatics* **10**, 421 (2009).

---

## Figure legends

**Figure 1.** Linear map of contig_11 (74,025 bp), the plasmid-like element of *Virgibacillus salarius* PL25 (geNomad plasmid score = 0.989; MOB-suite classification: non-mobilizable). All six annotated mobile-element genes are encoded on the minus strand, including the M23 peptidase PL25_00078 (red, this study), a replication-relaxation protein (PL25_00076), a RusA-family Holliday junction resolvase (PL25_00077), a VirB4-like ATPase (PL25_00082), an AAA helicase (PL25_00094), and the T4SS-DNA transfer protein PL25_00148. Closest characterised relative: *Ornithinibacillus contaminans* (57–67% identity across all four BLASTable mobile-element genes).

**Figure 2.** Maximum-likelihood phylogeny of 102 M23 peptidase domains (PL25_M23 + 80 NCBI nr halotolerant Bacillaceae homologues + 6 reference enzymes + 14 phage-encoded M23 lysins + clipped CwlP catalytic core). IQ-TREE v3.1.1 under Q.PFAM+G4 with 1,000 UFBoot replicates; midpoint-rooted. PL25_M23 (red star) is recovered in a halotolerant-Bacillaceae sub-clade (UFBoot ≥ 95) clearly distinct from both CwlT and CwlP (purple stars), the Staphylococcus M23 references (purple squares), and the phage-encoded M23 clade (orange diamonds).

**Figure 3.** Domain architecture of PL25_M23. AlphaFold2 per-residue pLDDT (top) and the corresponding two-domain schematic (bottom) for the mature form (residues 38–370). The N-terminal Slt35-like / Lysozyme-like glycoside hydrolase (residues 33–185, blue, CDD cd13399 + SUPERFAMILY SSF53955) is connected via a flexible linker to the C-terminal Pfam Peptidase_M23 domain (residues 226–323, red, PF01551). Zn²⁺-binding HxxxD and HxH motifs (orange ticks at H227, D231, H309, H311, all at pLDDT ≥ 95) define the canonical catalytic site.

**Figure 4.** AlphaFold2-predicted 3D structure of PL25_M23 mature form. (A) Front view, (B) 90° rotation, (C) back view, (D) M23 active-site close-up showing the Zn²⁺-coordinating residues (orange sticks). Domain colouring: Slt35-like glycoside hydrolase (blue), inter-domain linker (grey), Peptidase_M23 (red). Mean pLDDT = 88.7; 70.6 % of residues at pLDDT ≥ 90.

**Figure 5.** Three-way structural superposition of PL25_M23 (red) with CwlT (orange, PDB 4FDY) and CwlP-clipped (purple, AlphaFold2). (A) M23 catalytic domain only — PL25 vs CwlP RMSD = 0.63 Å (466 atoms); PL25 vs CwlT RMSD = 1.17 Å (65 atoms). (B) N-terminal glycoside-hydrolase domain only — PL25 Slt35-like vs CwlP SLT RMSD = 1.31 Å (401 atoms); PL25 vs CwlT GH25 RMSD = 7.92 Å (different fold). (C) Full catalytic core — inter-domain orientations differ between the three enzymes.

**Figure 6.** Synteny of contig_11 (PL25 plasmid) versus three closely related NCBI WGS contigs: *Ornithinibacillus contaminans* NZ_LDPV02000047.1, *Cytobacillus horneckiae* NZ_JARSFA010000044.1, and *Caldifermentibacillus hisashii* NZ_JBCITH010000002.1. 122 orthologous CDS detected by clinker; the conjugation operon (PL25 positions 2.5–22 kb) shows strong block synteny with 25+ orthologs in *O. contaminans* (45–77 % identity).

**Figure 7.** Integrated mobile-element architecture figure: (A) contig_11 gene map with M23 starred; (B) GC deviation from the chromosome mean (HGT signature, Δ = −5.6 %); (C) synteny ribbons to *O. contaminans*; (D) *O. contaminans* contig gene map.

**Figure 8.** Electrostatic comparison of PL25 M23 and reference cell-wall hydrolases. Surface electrostatic potentials are shown for PL25 M23, CwlP, CwlT, lysostaphin, and LytM, with red indicating negative potential, blue positive potential, and white neutral regions. PL25 M23 (red outline, bold) showed the lowest predicted pI (4.54), strongest negative net charge at pH 7 (−22.4), and high acidic residue content (D+E = 12.3%), consistent with an acidic surface profile associated with halophilic protein adaptation. Values were calculated for the modeled catalytic/mature regions used for structural comparison.

---

**Word counts:**
- Abstract: ~185 words
- Introduction: ~700 words
- Results: ~2,500 words (across 6 subsections)
- Discussion: ~1,150 words
- Methods: ~1,400 words
- **Total main text (Intro + Results + Discussion): ~4,350 words** ← within Nature Communications limit
