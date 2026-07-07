#
### Editorial Decision: **MAJOR REVISION** (Borderline Reject)

**Rationale:** The manuscript presents an interesting computational characterization of a prophage in an emerging ExPEC lineage. However, the title and framing are misleading—this is not a validation study. The lack of any wet-lab experiments, the ad hoc prophage detection method, and the overinterpretation of "active" and "functional" from purely computational evidence are significant issues. For FEMS Microbiology Letters, which publishes concise but experimentally grounded studies, this paper in its current form is not appropriate. The authors must either:
1. **Reframe** the paper as a computational characterization, remove "validation" from the title, and tone down claims; OR
2. **Add experimental validation** (prophage induction, hemolysis assay, MIC testing).

---

## PART 2 – COMPUTATIONAL EXPERIMENT AUDIT

### Experiment 1: Prophage Boundary Mapping

| Aspect | Assessment |
|--------|------------|
| **Why performed?** | Identify integrated prophage boundaries (attL/attR). |
| **Scientific justification?** | Yes, essential for prophage characterization. |
| **Addresses intended question?** | Partially. A 12-bp sliding window is ad hoc and may miss true boundaries. Established tools (PHASTER, Prophage Hunter) would be more reliable. |
| **State-of-the-art?** | No. Custom method is not benchmarked. |
| **Controls missing?** | No comparison to PHASTER or Prophage Hunter predictions. |
| **Validation missing?** | No long-read sequencing to confirm integration site. |
| **Alternative analyses expected?** | Reviewers will expect use of PHASTER or Prophage Hunter. |
| **Interpretation supported?** | Partially. Boundaries are plausible but not definitively validated. |

**Classification: Weak** (method is ad hoc and not benchmarked)

---

### Experiment 2: Functional Module Audit

| Aspect | Assessment |
|--------|------------|
| **Why performed?** | Determine if prophage is complete/defective. |
| **Scientific justification?** | Yes. |
| **Addresses intended question?** | Yes, confirms presence of all five modules. |
| **State-of-the-art?** | Adequate (manual curation based on annotations). |
| **Controls missing?** | No comparison to reference phages beyond DLP12. |
| **Validation missing?** | No induction experiment to confirm functional lytic machinery. |
| **Alternative analyses expected?** | Reviewers will expect confirmation of transcriptional activity (RNA-seq) or at least promoter prediction. |
| **Interpretation supported?** | **Weakly.** Presence of genes does not equal activity. "Active" is an overstatement without induction data. |

**Classification: Adequate** (but "active" claim is overinterpreted)

---

### Experiment 3: CRISPR Array Detection and Prophage Targeting

| Aspect | Assessment |
|--------|------------|
| **Why performed?** | Determine if host has immunological memory against prophage. |
| **Scientific justification?** | Yes. |
| **Addresses intended question?** | Yes. |
| **State-of-the-art?** | Yes (MinCED). |
| **Controls missing?** | No assessment of whether CRISPR array is functional (Cas genes not checked). |
| **Validation missing?** | No experimental confirmation of CRISPR interference. |
| **Alternative analyses expected?** | Reviewers will ask: are CRISPR-Cas genes intact? |
| **Interpretation supported?** | **Partially.** The absence of full-length spacers is clear. However, the conclusion that this explains stable maintenance is an overstatement—the prophage could also have anti-CRISPR genes or be maintained by other mechanisms. |

**Classification: Adequate**

---

### Experiment 4: CARD Resistance Phenotype Prediction

| Aspect | Assessment |
|--------|------------|
| **Why performed?** | Predict antibiotic resistance profile from genotype. |
| **Scientific justification?** | Yes. |
| **Addresses intended question?** | Yes. |
| **State-of-the-art?** | Yes (CARD database, ABRicate). |
| **Controls missing?** | No comparison to AMRFinderPlus or ResFinder (though previous paper may have this). |
| **Validation missing?** | No MIC testing to confirm predictions. |
| **Alternative analyses expected?** | Reviewers will expect resistance gene confirmation with multiple databases. |
| **Interpretation supported?** | **Partially.** Genotype-phenotype predictions are probabilistic, not definitive. Without MIC data, "Resistant" classifications are predictions, not validated phenotypes. |

**Classification: Good** (as a prediction; but "phenotype" in the title is misleading)





### Experiment 8: Comparative Synteny (clinker)

| Aspect | Assessment |
|--------|------------|
| **Why performed?** | Compare QA5221 prophage to DLP12. |
| **Scientific justification?** | Yes. |
| **Addresses intended question?** | Yes. |
| **State-of-the-art?** | Yes (clinker). |
| **Controls missing?** | No comparison to other lambdoid phages (e.g., lambda, HK97). |
| **Validation missing?** | No long-read sequencing to confirm synteny. |
| **Alternative analyses expected?** | Reviewers will expect broader comparison to multiple reference phages. |
| **Interpretation supported?** | Yes, the comparison to DLP12 effectively shows the QA5221 prophage is more complete. |

**Classification: Good**

---

## PART 3 – MISSING COMPUTATIONAL ANALYSES

### Critical Analyses

#### 1. Prophage Detection with Established Tools

| Parameter | Value |
|-----------|-------|
| **Scientific objective** | Validate prophage boundaries and completeness. |
| **Biological question** | Is the prophage accurately identified? |
| **Why reviewers expect it** | The custom sliding window method is ad hoc and not benchmarked. |
| **Expected improvement** | Critical; establishes confidence in prophage identification. |
| **Estimated runtime** | 1 hour |
| **Software** | PHASTER, Prophage Hunter, VirSorter2 |
| **CPU requirements** | 8 cores |
| **GPU requirements** | No |
| **RAM requirements** | 16 GB |
| **HPC required?** | No |
| **Difficulty** | Easy |
| **Estimated completion time** | 1 day |
| **Priority** | **Critical** |

---

#### 2. Prophage Induction Prediction

| Parameter | Value |
|-----------|-------|
| **Scientific objective** | Predict inducibility of the prophage. |
| **Biological question** | Is the prophage actually inducible, or is it a cryptic remnant? |
| **Why reviewers expect it** | "Active" requires evidence beyond gene presence. |
| **Expected improvement** | Critical; supports or refutes the "active" claim. |
| **Estimated runtime** | 1-2 hours |
| **Software** | Prophage induction predictors (e.g., PhageBoost), stress response gene analysis |
| **CPU requirements** | 8 cores |
| **GPU requirements** | No |
| **RAM requirements** | 16 GB |
| **HPC required?** | No |
| **Difficulty** | Moderate |
| **Estimated completion time** | 2 days |
| **Priority** | **Critical** |

---

#### 3. Check CRISPR-Cas Gene Integrity

| Parameter | Value |
|-----------|-------|
| **Scientific objective** | Determine if the CRISPR array is functional. |
| **Biological question** | Can the CRISPR system actually mount an immune response? |
| **Why reviewers expect it** | Presence of repeats/spacers is insufficient; Cas genes must be intact. |
| **Expected improvement** | High; validates CRISPR immune evasion conclusion. |
| **Estimated runtime** | 1 hour |
| **Software** | Prokka annotations, manual curation |
| **CPU requirements** | 4 cores |
| **GPU requirements** | No |
| **RAM requirements** | 8 GB |
| **HPC required?** | No |
| **Difficulty** | Easy |
| **Estimated completion time** | 1 day |
| **Priority** | **Critical** |

---

#### 4. Anti-CRISPR Gene Detection

| Parameter | Value |
|-----------|-------|
| **Scientific objective** | Identify if the prophage carries anti-CRISPR genes. |
| **Biological question** | How does the prophage evade CRISPR immunity? |
| **Why reviewers expect it** | Alternative explanation for CRISPR absence. |
| **Expected improvement** | Moderate; supports immune evasion mechanism. |
| **Estimated runtime** | 1 hour |
| **Software** | BLASTp against Acr database |
| **CPU requirements** | 4 cores |
| **RAM requirements** | 8 GB |
| **HPC required?** | No |
| **Difficulty** | Easy |
| **Estimated completion time** | 1 day |
| **Priority** | **High** |

---

#### 5. Broader Phage Synteny Comparison

| Parameter | Value |
|-----------|-------|
| **Scientific objective** | Place QA5221 prophage in context of other lambdoid phages. |
| **Biological question** | Is this a known phage variant or a novel prophage? |
| **Why reviewers expect it** | DLP12 comparison alone is insufficient. |
| **Expected improvement** | High; establishes novelty and relatedness. |
| **Estimated runtime** | 2-4 hours |
| **Software** | clinker, custom Python scripts |
| **CPU requirements** | 8 cores |
| **RAM requirements** | 32 GB |
| **HPC required?** | No |
| **Difficulty** | Moderate |
| **Estimated completion time** | 3 days |
| **Priority** | **High** |

---

#### 6. Phage Host Range Prediction

| Parameter | Value |
|-----------|-------|
| **Scientific objective** | Predict host range beyond the BLASTp inference. |
| **Biological question** | Is the phage strictly E. coli-specific? |
| **Why reviewers expect it** | Phage host range is biologically important. |
| **Expected improvement** | Moderate; supports zoonotic transmission context. |
| **Estimated runtime** | 1-2 hours |
| **Software** | BLASTp, PHACTS, or host prediction tools |
| **CPU requirements** | 8 cores |
| **RAM requirements** | 16 GB |
| **HPC required?** | No |
| **Difficulty** | Easy |
| **Estimated completion time** | 1 day |
| **Priority** | **Medium** |

---

#### 7. Pan-GWAS of Prophage Presence in ST354

| Parameter | Value |
|-----------|-------|
| **Scientific objective** | Determine if the prophage is lineage-specific or widespread. |
| **Biological question** | Is this prophage a marker of ST354? |
| **Why reviewers expect it** | Provides clinical and epidemiological context. |
| **Expected improvement** | High; strengthens significance. |
| **Estimated runtime** | 2-4 hours |
| **Software** | Panaroo, Scoary |
| **CPU requirements** | 16 cores |
| **RAM requirements** | 64 GB |
| **HPC required?** | No |
| **Difficulty** | Moderate |
| **Estimated completion time** | 3 days |
| **Priority** | **High** |

---

#### 8. Prophage Cargo Gene Evolutionary Analysis (dN/dS)

| Parameter | Value |
|-----------|-------|
| **Scientific objective** | Detect selection pressure on cargo genes. |
| **Biological question** | Are cargo genes under positive selection? |
| **Why reviewers expect it** | Supports adaptive significance. |
| **Expected improvement** | Moderate. |
| **Estimated runtime** | 1-2 hours |
| **Software** | PAML, HyPhy |
| **CPU requirements** | 8 cores |
| **RAM requirements** | 16 GB |
| **HPC required?** | No |
| **Difficulty** | Moderate |
| **Estimated completion time** | 2 days |
| **Priority** | **Medium** |

---

#### 9. Prophage-Encoded AMR Gene Analysis

| Parameter | Value |
|-----------|-------|
| **Scientific objective** | Check if the prophage carries any AMR genes. |
| **Biological question** | Does the prophage contribute to the resistome? |
| **Why reviewers expect it** | Prophage-mediated AMR is clinically relevant. |
| **Expected improvement** | Moderate; adds to resistance story. |
| **Estimated runtime** | 1 hour |
| **Software** | ABRicate, AMRFinderPlus |
| **CPU requirements** | 4 cores |
| **RAM requirements** | 8 GB |
| **HPC required?** | No |
| **Difficulty** | Easy |
| **Estimated completion time** | 1 day |
| **Priority** | **Medium** |

---

#### 10. Promoter and Regulatory Element Prediction

| Parameter | Value |
|-----------|-------|
| **Scientific objective** | Predict if the prophage genes are transcriptionally active. |
| **Biological question** | Is the prophage likely to be expressed? |
| **Why reviewers expect it** | Supports "active" claim. |
| **Expected improvement** | High; directly addresses the "active" question. |
| **Estimated runtime** | 1-2 hours |
| **Software** | BPROM, PromoScan, or manual curation |
| **CPU requirements** | 4 cores |
| **RAM requirements** | 8 GB |
| **HPC required?** | No |
| **Difficulty** | Moderate |
| **Estimated completion time** | 2 days |
| **Priority** | **High** |

---

## PART 4 – EXPERIMENTAL GAPS

### Major Claims and Support Levels

| Claim | Support | Justification |
|-------|---------|---------------|
| **Prophage is "active"** | **Weakly supported** | Gene presence ≠ activity. No induction experiment, no RNA-seq, no promoter prediction. "Active" is an overstatement. |
| **Prophage is "novel"** | **Partially supported** | Synteny differs from DLP12, but no comparison to other characterized phages. Broader comparison needed. |
| **Prophage carries enterohemolysin cytolysin** | **Fully supported** | ESMFold pLDDT 88.4 + Foldseek to characterized hemolysin. Strong structural evidence. |
| **CRISPR immune gap explains prophage stability** | **Partially supported** | No full-length spacers. But Cas gene integrity not checked; anti-CRISPR genes not ruled out. |
| **Resistance phenotype predicted** | **Fully supported** | CARD hits are high-confidence. However, these are predictions, not validated phenotypes. |
| **Biocide co-selection mechanism** | **Fully supported** | Physical linkage of qac genes with AMR genes in integron is well-supported. |
| **Phage is lambdoid Siphoviridae with E. coli host range** | **Fully supported** | Dual BLASTp + Foldseek strongly supports. |
| **ST354 is a zoonotic reservoir** | **Partially supported** | References cited, but no data showing zoonotic transmission in this study. |
| **Lineage-wide distribution of virulence markers** | **Weakly supported** | Mentioned but not shown in main text (relegated to Supplementary). |

### Unsupported Claims That Need Softening

1. **"Active prophage"** → Should be rephrased as "genomically intact prophage with complete lytic and lysogeny modules" or "putatively active prophage."

2. **"Validation"** in the title → Should be changed to "Characterization" or "Computational Prediction." There is no validation.

3. **"Resistant phenotype"** → Should be "predicted resistance profile" or "genotype-inferred resistance."

4. **"CRISPR immune evasion"** → Should be "CRISPR spacer absence" or "CRISPR immune gap." Evasion implies an active mechanism (e.g., anti-CRISPR), which is not shown.

5. **"Functional modules"** → Should be "predicted functional modules" unless transcriptional activity is confirmed.

---

## PART 5 – METHODOLOGICAL WEAKNESSES

### Ranked Weaknesses (Most Severe to Least Severe)

| Rank | Weakness | Severity |
|------|----------|----------|
| 1 | **"Validation" in title with NO wet-lab validation.** This is misleading. | **Critical** |
| 2 | **Ad hoc prophage detection method.** Sliding 12-bp window is not validated against known methods. | **Critical** |
| 3 | **No prophage induction experiment** (e.g., mitomycin C treatment). | **Critical** |
| 4 | **No RNA-seq or RT-PCR** to confirm transcriptional activity. | **Critical** |
| 5 | **CRISPR Cas gene integrity not checked.** Spacers without Cas = no function. | **High** |
| 6 | **No anti-CRISPR gene search** (alternative explanation for immune evasion). | **High** |
| 7 | **No comparison to other lambdoid phages** beyond DLP12. | **High** |
| 8 | **No MIC testing** to validate resistance phenotype predictions. | **High** |
| 9 | **No hemolysis assay** to validate enterohemolysin function. | **High** |
| 10 | **No long-read sequencing** to confirm prophage integration and synteny. | **High** |
| 11 | **Lineage-wide analyses relegated to Supplementary**—main text lacks context. | **Medium** |
| 12 | **No statistical analysis** of spacer presence/absence (binary comparisons). | **Medium** |
| 13 | **No code repository** for the custom prophage detection script. | **Medium** |
| 14 | **No promoter prediction** to assess transcriptional activity. | **Medium** |
| 15 | **No comparison of predicted vs. actual resistance phenotypes** in literature. | **Medium** |
| 16 | **No prophage prevalence analysis** in the ST354 cohort. | **Medium** |
| 17 | **No dN/dS analysis** of cargo genes. | **Medium** |
| 18 | **No search for phage-related hypothetical proteins** in the prophage beyond the 5 candidates. | **Medium** |
| 19 | **No assessment of false positive rates** for Foldseek or BLASTp. | **Low** |
| 20 | **Figures lack scale bars** in synteny diagram. | **Low** |

---

## PART 6 – REVIEWER CRITICISMS

### Major Criticisms (20)

1. **"Active" prophage claim is unsupported**
   - *Why:* Gene presence ≠ function. No induction assay.
   - *Justified:* Yes, this is a critical overstatement.
   - *How to fix:* Remove "active" or add caveats. Run mitomycin C induction.
   - *Requires:* Writing OR wet-lab experiment.

2. **Title promises "Validation" but provides none**
   - *Why:* No wet-lab experiments.
   - *Justified:* Yes, misleading.
   - *How to fix:* Change title to "Computational Characterization."
   - *Requires:* Writing only.

3. **Prophage detection method is ad hoc**
   - *Why:* 12-bp sliding window is not a standard method.
   - *Justified:* Yes.
   - *How to fix:* Use PHASTER or Prophage Hunter.
   - *Requires:* New computational analysis.

4. **No Cas gene integrity check**
   - *Why:* CRISPR immunity requires intact Cas proteins.
   - *Justified:* Yes.
   - *How to fix:* Annotate Cas genes.
   - *Requires:* New computational analysis.

5. **No anti-CRISPR gene search**
   - *Why:* Alternative explanation for immune evasion.
   - *Justified:* Yes.
   - *How to fix:* BLASTp against Acr database.
   - *Requires:* New computational analysis.

6. **No RNA-seq or RT-PCR to confirm expression**
   - *Why:* "Active" requires expression.
   - *Justified:* Yes.
   - *How to fix:* Add RNA-seq data or acknowledge limitation.
   - *Requires:* Wet-lab OR writing.

7. **No MIC testing to validate resistance predictions**
   - *Why:* CARD predictions are not validated.
   - *Justified:* Yes.
   - *How to fix:* Perform MIC testing or add caveats.
   - *Requires:* Wet-lab OR writing.

8. **No hemolysis assay for enterohemolysin**
   - *Why:* Structural prediction is not functional validation.
   - *Justified:* Yes.
   - *How to fix:* Perform hemolysis assay or add caveats.
   - *Requires:* Wet-lab OR writing.

9. **No comparison to other lambdoid phages**
   - *Why:* DLP12 comparison alone is insufficient.
   - *Justified:* Yes.
   - *How to fix:* Add comparisons to lambda, HK97, etc.
   - *Requires:* New computational analysis.

10. **Lineage-wide context is missing from main text**
    - *Why:* ST354 context is important.
    - *Justified:* Yes.
    - *How to fix:* Move key lineage-wide results to main text.
    - *Requires:* Writing.

11. **No prophage prevalence analysis**
    - *Why:* Is this prophage unique to QA5221 or common in ST354?
    - *Justified:* Yes.
    - *How to fix:* Analyze prophage presence in the 32-genome cohort.
    - *Requires:* New computational analysis.

12. **No promoter prediction for prophage genes**
    - *Why:* Supports transcriptional activity.
    - *Justified:* Yes.
    - *How to fix:* Use BPROM or PromoScan.
    - *Requires:* New computational analysis.

13. **No dN/dS analysis of cargo genes**
    - *Why:* Detects positive selection.
    - *Justified:* Yes.
    - *How to fix:* Perform dN/dS if homologs are available.
    - *Requires:* New computational analysis.

14. **No TEM imaging to confirm phage morphology**
    - *Why:* Provides visual confirmation.
    - *Justified:* Partially (beyond scope).
    - *How to fix:* Acknowledge as future direction.
    - *Requires:* Writing.

15. **No long-read sequencing**
    - *Why:* Confirm integration site and synteny.
    - *Justified:* Yes.
    - *How to fix:* Acknowledge limitation; recommend long-read sequencing.
    - *Requires:* Writing.

16. **Resistance phenotype predictions are overinterpreted**
    - *Why:* Genotype ≠ phenotype.
    - *Justified:* Yes.
    - *How to fix:* Use "predicted" language throughout.
    - *Requires:* Writing.

17. **No multiple comparison correction**
    - *Why:* Statistical rigor.
    - *Justified:* Yes.
    - *How to fix:* Apply Bonferroni or FDR correction.
    - *Requires:* Writing.

18. **Figure 3 synteny diagram is cluttered**
    - *Why:* Hard to read.
    - *Justified:* Yes.
    - *How to fix:* Simplify or add zoomed panels.
    - *Requires:* Figures.

19. **No code repository**
    - *Why:* Reproducibility.
    - *Justified:* Yes.
    - *How to fix:* Add GitHub link.
    - *Requires:* Writing.

20. **"Functional" modules claim is unsupported**
    - *Why:* Gene presence ≠ function.
    - *Justified:* Yes.
    - *How to fix:* Use "predicted functional modules."
    - *Requires:* Writing.

---

### Minor Criticisms (30)

1. **"Extraintestinal" spelled wrong in some places** → Fix.
2. **DLP12 accession wrong (NC_001799.1 is bacteriophage lambda)** → Fix.
3. **Methods refer to Supplementary Material for baseline data** → But Supplementary not included.
4. **No line numbers** → Add for review.
5. **Reference formatting inconsistent** → Fix.
6. **Figure 1 missing scale** → Add.
7. **Figure 2 donut chart not described in Results** → Add description.
8. **Spacer sequences not listed** → Add in Supplementary.
9. **No coordinates for CRISPR array** → Add.
10. **No discussion of prophage excision sites** → Add.
11. **No discussion of phage termini** → Add.
12. **"Candidate 1-5" names are generic** → Add functional names.
13. **No discussion of phage packaging mechanism** → Add.
14. **No discussion of phage integration mechanism** → Add.
15. **No discussion of phage receptor specificity** → Add.
16. **CARD version not specified** → Add.
17. **Foldseek version not specified** → Add.
18. **ESMFold version not specified** → Add.
19. **No discussion of phage-host coevolution** → Add.
20. **No discussion of clinical implications** → Expand.
21. **No discussion of limitations of bioinformatics predictions** → Add.
22. **No discussion of database coverage** → Add.
23. **No discussion of false positives in Foldseek** → Add.
24. **No discussion of phage annotation accuracy** → Add.
25. **Abbreviations not defined** → Fix.
26. **Figure 5 "I" vs "R" labels confusing** → Clarify.
27. **No discussion of biocide resistance clinical relevance** → Add.
28. **No discussion of integron mobility** → Add.
29. **No discussion of phage cargo evolution** → Add.
30. **Conclusion too brief** → Expand.

---

## PART 7 – COMPUTATIONAL COST–BENEFIT ANALYSIS

| Analysis | Scientific Benefit | Reviewer Impact | Acceptance Improvement | Cost | Effort | Time | Priority |
|----------|-------------------|-----------------|----------------------|------|--------|------|----------|
| PHASTER/Prophage Hunter validation | **Very High** | **Critical** | **High** | Low | Easy | 1 day | **Critical** |
| CRISPR Cas integrity check | **High** | **Critical** | **High** | Low | Easy | 1 day | **Critical** |
| Promoter prediction | **High** | **High** | **High** | Low | Moderate | 2 days | **Critical** |
| Anti-CRISPR gene search | **High** | **High** | **High** | Low | Easy | 1 day | **High** |
| Broader phage synteny comparison | **High** | **High** | **High** | Low | Moderate | 3 days | **High** |
| Prophage prevalence in ST354 | **High** | **High** | **High** | Low | Moderate | 3 days | **High** |
| Pan-GWAS of prophage presence | **High** | **High** | **High** | Medium | Moderate | 3 days | **High** |
| Prophage induction prediction | **High** | **High** | **High** | Low | Moderate | 2 days | **High** |
| Code repository | **High** | **High** | **High** | Low | Easy | 1 day | **High** |
| Reframe as computational characterization | **Very High** | **Critical** | **High** | Low | Easy | 1 day | **Critical** |
| dN/dS analysis | Medium | Medium | Medium | Low | Moderate | 2 days | **Medium** |
| Prophage AMR gene search | Medium | Medium | Medium | Low | Easy | 1 day | **Medium** |
| Host range prediction | Medium | Medium | Low | Low | Easy | 1 day | **Medium** |
| TEM imaging | High | High | High | Very High | Very Hard | Months | **Low** (wet-lab) |
| Long-read sequencing | High | High | High | Very High | Very Hard | Months | **Low** (wet-lab) |
| MIC testing | High | High | High | High | Moderate | 1 week | **Low** (wet-lab) |

---

## PART 8 – WRITING VS SCIENCE

### Criticisms Solvable by Writing Only

| Criticism | Solution |
|-----------|----------|
| "Validation" in title | Change to "Computational Characterization" |
| "Active" overstatement | Use "predicted active" or "genomically intact" |
| "Functional" modules | Use "predicted functional modules" |
| "Resistant phenotype" | Use "predicted resistance profile" |
| "CRISPR immune evasion" | Use "CRISPR spacer absence" or "immune gap" |
| No line numbers | Add |
| Reference formatting | Fix |
| Abbreviation inconsistencies | Fix |
| No limitations section | Add |
| No clinical implications | Expand |
| No future directions | Add |

### Criticisms Requiring New Computational Analysis

| Criticism | Analysis Required |
|-----------|-------------------|
| Ad hoc prophage detection | PHASTER/Prophage Hunter validation |
| No Cas gene check | Annotate Cas genes |
| No anti-CRISPR search | BLASTp against Acr database |
| No promoter prediction | BPROM or PromoScan |
| No broader phage comparison | More synteny comparisons |
| No prophage prevalence | Analyze 32-genome cohort |
| No Pan-GWAS | Perform Pan-GWAS |
| No dN/dS | Perform dN/dS |

### Criticisms Requiring Wet-Lab Validation

| Criticism | Solution |
|-----------|----------|
| Prophage not confirmed as inducible | Mitomycin C induction |
| Enterohemolysin not confirmed as functional | Hemolysis assay |
| Resistance phenotype not validated | MIC testing |
| CRISPR not confirmed as functional | CRISPR interference assay |
| Synteny not confirmed | Long-read sequencing |
| Phage morphology not confirmed | TEM imaging |

---

## PART 9 – FINAL PRIORITY ROADMAP

### Tier 1 – Essential Before Submission

| Item | Reviewer Impact | Difficulty | Time | HPC? | Quality Increase |
|------|----------------|------------|------|------|------------------|
| **Reframe manuscript as computational characterization** | Critical | Easy | 1 day | No | Very High |
| **Remove "validation" from title** | Critical | Easy | 1 hour | No | Very High |
| **Run PHASTER/Prophage Hunter** | Critical | Easy | 1 day | No | High |
| **Check CRISPR Cas gene integrity** | Critical | Easy | 1 day | No | High |
| **Add promoter prediction** | High | Moderate | 2 days | No | High |
| **Add anti-CRISPR gene search** | High | Easy | 1 day | No | High |
| **Add broader phage synteny comparison** | High | Moderate | 3 days | No | High |
| **Tone down "active" language** | Critical | Easy | 2 hours | No | High |
| **Add limitations section** | Critical | Easy | 1 day | No | High |

**Estimated total time:** 1-2 weeks

---

### Tier 2 – Strongly Recommended

| Item | Reviewer Impact | Difficulty | Time | HPC? | Quality Increase |
|------|----------------|------------|------|------|------------------|
| **Prophage prevalence in ST354** | High | Moderate | 3 days | No | High |
| **Pan-GWAS** | High | Moderate | 3 days | No | High |
| **Prophage induction prediction** | High | Moderate | 2 days | No | High |
| **Code repository** | High | Easy | 1 day | No | High |
| **Move lineage-wide results to main text** | Medium | Easy | 1 day | No | Medium |
| **Improve Figure 3 (synteny)** | Medium | Easy | 1 day | No | Medium |

**Estimated total time:** 2 weeks

---

### Tier 3 – Optional Improvements

| Item | Reviewer Impact | Difficulty | Time | HPC? | Quality Increase |
|------|----------------|------------|------|------|------------------|
| dN/dS analysis | Medium | Moderate | 2 days | No | Medium |
| Prophage AMR gene search | Medium | Easy | 1 day | No | Medium |
| Host range prediction | Medium | Easy | 1 day | No | Low |
| Improve figures | Medium | Easy | 2 days | No | Medium |

**Estimated total time:** 1 week

---

### Tier 4 – Unnecessary (for this submission)

| Item | Reason |
|------|--------|
| Mitomycin C induction | Best as follow-up wet-lab study |
| Hemolysis assay | Best as follow-up wet-lab study |
| MIC testing | Best as follow-up wet-lab study |
| TEM imaging | Best as follow-up wet-lab study |
| Long-read sequencing | Best as follow-up wet-lab study |

---

## PART 10 – FINAL EDITOR REPORT

### 1. What are the 20 most important weaknesses?

1. **Misleading title** — "Validation" implies wet-lab experiments; none are performed.
2. **No wet-lab validation** — Title promises validation but delivers only computational predictions.
3. **Ad hoc prophage detection method** — 12-bp sliding window is not validated against established tools.
4. **"Active" prophage claim unsupported** — No induction experiment, no RNA-seq, no promoter prediction.
5. **CRISPR Cas gene integrity not checked** — Spacers without Cas proteins are non-functional.
6. **No anti-CRISPR gene search** — Alternative explanation for CRISPR absence.
7. **No broader phage comparison** — DLP12 comparison alone is insufficient.
8. **No promoter prediction** — Supports transcriptional activity assessment.
9. **No prophage prevalence analysis** — Is this prophage unique or common in ST354?
10. **No Pan-GWAS** — Links prophage to clinical phenotypes.
11. **No dN/dS analysis** — Detects selection on cargo genes.
12. **Resistance phenotype predictions not validated** — No MIC testing.
13. **Enterohemolysin not functionally validated** — No hemolysis assay.
14. **Lineage-wide context relegated to Supplementary** — Weakens main text.
15. **No code repository** — Hinders reproducibility.
16. **No statistical analysis** — No multiple comparison correction.
17. **"Functional" modules claim** — Gene presence ≠ function.
18. **No TEM imaging** — No visual confirmation of phage morphology.
19. **No long-read sequencing** — Cannot confirm integration site.
20. **Conclusion too brief** — Does not synthesize findings effectively.

---

### 2. Which missing computational analyses should definitely be performed?

**Essential:**
1. ✅ PHASTER/Prophage Hunter validation
2. ✅ CRISPR Cas gene integrity check
3. ✅ Promoter prediction
4. ✅ Anti-CRISPR gene search
5. ✅ Broader phage synteny comparison
6. ✅ Reframe manuscript (remove "validation," tone down "active")
7. ✅ Add limitations section

**Strongly recommended:**
8. ✅ Prophage prevalence in ST354
9. ✅ Pan-GWAS
10. ✅ Prophage induction prediction
11. ✅ Code repository

---

### 3. Which suggested analyses are unnecessary?

1. **TEM imaging** → Best as follow-up wet-lab study.
2. **Long-read sequencing** → Best as follow-up study.
3. **Mitomycin C induction** → Best as follow-up wet-lab study.
4. **Hemolysis assay** → Best as follow-up wet-lab study.
5. **MIC testing** → Best as follow-up wet-lab study.

---

### 4. Are the molecular dynamics simulations sufficient?

**N/A.** This manuscript does not contain MD simulations. This is a notable omission for the enterohemolysin candidate, given that the companion paper (Paper 1) included MD simulations that were criticized. The absence of MD here is actually appropriate, as the focus is on prophage characterization rather than protein dynamics. However, it means the enterohemolysin is only structurally characterized by ESMFold and Foldseek, not by MD.

---

### 5. Are the structural analyses sufficient?

**Partially.** ESMFold pLDDT (88.4) and Foldseek hits provide strong structural evidence. However:
- No comparison to AlphaFold3.
- No MD simulations (though not essential for this paper's scope).
- No experimental validation (CD/crystallography).

**Recommendation:** The structural analyses are adequate for a computational characterization but should be acknowledged as predictive.

---

### 6. Is the pangenomic analysis sufficient?

**Partially.** The paper states that baseline pangenomic analysis is in the Supplementary Material. This is problematic because:
- The main text lacks lineage-wide context.
- Reviewers may not read the Supplementary thoroughly.

**Recommendation:** Move key pangenomic results (e.g., open pangenome, accessory partition sizes) to the main text.

---

### 7. Is the protein language model analysis sufficient?

**N/A.** This paper does not use PLM analysis (ESM-2, UMAP, etc.). The structural prediction uses ESMFold, which is a PLM-based predictor, but there is no latent-space analysis.

---

### 8. Are the evolutionary analyses sufficient?

**Partially.** The phage taxonomic classification and synteny comparison are well done. However:
- No dN/dS analysis.
- No broader evolutionary context (e.g., phylogenetic tree of the prophage).

**Recommendation:** Add dN/dS analysis if homologs are available.

---

### 9. Are the resistance claims sufficiently supported?

**Partially.** The CARD-based resistance profile is a prediction, not a validated phenotype. The physical linkage of biocide and AMR genes is well-supported and mechanistically plausible.

**Recommendation:** Use "predicted resistance profile" language. Add caveats about lack of MIC validation.

---

### 10. Are the virulence claims sufficiently supported?

**Partially.** The enterohemolysin structural prediction is strong (pLDDT 88.4, Foldseek to characterized hemolysin). However:
- No hemolysis assay.
- No confirmation that the protein is expressed.
- No functional characterization.

**Recommendation:** Use "putative" or "predicted" language. Add caveats.

---

### 11. Which claims should be softened?

| Current Claim | Proposed Revision |
|---------------|-------------------|
| "Validation" (title) | "Computational Characterization" |
| "Active prophage" | "Putatively active prophage" or "genomically intact prophage" |
| "Functional modules" | "Predicted functional modules" |
| "Resistant phenotype" | "Predicted resistance profile" |
| "CRISPR immune evasion" | "CRISPR spacer absence" or "CRISPR immune gap" |
| "Pore-forming cytolysin" | "Predicted pore-forming cytolysin" |
| "Functional enterohemolysin" | "Putative enterohemolysin" |
| "Validated" | "Characterized" or "Computed" |

---

### 12. Which figures should be improved?

| Figure | Issue | Improvement |
|--------|-------|-------------|
| **Figure 1** (Module map) | Good | Add scale bar |
| **Figure 2** (Donut chart) | Not described in text | Add description |
| **Figure 3** (Synteny) | Cluttered | Simplify or add zoomed panels |
| **Figure 4** (CRISPR) | Good | Add spacer sequences |
| **Figure 5** (Phenotype) | Good | Clarify "I" vs "R" labels |
| **Figure 6** (pLDDT) | Good | Add 3D model visualization |
| **Figure 7** (Taxonomy) | Good | Add confidence intervals |

---

### 13. Which supplementary analyses should be added?

**Must add:**
- Supplementary Table: PHASTER/Prophage Hunter comparison
- Supplementary Table: CRISPR Cas genes
- Supplementary Table: Anti-CRISPR gene hits
- Supplementary Table: Promoter predictions
- Supplementary Figure: Broader phage synteny comparison
- Supplementary Table: Prophage prevalence in ST354
- Supplementary Figure: Phage genome map with gene annotations

**Strongly recommended:**
- Supplementary Table: dN/dS results
- Supplementary Table: Pan-GWAS results
- Supplementary Table: Phage host range predictions
- Supplementary File: Code repository

---

### 14. What is the single biggest weakness of the manuscript?

**The title is misleading.** The manuscript claims to perform "Validation" but contains **no wet-lab experiments whatsoever**. This is a computational characterization study, not a validation study. The title must be changed to accurately reflect the content. Additionally, the "active prophage" claim is unsupported without induction experiments or transcriptional data.

**Recommendation:** Change title to "Computational Characterization of a Putatively Active Prophage and Resistance Hotspots in the Extraintestinal Pathogenic Escherichia coli ST354 Lineage" or similar.

---

### 15. If you were the handling editor, what revisions would you require before sending the paper for external review?

**Required revisions (must be completed before submission):**

1. ✅ **Change the title** — Remove "Validation"; replace with "Computational Characterization" or similar.
2. ✅ **Remove "active" and "functional" overstatements** — Use "putatively active" and "predicted functional" throughout.
3. ✅ **Run PHASTER or Prophage Hunter** to validate prophage boundaries.
4. ✅ **Check CRISPR Cas gene integrity** — Confirm the array is functional.
5. ✅ **Add anti-CRISPR gene search** — Rule out alternative immune evasion.
6. ✅ **Add promoter prediction** — Support transcriptional activity inference.
7. ✅ **Add broader phage synteny comparison** — Compare to lambda, HK97, etc.
8. ✅ **Add prophage prevalence analysis** — Is it common in ST354?
9. ✅ **Add limitations section** — Acknowledge: no wet-lab validation, no induction, no MIC, no hemolysis assay.
10. ✅ **Add code repository** — For the custom prophage detection script.
11. ✅ **Move key lineage-wide results to main text** — Not just Supplementary.
12. ✅ **Tone down resistance claims** — Use "predicted" language.
13. ✅ **Tone down virulence claims** — Use "putative" language.

**Recommended revisions (optional but beneficial):**
14. Add dN/dS analysis.
15. Add Pan-GWAS.
16. Add prophage induction prediction.
17. Add host range prediction.

---

## FINAL RANKED ACTION PLAN

### Highest Impact Actions (Tier 1)

| Rank | Action | Expected Impact | Time | Difficulty |
|------|--------|----------------|------|------------|
| 1 | **Change title** (remove "Validation") | Critical | 1 hour | Easy |
| 2 | **Reframe as computational characterization** | Critical | 1 day | Easy |
| 3 | **Tone down "active" and "functional" claims** | Critical | 2 hours | Easy |
| 4 | **Run PHASTER/Prophage Hunter** | High | 1 day | Easy |
| 5 | **Check CRISPR Cas gene integrity** | High | 1 day | Easy |
| 6 | **Add anti-CRISPR gene search** | High | 1 day | Easy |
| 7 | **Add promoter prediction** | High | 2 days | Moderate |
| 8 | **Add broader phage synteny comparison** | High | 3 days | Moderate |
| 9 | **Add limitations section** | High | 1 day | Easy |
| 10 | **Add code repository** | High | 1 day | Easy |

### Strongly Recommended Actions (Tier 2)

| Rank | Action | Expected Impact | Time | Difficulty |
|------|--------|----------------|------|------------|
| 11 | **Prophage prevalence in ST354** | High | 3 days | Moderate |
| 12 | **Pan-GWAS** | High | 3 days | Moderate |
| 13 | **Prophage induction prediction** | High | 2 days | Moderate |
| 14 | **Move lineage-wide results to main text** | Medium | 1 day | Easy |

### Optional Actions (Tier 3)

| Rank | Action | Expected Impact | Time | Difficulty |
|------|--------|----------------|------|------------|
| 15 | dN/dS analysis | Medium | 2 days | Moderate |
| 16 | Host range prediction | Medium | 1 day | Easy |
| 17 | Improve figures | Medium | 2 days | Easy |

### Unnecessary Actions (Tier 4)

| Rank | Action | Reason |
|------|--------|--------|
| 18 | Mitomycin C induction | Best as follow-up wet-lab study |
| 19 | Hemolysis assay | Best as follow-up wet-lab study |
| 20 | MIC testing | Best as follow-up wet-lab study |
| 21 | TEM imaging | Best as follow-up wet-lab study |
| 22 | Long-read sequencing | Best as follow-up study |

---

## SUMMARY OF REQUIRED ACTIONS

To make the manuscript publication-ready for **FEMS Microbiology Letters**, the following **must** be done:

### Critical (Must Do)
1. ✅ Change title (remove "Validation")
2. ✅ Reframe as computational characterization
3. ✅ Tone down "active" and "functional" language
4. ✅ Run PHASTER/Prophage Hunter
5. ✅ Check CRISPR Cas genes
6. ✅ Add anti-CRISPR search
7. ✅ Add promoter prediction
8. ✅ Add broader phage synteny comparison
9. ✅ Add limitations section
10. ✅ Add code repository

### Important (Strongly Recommended)
11. ✅ Prophage prevalence in ST354
12. ✅ Pan-GWAS
13. ✅ Prophage induction prediction
14. ✅ Move lineage-wide results to main text

### Optional (If Time Permits)
15. ✅ dN/dS analysis
16. ✅ Host range prediction

**Estimated total time for Critical + Important actions:** 2-3 weeks

**Recommendation:** The manuscript has interesting findings but is not publication-ready in its current form. The title is misleading, the computational methods are ad hoc in places, and the lack of wet-lab validation means the claims must be significantly toned down. With the revisions outlined above, the manuscript would be suitable for FEMS Microbiology Letters as a **Computational Characterization** study. Without these revisions, the manuscript would likely be rejected for overclaiming and methodological weaknesses.