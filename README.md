# kechris-lab

### Introduction

Chronic obstructive pulmonary disease (COPD) is a leading cause of
illness and death in the United States [[\[1\]](#loc-9)]{#loc-1
style="text-decoration: underline" role="doc-biblioref"}. It comprises
two conditions, chronic bronchitis and emphysema which frequently
co-occur [[\[1\]](#loc-9)]{style="text-decoration: underline"
role="doc-biblioref"}, and is a heterogenous disease with respiratory
symptoms that lead to airflow limitation due to abnormalities of the
airway and/or aveoli [[\[2\]](#loc-10)]{#loc-2
style="text-decoration: underline" role="doc-biblioref"}. Inhalation of
noxious particles or gas is the primary etiologic driver, with cigarette
smoke being the most important environmental risk factor
[[\[1\]](#loc-9)]{style="text-decoration: underline"
role="doc-biblioref"}. Adjusting for age and sex, smokers are 3.18 times
more likely to develop COPD than non-smokers [[\[3\]](#loc-11)]{#loc-3
style="text-decoration: underline" role="doc-biblioref"} yet fewer than
50% of smokers develop the disease [[\[4\]](#loc-12)]{#loc-4
style="text-decoration: underline" role="doc-biblioref"}. Males are also
1.48 times more likely to develop the disease compared to females, and
the odds of developing COPD increase with age
[[\[3\]](#loc-11)]{style="text-decoration: underline"
role="doc-biblioref"}. COPD's incomplete penetrance suggests that
downstream responses to smoking, rather than exposure alone, may mediate
progression to clinical COPD.

Although the lung is the main organ affected, the clinical burden of
COPD extends beyond the respiratory system. Comorbidities such as muscle
wasting, cardiovascular disease, osteoporosis, and depression are
consistent with generalized metabolic disturbance
[[\[5\]](#loc-13)]{#loc-5 style="text-decoration: underline"
role="doc-biblioref"}. Previous studies have found hundreds of
metabolites associated with lung function measures like Forced
Expiratory Volume (FEV) and Forced Vital Capacity (FVC), and
multi-metabolite signatures can predict COPD severity
[[\[5\]](#loc-13)]{#loc-6 style="text-decoration: underline"
role="doc-biblioref"},
[[\[6\]](#loc-14)]{style="text-decoration: underline"
role="doc-biblioref"}. Nicotine has a short half-life so metabolomic
profiles of current smokers typically reflect its more stable
metabolite, cotinine (and related derivatives), which serves as a
reliable biomarker of recent tobacco exposure [[\[7\]](#loc-15)]{#loc-7
style="text-decoration: underline" role="doc-biblioref"}. Many analyses,
however, miss disease-relevant coordination among metabolites and focus
on single metabolites. Smoking likely reorganizes the inherently
networked metabolite relationships through processes like inflammation
and oxidative stress, revealing metabolite relationships that may
persist or rewire after smoking cessation.

What remains poorly understood is whether (and how) metabolite
interactions differ between current, and former smokers, and how these
cessation-associated changes differ between individuals with and without
COPD. *I hypothesize that active smoking induces alterations in the
structure of metabolic interactions, and smoking cessation partially
reverses some relationships.* To test this, I will use a dataset with
smoking status, COPD staging, and approximately 1,000 metabolites
measured on the Metabolon platform to infer condition-specific
metabolite networks using a graphical lasso model that finds shared and
condition-specific structure via fusion and condition-adaptive penalties
(Condition-Adaptive Fused Graphical Lasso
(CFGL))[[\[8\]](#loc-16)]{#loc-8 style="text-decoration: underline"
role="doc-biblioref"}. This approach extends gene co-expression concepts
to metabolomics by modeling conditional dependencies among metabolites.
The condition-specific network structure will reveal how biological
processes are differentially active across smoking status and disease
status.

The results of this research will yield targetable pathways for
prevention or stage-tailored intervention in smoking individuals with or
at-risk for COPD.

**References**

::: {.section role="doc-bibliography"}
- [[[[\[1\]](#loc-1)]{style="text-decoration: underline"
  role="doc-backlink"}]{.prefix} S. Antwi, S. E. Steck, and K. Heidari,
  "Association between prevalence of chronic obstructive pulmonary
  disease and health-related quality of life, South Carolina, 2011,"
  *Prev. Chronic Dis.*, vol. 10, no. 130192, p. E215, Dec.
  2013.]{#loc-9}
- [[[[\[2\]](#loc-2)]{style="text-decoration: underline"
  role="doc-backlink"}]{.prefix} Y. Zhuang *et al.*, "Deep learning on
  graphs for multi-omics classification of COPD," *PLoS One*, vol. 18,
  no. 4, p. e284563, Apr. 2023.]{#loc-10}
- [[[[\[3\]](#loc-3)]{style="text-decoration: underline"
  role="doc-backlink"}]{.prefix} W. H. Thompson and S. St-Hilaire,
  "Prevalence of chronic obstructive pulmonary disease and tobacco use
  in veterans at Boise Veterans Affairs Medical Center," *Respir. Care*,
  vol. 55, no. 5, pp. 555--560, May 2010.]{#loc-11}
- [[[[\[4\]](#loc-4)]{style="text-decoration: underline"
  role="doc-backlink"}]{.prefix} A. Agustí *et al.*, "Global initiative
  for chronic obstructive lung disease 2023 report: GOLD executive
  summary," *Eur. Respir. J.*, vol. 61, no. 4, p. 2300239, Apr.
  2023.]{#loc-12}
- [[[[\[5\]](#loc-5)]{style="text-decoration: underline"
  role="doc-backlink"}]{.prefix} S. Godbole *et al.*, "A metabolomic
  severity score for airflow obstruction and emphysema," *Metabolites*,
  vol. 12, no. 5, p. 368, Apr. 2022.]{#loc-13}
- [[[[\[6\]](#loc-6)]{style="text-decoration: underline"
  role="doc-backlink"}]{.prefix} Global Initiative for Chronic
  Obstructive Lung Disease, *Pocket guide to COPD diagnosis, management,
  and prevention*. Global Initiative for Chronic Obstructive Lung
  Disease, 2025.]{#loc-14}
- [[[[\[7\]](#loc-7)]{style="text-decoration: underline"
  role="doc-backlink"}]{.prefix} M. A. Miller *et al.*, "Gene and
  metabolite time-course response to cigarette smoking in mouse lung and
  plasma," *PLoS One*, vol. 12, no. 6, p. e178281, June 2017.]{#loc-15}
- [[[[\[8\]](#loc-8)]{style="text-decoration: underline"
  role="doc-backlink"}]{.prefix} Y. Lyu *et al.*, "Condition-adaptive
  fused graphical lasso (CFGL): An adaptive procedure for inferring
  condition-specific gene co-expression network," *PLoS Comput. Biol.*,
  vol. 14, no. 9, p. e1006436, Sept. 2018.]{#loc-16}
:::

### Background and Significance

Chronic obstructive pulmonary disease (COPD) is a leading cause of
illness and death in the United States
[[\[1\]](#loc-9)]{style="text-decoration: underline"
role="doc-biblioref"}. It comprises two conditions, chronic bronchitis
and emphysema which frequently co-occur
[[\[1\]](#loc-9)]{style="text-decoration: underline"
role="doc-biblioref"}, and is a heterogenous disease with respiratory
symptoms that lead to airflow limitation due to abnormalities of the
airway and/or aveoli
[[\[2\]](#loc-10)]{style="text-decoration: underline"
role="doc-biblioref"}. Inhalation of noxious particles or gas is the
primary etiologic driver, with cigarette smoke being the most important
environmental risk factor
[[\[1\]](#loc-9)]{style="text-decoration: underline"
role="doc-biblioref"}. Adjusting for age and sex, smokers are 3.18 times
more likely to develop COPD than non-smokers
[[\[3\]](#loc-11)]{style="text-decoration: underline"
role="doc-biblioref"} yet fewer than 50% of smokers develop the disease
[[\[4\]](#loc-12)]{style="text-decoration: underline"
role="doc-biblioref"}. Males are also 1.48 times more likely to develop
the disease compared to females, and the odds of developing COPD
increase with age [[\[3\]](#loc-11)]{style="text-decoration: underline"
role="doc-biblioref"}. COPD's incomplete penetrance suggests that
downstream responses to smoking, rather than exposure alone, may mediate
progression to clinical COPD.

Although the lung is the main organ affected, the clinical burden of
COPD extends beyond the respiratory system. Comorbidities such as muscle
wasting, cardiovascular disease, osteoporosis, and depression are
consistent with generalized metabolic disturbance
[[\[5\]](#loc-13)]{style="text-decoration: underline"
role="doc-biblioref"}. Previous studies have found hundreds of
metabolites associated with lung function measures like Forced
Expiratory Volume (FEV) and Forced Vital Capacity (FVC), and
multi-metabolite signatures can predict COPD severity
[[\[5\]](#loc-13)]{style="text-decoration: underline"
role="doc-biblioref"},
[[\[6\]](#loc-14)]{style="text-decoration: underline"
role="doc-biblioref"}. Nicotine has a short half-life so metabolomic
profiles of current smokers typically reflect its more stable
metabolite, cotinine (and related derivatives), which serves as a
reliable biomarker of recent tobacco exposure
[[\[7\]](#loc-15)]{style="text-decoration: underline"
role="doc-biblioref"}. Many analyses, however, miss disease-relevant
coordination among metabolites and focus on single metabolites. Smoking
likely reorganizes the inherently networked metabolite relationships
through processes like inflammation and oxidative stress, revealing
metabolite relationships that may persist or rewire after smoking
cessation.

### Hypothesis to be tested

*I hypothesize that active smoking induces alterations in the structure
of metabolic interactions, and smoking cessation partially reverses some
relationships.*

- Smoking purturbs metabolic programs through inflammatory and oxidative
  processes, producing condition-specific metabolite-metabolite
  relationships.

- Active smokers may have a distinct metabolic response that may not be
  fully reversed by smoking cessation.

### Technical approach

I will model metabolomic regulation by estimating a sparse network,
where edges represent unique metabolite associations by using
Condition-Adaptive Fused Graphical Lasso, and I will identify
conditional dependencies and condition-specific network edges that
reflect differences between cotinine-positive and cotinine-negative
metabolite networks.

### Evaluation approach

- I will evaluate edge stability by bootstrapping network estimation and
  report results for highly connected nodes as confidence intervals on
  the edge weights. I will then calculate Jaccard similarity between the
  edges of replicate networks.

- I will validate my results by assessing how metabolite edges reflect
  known metabolic pathways. For example, I will expect to see a
  subnetwork/s related to nicotine metabolism, innflamation, and
  oxidative stress.

### Specific Aims / Training goals

**Specific Aim 1:** *Construct metabolite networks stratified by smoking
status*

I will estimate sparse metabolite networks for current and former
smokers using CFGL by modeling both shared and condition-specific
(current/former smoking) metabolite associations. By comparing these
networks within controls with normal lung function, this aim will
identify which metabolic interactions represent baseline function versus
smoking-associated perturbations.

Results:

1.  Precision matrices for each condition (current/former smoking)

2.  Sets of shared and condition-specific edges

**Specific Aim 2:** *Interpret network topology to identify biologically
meaningful smoking-associated metabolic programs*

I will characterize smoking-specific networks by identifying metabolite
hubs and use Metabolon annotations to find pathway-level enrichment. By
testing whether known smoking related metabolites in nicotine metabolism
associate in plausible ways, this aim evaluates whether network
differences reflect coherent biology rather than chance artifacts.

Results:

1.  Metabolite hubs and pathway annotated sub-networks.

2.  Comparison of pathways across conditions.

3.  Coherence of networks with known smoking exposure biology.

**Specific Aim 3 (time permitting):** *Assess how COPD modifies
smoking-associated metabolic network structures*

I will extend the network framework to incorporate COPD status. This
will reveal whether COPD alters how smoking-associated metabolite
interactions persist, change magnitude, or reorganize. This aim will
test whether network structures change between smoking and COPD
conditions to see if smoking related changes are further perturbed by
COPD status.

Results:

1.  Differences in network associations for each condition:
    COPD+/current-smoking, COPD+/former-smoking, COPD-/current-smoking,
    and COPD-/former-smoking.

2.  Identification of metabolite relationships that persist between
    COPD+/- and between current/former smoking.

### Future Directions

Alternative approaches could include differential network analysis and
module based analysis (such as SmCCNet) to identify connected
subnetworks within the larger global networks.