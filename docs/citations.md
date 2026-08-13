# VUR Model Citation Ledger

This file records where model default values came from and how they are used.

## Parameter Mapping

| ID | Parameter(s) Used | Value(s) Used in Code | Source |
|---|---|---|---|
| C1 | Ureter diameter trend by age | Infancy near `3.2 mm`; gradual increase to adolescent range | Shashi et al., Pediatric Radiology 2022. [PubMed](https://pubmed.ncbi.nlm.nih.gov/35386015/) |
| C2 | Ureter length by age | Approximation `ureter length (cm) = age + 12` | Zarzour et al., J Pediatr Urol 2019. [PubMed](https://pubmed.ncbi.nlm.nih.gov/31324475/) |
| C3 | Intravesical ureter / UVJ caliber and tunnel ratio context | For `1-3 years`, intravesical ureter diameter about `1.4 mm`; common anti-reflux ratio context near `5:1` | Meyer et al., translational review. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC3126077/) |
| C4 | Female urethral length in early childhood | `6-12 mo: 2.50 cm`, `12-24 mo: 2.31 cm`, `24-36 mo: 2.59 cm` | Sheldon et al., J Pediatr Surg 2019. [PubMed](https://pubmed.ncbi.nlm.nih.gov/30503195/) |
| C5 | Male urethral length age trend | Regression-based age trend used for male pediatric scaling (`~8.7 + 0.55*age years`, 1-15 y) | Dudhani et al., J Pediatr Urol 2023. [PubMed](https://pubmed.ncbi.nlm.nih.gov/37419832/) |
| C6 | Pediatric filling pressure target | End-filling detrusor pressure typically expected around `<=10 cmH2O` at expected capacity; detrusor overactivity commonly defined by rise `>15 cmH2O` | ICCS terminology/urodynamic standardization report. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC2245874/) |
| C7 | Infant voiding pressure high range | Maximum detrusor rise in infancy commonly around `95-120 cmH2O` | Holmdahl et al., Br J Urol 1988. [PubMed](https://pubmed.ncbi.nlm.nih.gov/2834578/) |
| C8 | High-pressure pattern in infant VUR | Max detrusor pressure higher in VUR cohort (`~161 cmH2O`) vs non-VUR infant cohort (`~117 cmH2O`) | Yeung et al., Br J Urol 1997. [PubMed](https://pubmed.ncbi.nlm.nih.gov/9258186/) |
| C9 | Sex-stratified normal voiding pressure ranges (pediatric) | Girls `~30-65 cmH2O`, boys `~55-80 cmH2O` used as non-VUR anchor range for older infants/children | Tekgül et al., pediatric neuro-urology review. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC5605603/) |
| C10 | Urethral caliber proxy for flow resistance | Infant/child meatal caliber review supports size in the `~10-12 Fr` range (used as urethral diameter proxy around `3.8 mm`) | Gadelkareem et al., Minerva Urol Nefrol 2021. [PubMed](https://pubmed.ncbi.nlm.nih.gov/34782235/) |
| C11 | VUR CFD phase/boundary context | Filling vs voiding concept, pressure-phase setup, and reflux under voiding pressure conditions | Kim et al., Comput Biol Med 2022. [DOI](https://doi.org/10.1016/j.compbiomed.2022.105456) |
| C12 | Ureter peristalsis and urine fluid assumptions | Urine density/viscosity assumptions and peristaltic-wave modeling context | Razavi & Jouybar, Bio-Med Mater Eng 2018. [DOI](https://doi.org/10.3233/BME-181026) |
| C13 | Full bladder capacity formula family (classic) | `EBC (mL) = (age + 2) x 30` for older children (Koff-style equation) | Austin et al., J Urol 2019 (capacity references in guideline review). [PubMed](https://pubmed.ncbi.nlm.nih.gov/31173172/) |
| C14 | Full bladder capacity formula family (Kaefer) | `<2y: (2 x age + 2) oz`; `>=2y: (age/2 + 6) oz` (converted to mL in code) | Kaefer et al., J Urol 1997. [PubMed](https://pubmed.ncbi.nlm.nih.gov/9102531/) |
| C15 | Full bladder capacity formula family (infant/weight) | Infant estimate `~7 x weight (kg)` mL, combined with age-based transitions | Holmdahl et al., J Urol 1996. [PubMed](https://pubmed.ncbi.nlm.nih.gov/8693709/) |
| C16 | International VUR grading morphology | Grade progression includes increasing ureteral dilation and tortuosity in higher grades (esp. IV-V) | Lebowitz et al., Radiology 1985. [PubMed](https://pubmed.ncbi.nlm.nih.gov/3909456/) |
| C17 | Ureter wall compliance / stiffness context | Fibrotic/scarred ureter mechanics and altered collagen-elastic response linked to reduced compliance | Knudsen et al., Neurourol Urodyn 1994. [PubMed](https://pubmed.ncbi.nlm.nih.gov/7887788/) |
| C18 | Peristalsis propagation velocity context | Average ureteral peristaltic wave speed around `~2 cm/s` used for baseline wave dynamics | Griffiths et al. (as summarized in Razavi & Jouybar references). [DOI](https://doi.org/10.3233/BME-181026) |
| C19 | Renal pelvic pressure as obstruction-severity signal in pediatric CFD workflows | Supports using pressure-centric outcome tracks (`peak_renal_pelvis_pressure_pa`, obstruction index trend) for scenario ranking; no direct fixed constant imported | Nishimura et al., Front Urol 2025. [DOI](https://doi.org/10.3389/fruro.2025.1634278) |
| C20 | Pediatric upper-tract CFD feature context | Supports pressure-flow morphology outputs for post-reconstruction comparison and model extension planning; no direct fixed constant imported | Yang et al., Comput Methods Programs Biomed 2026. [DOI](https://doi.org/10.1016/j.cmpb.2025.109077) |
| C21 | Bladder wall deformation effect on UVJ closure/resistance | FE modeling shows UVJ closure tracks bladder wall deformation and changing UVJ geometry with storage, not pressure alone; reference state near 10% capacity | Kalayeh et al., Neurourol Urodyn 2020. [PubMed](https://pubmed.ncbi.nlm.nih.gov/33017072/) |
| C22 | BBD impact on rUTI risk in VUR | In toilet-trained children, combined BBD + VUR has markedly higher recurrent UTI risk than either alone | Shaikh et al., Pediatrics 2016 (RIVUR/CUTIE analysis). [PubMed](https://pubmed.ncbi.nlm.nih.gov/26647376/) |
| C23 | BBD prevalence and recurrent UTI effect size in primary VUR | Meta-analysis supports BBD as common in VUR and associated with higher recurrent UTI risk | Madaan et al., Front Pediatr 2020. [DOI](https://doi.org/10.3389/fped.2020.00084) |
| C24 | Clinical impact of BBD on reflux resolution and treatment outcomes | Guideline synthesis notes lower spontaneous resolution and lower endoscopic success when BBD is present, motivating explicit BBD modeling in risk stratification | AUA Vesicoureteral Reflux Guideline/Topic Summary. [AUA](https://www.auanet.org/guidelines-and-quality/guidelines/vesicoureteral-reflux-topics) |

## Notes On Use

- The model is a reduced-order simulator, not a direct VCUG image classifier.
- Where direct toddler (`18-24 months`) values were unavailable, interpolation/extrapolation from adjacent pediatric cohorts was used.
- These values are meant as priors for calibration; local institutional data should overwrite defaults when available.
