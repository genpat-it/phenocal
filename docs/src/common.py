"""Shared pieces of the two GitHub Pages: style, head, citations and references.

Citations are written in the templates as @@CITE:key@@ (or @@CITE:key1,key2@@);
`cite()` numbers them in order of first appearance and `references()` renders
the list. Every DOI below was checked against Crossref.
"""
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))

REFS = {
    "phenocal": ('de Ruvo A, Castelli P, Di Pasquale A, Radomski N. PhenoCal: natural-control label '
                 'calibration for cross-cohort machine learning. In: <i>12th International Artificial '
                 'Intelligence Symposium (AIS 2026)</i>. Lecture Notes in Computer Science. Springer '
                 'Nature; 2026. Reproducibility archive: '
                 '<a href="https://doi.org/10.5281/zenodo.21461769">doi:10.5281/zenodo.21461769</a>.'),
    "chewbbaca": ('Silva M, Machado MP, Silva DN, Rossi M, Moran-Gilad J, Santos S, Ramirez M, Carriço JA. '
                  'chewBBACA: a complete suite for gene-by-gene schema creation and strain identification. '
                  '<i>Microbial Genomics</i>. 2018;4(3):e000166. '
                  '<a href="https://doi.org/10.1099/mgen.0.000166">doi:10.1099/mgen.0.000166</a>.'),
    "cgmlstdists": ('Seemann T. cgmlst-dists: calculate distance matrix from cgMLST allele call tables. '
                    'Software. <a href="https://github.com/tseemann/cgmlst-dists">'
                    'github.com/tseemann/cgmlst-dists</a>.'),
    "reportree": ('Mixão V, Pinto M, Sobral D, Di Pasquale A, Gomes J, Borges V. ReporTree: a '
                  'surveillance-oriented tool to strengthen the linkage between pathogen genetic clusters '
                  'and epidemiological data. <i>Genome Medicine</i>. 2023;15(1):43. '
                  '<a href="https://doi.org/10.1186/s13073-023-01196-1">doi:10.1186/s13073-023-01196-1</a>.'),
    "moura": ('Moura A, Criscuolo A, Pouseele H, et al. Whole genome-based population biology and '
              'epidemiological surveillance of <i>Listeria monocytogenes</i>. <i>Nature Microbiology</i>. '
              '2016;2:16185. <a href="https://doi.org/10.1038/nmicrobiol.2016.185">'
              'doi:10.1038/nmicrobiol.2016.185</a>.'),
    "efron": ('Efron B. Bootstrap methods: another look at the jackknife. <i>The Annals of Statistics</i>. '
              '1979;7(1):1–26. <a href="https://doi.org/10.1214/aos/1176344552">'
              'doi:10.1214/aos/1176344552</a>.'),
    "silverman": ('Silverman BW. <i>Density Estimation for Statistics and Data Analysis</i>. London: '
                  'Chapman &amp; Hall; 1986.'),
    "combat": ('Johnson WE, Li C, Rabinovic A. Adjusting batch effects in microarray expression data using '
               'empirical Bayes methods. <i>Biostatistics</i>. 2007;8(1):118–127. '
               '<a href="https://doi.org/10.1093/biostatistics/kxj037">doi:10.1093/biostatistics/kxj037</a>.'),
    "leek": ('Leek JT, Scharpf RB, Bravo HC, et al. Tackling the widespread and critical impact of batch '
             'effects in high-throughput data. <i>Nature Reviews Genetics</i>. 2010;11(10):733–739. '
             '<a href="https://doi.org/10.1038/nrg2825">doi:10.1038/nrg2825</a>.'),
    "limma": ('Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK. limma powers differential '
              'expression analyses for RNA-sequencing and microarray studies. <i>Nucleic Acids Research</i>. '
              '2015;43(7):e47. <a href="https://doi.org/10.1093/nar/gkv007">doi:10.1093/nar/gkv007</a>.'),
    "frenay": ('Frénay B, Verleysen M. Classification in the presence of label noise: a survey. <i>IEEE '
               'Transactions on Neural Networks and Learning Systems</i>. 2014;25(5):845–869. '
               '<a href="https://doi.org/10.1109/TNNLS.2013.2292894">doi:10.1109/TNNLS.2013.2292894</a>.'),
    "shimodaira": ('Shimodaira H. Improving predictive inference under covariate shift by weighting the '
                   'log-likelihood function. <i>Journal of Statistical Planning and Inference</i>. '
                   '2000;90(2):227–244. <a href="https://doi.org/10.1016/S0378-3758(00)00115-4">'
                   'doi:10.1016/S0378-3758(00)00115-4</a>.'),
    "palma": ('Palma F, Radomski N, Guérin A, et al. Genomic elements located in the accessory repertoire '
              'drive the adaptation to biocides in <i>Listeria monocytogenes</i> strains from different '
              'ecological niches. <i>Food Microbiology</i>. 2022;106:103757. '
              '<a href="https://doi.org/10.1016/j.fm.2021.103757">doi:10.1016/j.fm.2021.103757</a>.'),
    "kragh": ('Kragh ML, Ivanova M, Dalgaard P, Leekitcharoenphon P, Truelstrup Hansen L. Sensitivity of '
              '240 <i>Listeria monocytogenes</i> isolates to common industrial biocides is more dependent on '
              'the presence of residual organic matter or biofilm than on genetic determinants. '
              '<i>Food Control</i>. 2024;158:110244. <a href="https://doi.org/10.1016/j.foodcont.2023.110244">'
              'doi:10.1016/j.foodcont.2023.110244</a>.'),
    "ivanova": ('Ivanova M, Kragh ML, Szarvas J, et al. Large-scale phenotypic and genomic analysis of '
                '<i>Listeria monocytogenes</i> reveals diversity in the sensitivity to quaternary ammonium '
                'compounds but not to peracetic acid. <i>Applied and Environmental Microbiology</i>. '
                '2025;91:e01829-24. <a href="https://doi.org/10.1128/aem.01829-24">doi:10.1128/aem.01829-24</a>.'),
    "cooper": ('Cooper AL, Carrillo CD, Deschênes M, Blais BW. Genomic markers for quaternary ammonium '
               'compound resistance as a persistence indicator for <i>Listeria monocytogenes</i> '
               'contamination in food manufacturing environments. <i>Journal of Food Protection</i>. '
               '2021;84(3):389–398. <a href="https://doi.org/10.4315/JFP-20-328">doi:10.4315/JFP-20-328</a>.'),
    "he": ('He Y, Xu T, Li S, Mann DA, Britton B, Oliver HF, den Bakker HC, Deng X. Integrative assessment '
           'of reduced <i>Listeria monocytogenes</i> susceptibility to benzalkonium chloride in produce '
           'processing environments. <i>Applied and Environmental Microbiology</i>. 2022;88(21):e01269-22. '
           '<a href="https://doi.org/10.1128/aem.01269-22">doi:10.1128/aem.01269-22</a>.'),
}

AUTHORS = [  # name, affiliations, ORCID, equal contribution
    ("Andrea de Ruvo", "1,3", "0000-0001-6836-3455", True),
    ("Pierluigi Castelli", "2,3", "0000-0001-6518-1752", True),
    ("Adriano Di Pasquale", "3", "0000-0002-9328-3972", False),
    ("Nicolas Radomski", "3", "0000-0002-7480-4197", False),
]
AFFILIATIONS = [
    "Gran Sasso Science Institute (GSSI), L'Aquila, Italy",
    "Università degli Studi di Teramo (UNITE), Teramo, Italy",
    "Istituto Zooprofilattico Sperimentale dell'Abruzzo e del Molise “G. Caporale” (IZSAM), Teramo, Italy",
]

BIBTEX = """@inproceedings{deruvo2026phenocal,
  title     = {{PhenoCal}: Natural-Control Label Calibration for
               Cross-Cohort Machine Learning},
  author    = {de Ruvo, Andrea and Castelli, Pierluigi and
               Di Pasquale, Adriano and Radomski, Nicolas},
  booktitle = {12th International Artificial Intelligence Symposium (AIS 2026)},
  series    = {Lecture Notes in Computer Science},
  publisher = {Springer Nature},
  year      = {2026}
}"""


def authors_html():
    orcid = ('<a class="orcid" href="https://orcid.org/{0}" title="ORCID {0}" aria-label="ORCID {0}">'
             '<svg viewBox="0 0 16 16" width="13" height="13" aria-hidden="true"><circle cx="8" cy="8" r="8" '
             'fill="#a6ce39"/><text x="8" y="11.6" text-anchor="middle" font-size="9" font-weight="700" '
             'fill="#fff" font-family="system-ui,sans-serif">iD</text></svg></a>')
    names = ", ".join(f'<span class="au">{n}<sup>{a}{"*" if eq else ""}</sup>{orcid.format(o)}</span>'
                      for n, a, o, eq in AUTHORS)
    affs = "".join(f"<li><sup>{i + 1}</sup>{esc_text(a)}</li>" for i, a in enumerate(AFFILIATIONS))
    return (f'<p class="authors">{names}</p><ol class="affil">{affs}</ol>'
            '<p class="affil-note"><sup>*</sup>These authors contributed equally.</p>')


def esc_text(s):
    return s.replace("&", "&amp;")


def cite(html):
    """Number @@CITE:...@@ markers by first appearance; return html and the ordered keys."""
    order = []

    def repl(m):
        keys = [k.strip() for k in m.group(1).split(",")]
        nums = []
        for k in keys:
            if k not in REFS:
                raise SystemExit(f"unknown reference {k}")
            if k not in order:
                order.append(k)
            nums.append(order.index(k) + 1)
        links = ", ".join(f'<a href="#ref-{n}">{n}</a>' for n in sorted(nums))
        return f'<span class="cite">[{links}]</span>'

    html = re.sub(r"@@CITE:([^@]+)@@", repl, html)
    return html, order


def references(order):
    items = "".join(f'<li id="ref-{i + 1}">{REFS[k]}</li>' for i, k in enumerate(order))
    return f'<ol class="refs">{items}</ol>'


def head(title, description):
    css = open(os.path.join(HERE, "style.css")).read()
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<meta name="description" content="{description}">
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.min.css">
<script defer src="https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.min.js"></script>
<script defer src="https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/contrib/auto-render.min.js"
  onload="renderMathInElement(document.body,{{delimiters:[{{left:'$$',right:'$$',display:true}},{{left:'\\\\(',right:'\\\\)',display:false}}]}})"></script>
<style>
{css}</style>
</head>
"""


FOOTER = """<footer>
phenocal is free software under the MIT licence, developed at the Istituto Zooprofilattico
Sperimentale dell'Abruzzo e del Molise “G. Caporale” (IZSAM) and the Gran Sasso Science
Institute (GSSI). <a href="https://github.com/genpat-it/phenocal">github.com/genpat-it/phenocal</a>
</footer>
"""
