---
layout: inner
title: Publications
permalink: /publications/

publications:
  - doi: 10.1093/mnras/stz357
    inpress: true
    ads: false
    title: >-
      Testing the stability of supersonic ionized Bondi accretion flows 
      with radiation hydrodynamics
    authors:
      - name: Vandenbroucke B.
      - name: Sartorio N. S.
      - name: Wood K.
      - name: Lund K.
      - name: Falceta-Gon&ccedil;alves D.
      - name: Haworth T. J.
      - name: Bonnell I.
      - name: Keto E.
      - name: Tootill D.
    journal: Monthly Notices of the Royal Astronomical Society
    volume: false
    page: false
    month: 2
    year: 2019
    description: >-
      Paper about the numerical aspects of spherically symmetric Bondi
      accretion with radiation. At the end of the paper, the effect of 1D
      stability issues on 3D RHD simulations is illustrated using CMacIonize.
  - doi: 10.1093/mnras/sty554
    ads: 2018MNRAS.476.4032V
    title: >-
      Radiative transfer calculations of the diffuse ionized gas in disc 
      galaxies with cosmic ray feedback
    authors:
      - name: Vandenbroucke B.
      - name: Wood K.
      - name: Girichidis P.
      - name: Hill A. S.
      - name: Peters T.
    journal: Monthly Notices of the Royal Astronomical Society
    volume: 476
    page: 4032
    month: 5
    year: 2018
    description: >-
      Study of the diffuse ionized gas in disc galaxies that used 
      CMacIonize to post-process snapshots from the SILCC project.
  - doi: 10.1016/j.ascom.2018.02.005
    ads: 2018A%26C....23...40V
    title: >-
      The Monte Carlo photoionization and moving-mesh radiation 
      hydrodynamics code CMacIonize
    authors:
      - name: Vandenbroucke B.
      - name: Wood K.
    journal: Astronomy and Computing
    volume: 23
    page: 40
    month: 4
    year: 2018
    description: CMacIonize code paper.
---

This is a list of publications about CMacIonize, or making use of 
CMacIonize. If you have used CMacIonize in your work and want to add it 
to this list, please contact us!

{% assign sorted_pubs = page.publications | sort:"month" %}
{% assign sorted_pubs = sorted_pubs | sort:"year" %}

{% for publication in sorted_pubs reversed %}
{% if publication.doi %}[{{ publication.title }}](http://dx.doi.org/{{ publication.doi }}){% elsif publication.url %}[{{ publication.title }}]({{ publication.url }}){% else %}{{ publication.title }}{% endif %}{% if publication.ads %} ([ADS record](http://adsabs.harvard.edu/abs/{{ publication.ads }})){% endif %}
: *{% for name in publication.authors %}{{ name.name }}, {% endfor %}{{ publication.year }},
    {% if publication.journal %}{{ publication.journal }}, {% if publication.inpress %}in press{% else %}{{ publication.volume }}, {{ publication.page }}{% endif %}*  {% elsif publication.bookinfo %}{{ publication.bookinfo }}*  {% endif %}
    {{ publication.description }}
{% endfor %}
