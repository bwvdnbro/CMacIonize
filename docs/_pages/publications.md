---
layout: inner
title: Publications
permalink: /publications/

publications:
  - doi: 10.1093/mnras/sty554
    title: Radiative transfer calculations of the diffuse ionized gas in disc galaxies with cosmic ray feedback
    authors:
      - name: Vandenbroucke, B.
      - name: Wood, K.
      - name: Girichidis, P.
      - name: Hill, A. S.
      - name: Peters, T.
    journal: Monthly Notices of the Royal Astronomical Society
    volume: 476
    page: 4032
    month: 5
    year: 2018
    description: Study of the diffuse ionized gas in disc galaxies that used CMacIonize to post-process snapshots from the SILCC project.
  - doi: 10.1016/j.ascom.2018.02.005
    title: The Monte Carlo photoionization and moving-mesh radiation hydrodynamics code CMacIonize
    authors:
      - name: Vandenbroucke, B.
      - name: Wood, K.
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

{% assign sorted_pubs = (page.publications | sort:"month") %}
{% assign sorted_pubs = (sorted_pubs | sort:"year") %}

{% for publication in sorted_pubs reversed %}
{% if publication.doi %}
[{{ publication.title }}](http://dx.doi.org/{{ publication.doi }})
{% elsif publication.url %}
[{{ publication.title }}]({{ publication.url }})
{% else %}
{{ publication.title }}
{% endif %}
: *{% for name in publication.authors %}{{ name.name }}, {% endfor %}{{ publication.year }},
    {% if publication.journal %}{{ publication.journal }}, {{ publication.volume }}, {{ publication.page }}*  {% elsif publication.bookinfo %}{{ publication.bookinfo }}*  {% endif %}
    {{ publication.description }}
{% endfor %}
