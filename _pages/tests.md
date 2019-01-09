---
layout: inner
title: Test problems
permalink: /tests/

tests:
  - title: Kelvin-Helmholtz test with and without bulk velocity
    path: kh_bulk
    codes:
      - name: MPI-AMRVAC
        path: amrvac
        files:
          - name: README.txt
            description: Short description of how to run the test
          - name: amrvacusr.t
            description: Setup code for the test
          - name: kh_fixed_fd.par
            description: Parameter file for 100x100 finite difference test without bulk velocity
          - name: kh_fixed_hllc.par
            description: Parameter file for 100x100 finite volume test without bulk velocity
          - name: kh_fixed_hllc_hires.par
            description: Parameter file for 400x400 finite volume test without bulk velocity
          - name: kh_v100_fd.par
            description: Parameter file for 100x100 finite difference test with bulk velocity v = 100
          - name: kh_v100_hllc.par
            description: Parameter file for 100x100 finite volume test with bulk velocity v = 100
          - name: kh_v100_hllc_hires.par
            description: Parameter file for 400x400 finite volume test with bulk velocity v = 100
      - name: Shadowfax
        path: shadowfax
        files:
          - name: README.txt
            description: Short description of how to run the test
          - name: kh_fixed.ini
            description: Parameter file for 100x100 setup without bulk velocity
          - name: kh_fixed_hires.ini
            description: Parameter file for 400x400 setup without bulk velocity
          - name: kh_v100.ini
            description: Parameter file for 100x100 setup with bulk velocity v = 100
          - name: kh_v100_hires.ini
            description: Parameter file for 400x400 setup with bulk velocity v = 100
          - name: make_IC_fixed.py
            description: Initial condition generating script for 100x100 setup without bulk velocity
          - name: make_IC_fixed_hires.py
            description: Initial condition generating script for 400x400 setup without bulk velocity
          - name: make_IC_v100.py
            description: Initial condition generating script for 100x100 setup with bulk velocity v = 100
          - name: make_IC_v100_hires.py
            description: Initial condition generating script for 400x400 setup with bulk velocity v = 100
---

This page contains a number of initial condition and parameter files that can be
used to run test problems that are somewhat more involved than the test problems
included in the default testsuite that is part of Shadowfax.

We provide initial conditions for both Shadowfax and public competitor codes, so
that anyone can repeat the tests we carried out with the code. For a detailed
description of the results and comparison between Shadowfax and other methods,
please refer to [the relevant publications][publications_link].

We have compared Shadowfax with MPI-AMRVAC[^amrvac] and SWIFT[^swift].

Overview of the test problems:
{% for test in page.tests %}
- [{{ test.title }}](#{{ test.path }})
{% endfor %}

{% for test in page.tests %}
## {{ test.title }} {#{{ test.path }}}

{% for code in test.codes %}
### {{ code.name }}

<div class="row">

{% for file in code.files %}
<div class="col-md-3">
  <div class="panel panel-primary">
    <div class="panel-heading">
     <a class="list-group-item active" href="{{ site.data.global.url }}/tests/{{ test.path }}/{{ code.path }}/{{ file.name }}">
     <h3 class="panel-title">{{ file.name }}</h3>
     </a>
    </div>
    <div class="panel-body text-left">
      {{ file.description }}
    </div>
  </div>
</div>
{% endfor %}

</div>

{% endfor %}

{% endfor %}

[^amrvac]: [ascl:1208.014][amrvac_link]  
    Keppens, R., Meliani, Z., van Marle, A. J., Delmont, P., Vlasis, A.,
    van der Holst, B., 2012, Journal of Computational Physics, 231, 718

[^swift]: <https://gitlab.cosma.dur.ac.uk/swift/swiftsim>

[publications_link]: {{ site.data.global.url }}/publications
[amrvac_link]: http://ascl.net/1208.014
